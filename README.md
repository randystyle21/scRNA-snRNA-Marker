### Welcome to the Github repositroy for 
Single Nucleus RNA Sequencing of Human Pancreatic Islets Identifies Novel Gene Sets and Distinguishes β-Cell Subpopulations with Dynamic Transcriptome Profiles. 
Here, I will 

### Before we begin the actual analysis of the data, load necessary packages. 
```
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(cowplot)
library(SingleCellExperiment)
library(ensembldb)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(presto)
library(raster)
library(rgdal)
library(classInt)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(SoupX)
library(AnnotationHub)
library(msigdb)
library(escape)
```
#### SoupX / Decontamination
# Before begin to make the actual analysis, begin with the ambient mRNA decontamination. 

```c
  toc <- Seurat::Read10X_h5('path for H5 file/filtered_feature_bc_matrix.h5')
  tod <- Seurat::Read10X_h5('path for H5 file/raw_feature_bc_matrix.h5')
  ```
# if you used GRCh38-mm10 genome reference due to the human cell graft to mice
Then you may want to consider omitting GRCh38_tag in the genes. All human genome will appear as capital letters, and mouse gene will be tagged with mm10, you don't really have to keep GRCh38 prefix.
If you want to keep it or used human/mouse only library then, you don't have to run this block. 
```c
  tod@Dimnames[[1]] -> alpha
  gsub(alpha,pattern='GRCh38_', replacement = '') -> alpha
  tod@Dimnames[[1]] <- alpha
  
  toc@Dimnames[[1]] -> alpha
  gsub(alpha,pattern='GRCh38_', replacement = '') -> alpha
  toc@Dimnames[[1]] <- alpha
```
## Now continue to create the ambient mRNA calibrated data object and write it to new matrix file
```
  sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
  sc <- estimateSoup(sc)
  
  srat = CreateSeuratObject(sc$toc)
  srat = NormalizeData(srat)
  srat = ScaleData(srat)
  srat = FindVariableFeatures(srat)
  srat = RunPCA(srat,pcs.compute=30)
  srat = RunTSNE(srat,dims.use=seq(30))
  srat <- FindNeighbors(srat, dims=1:20)
  srat = FindClusters(srat,dims.use=seq(30),resolution=1)
  metadata = as.data.frame(srat@reductions$tsne@cell.embeddings)
  colnames(metadata) = c('RD1','RD2')
  
  metadata$Cluster = srat@meta.data$seurat_clusters
  metadata$Annotation <- srat@meta.data$seurat_clusters
  
  toc = sc$toc
  scNoDrops = SoupChannel(toc, toc, calcSoupProfile = FALSE)
  # Calculate soup profile
  soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
  scNoDrops = setSoupProfile(scNoDrops, soupProf)
  
  ##add metadata for visuality checks
  sc = setClusters(sc, setNames(metadata$Cluster, rownames(metadata)))
  sc = setDR(sc, metadata[colnames(sc$toc), c("RD1", "RD2")])
  ```
  
  ## Assign the contamination fraction. 
  We have used 20%, yet you can adjust according to your sample condition
 ```c 
  sc = setContaminationFraction(sc, 0.2)
  out = adjustCounts(sc)
  DropletUtils:::write10xCounts('.desired_filename', out) 
 ```
  if you simply want to use autoestimate, you can use
  ```c
  sc = autoEstCont(sc)
  ```
  for more detailed algorithm and tutorial, please visit SoupX's github <https://github.com/constantAmateur/SoupX> 
  or demo tutorial written by Matthew Daniel Young at <https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html>.
  
  ## Creation of the Seurat Data Object
 ```c 
  in_data_dir <- ("path to the directories for Soupxed Matrix")
  samples <- dir(in_data_dir)
  
   
  
  seurat_list <- lapply(samples, function(sample){
    cur_data <- Read10X(paste0(in_data_dir,"/",sample))
    cur_seurat <- CreateSeuratObject(
      counts = cur_data,
      min.cells=3,
      min.features=200,
    )
    cur_seurat$SampleID <- sample
    return(cur_seurat)
  })
  
  merged_seurat <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])

  ```
  
  # Filtering 
  Filtering hyperparameter is a tricky issue, especially when you deal with the GRCh38-mm10. All other parameters are followed commonly used industry standard. 
  Yet, the mouse ratio can be tailored made for your purposes of research and how tolerant you would be for mouse gene contamination. We have used 10%. 
  
  ```c
   
  merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
  merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
  merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
  merged_seurat@meta.data -> metadata
  
  merged_seurat$mouseRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "mm10--")
  
  
  filtered_seurat <- subset(x = merged_seurat, 
                            subset= (nCount_RNA >= 500) & 
                              (nFeature_RNA >= 250) & 
                              (log10GenesPerUMI > 0.80) & 
                              (mitoRatio < 0.20) & 
                              (mouseRatio < 10))
  
  counts <- GetAssayData(object = filtered_seurat, slot = "counts")
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  filtered_counts <- counts[keep_genes, ]
  filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
  
  ```
  ## Doublet Removal 
  I believe Doubletfinder is a very convenient tool to control the multiplet in scRNA data. You can check out more detailed code from <https://github.com/chris-mcginnis-ucsf/DoubletFinder>.
  Yet, I will just share replica of their tutorial for convenience. 
  
  ```c
  library(ROCR)
  library(Splatter)
  library(Matrix)
  library(KernSmooth)
  library(ROCR)
  library(parallel)
  library(AnnotationHub)
  library(ensembldb)
  library(DoubletFinder)
  
  
  ###pre-process
  filtered_seurat <- SCTransform(filtered_seurat)
  filtered_seurat <- RunPCA(filtered_seurat)
  filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:10)
  
  
  ## pK Identification (no ground-truth) ------------------------------------------------------------------------------------------
  sweep.res.list_filtered_seurat <- paramSweep_v3(filtered_seurat, PCs = 1:10, sct = TRUE)
  sweep.stats_filtered_seurat <- summarizeSweep(sweep.res.list_filtered_seurat, GT = FALSE)
  pk_obj <- find.pK(sweep.stats_filtered_seurat)
  
  homotypic.prop <- modelHomotypic(filtered_seurat@meta.data$seurat_clusters)           ## ex: annotations <- filtered_seurat@meta.data$ClusteringResults
  ```
  
  #Parameter set up can be done from here. You can adjust the formation estimation according to your data set.
  ```c
  nExp_poi <- round(0.20*nrow(filtered_seurat@meta.data))  ## Assuming 20% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  filtered_seurat <- doubletFinder_v3(filtered_seurat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
  filtered_seurat@meta.data -> metaclean
  ```
  # Change the column name for the Doubletfinder for convenience. 
  your_number is unique to your data set, so change it with your number accordingly. Also, it's a good point that if you want to visualize the Doubletfinder result here before we remove them. 
  ```c
  names(metaclean)[names(metaclean) == 'DF.classifications_0.25_0.09_your_number'] <- 'Doubletfinder'
  filtered_seurat@meta.data <- metaclean
  DimPlot(filtered_seurat, group.by="Doubletfinder")
  
  filtered_seurat <- subset(filtered_seurat, Doubletfinder=="Singlet")
  ```
  
  # Integration 
  So, people performs phase control by cell cycle scoring. For detialed coding and logic behind of this, please refer to Seurat's original vignette <https://satijalab.org/seurat/articles/cell_cycle_vignette.html>. 
  You can also download the file that you need for scoring from the link.  
  
  
  ```c
    #cell phase scoring
  seurat_phase <- NormalizeData(filtered_seurat)
  # Load cell cycle markers
  load("your_cellcycle_data_file.rda")
  
  # Score cells for cell cycle
  seurat_phase <- CellCycleScoring(seurat_phase, 
                                   g2m.features = g2m_genes, 
                                   s.features = s_genes)
  # Identify the most variable genes
  seurat_phase <- FindVariableFeatures(seurat_phase, 
                                       selection.method = "vst",
                                       nfeatures = 2000, 
                                       verbose = FALSE)
  
  # Scale the counts
  seurat_phase <- ScaleData(seurat_phase)
  # Perform PCA
  seurat_phase <- RunPCA(seurat_phase)
  
  # Plot the PCA colored by cell cycle phase
  DimPlot(seurat_phase,
          reduction = "pca", group.by="Phase")
  options(future.globals.maxSize = 256000 * 1024^2)
  split_seurat <- SplitObject(filtered_seurat, split.by = "SampleID")
  
  split_seurat <- split_seurat[c("your data object names")]
  
  for (i in 1:length(split_seurat)) {
    split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
    split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
    split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
  }
  
  #Select Anchors & Integration Features 
  integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                              nfeatures = 3000) 
  split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                     anchor.features = integ_features)
  integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                          normalization.method = "SCT", 
                                          anchor.features = integ_features)
  seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
  
  
  ```
  
## Clustering
So, now we can create the cluster with different Louvain resolution. You can simply assign lower on for convenience, or if you are looking for very small subpopultion, then you may go higher.
```c
  seurat_integrated <- RunPCA(object = seurat_integrated)
  PCAPlot(seurat_integrated, split.by = "seq_type")  
  seurat_integrated <- RunUMAP(seurat_integrated, 
                               dims = 1:40,
                               reduction = "pca", return.model=TRUE)
  
  
  DimPlot(seurat_integrated)          
  
  
  seurat_integrated <- FindNeighbors(object = seurat_integrated, dims = 1:40)
  seurat_integrated <- FindClusters(object = seurat_integrated, resolution = c(0.2,0.8, 1.4,1.8,2.0))
  seurat_integrated@meta.data %>% 
    View()
  Idents(object = seurat_integrated) <- "integrated_snn_res.your_desired_resolution"
  
  
## Azimuth Reference Clustering & Label 
Azimuth is an excellent tool to project your data to publically available reference, especially when you created the data with Seurat. Yet, the online based app didn't allow to upload large sized integrated data at the moment of this research. Instead, I have downloaded the analysis script template which is provided when you do the analysis with the app. You can download it with the same way from <https://app.azimuth.hubmapconsortium.org/app/human-pancreas>
```c

  ############# Azimuth
  library(Seurat)
  library(glmGamPoi)
  library(Azimuth)
  # Load helper functions from Azimuth
  
  source("https://raw.githubusercontent.com/satijalab/azimuth/master/R/helpers.R")
  
  # Download the Azimuth reference and extract the archive
  
  # Load the reference
  ```
  
  This loading reference function may or maynot work with your own platform. It functions well with Azimuth app, but this link doesn't work time to time. 
  ```c
  # Change the file path based on where the reference is located on your system.
  reference <- Azimuth:::LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pancreas")

   
  # Preprocess with SCTransform
  seurat_integrated <- SCTransform(
    object = seurat_integrated,
    assay = "RNA",
    new.assay.name = "refAssay",
    residual.features = rownames(x = reference$map),
    reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
    method = 'glmGamPoi',
    n_cells=2000,
    n_genes = 2000,
    do.correct.umi = FALSE,
    do.scale = FALSE,
    do.center = TRUE
  )
  
  
  
  features = intersect(rownames(x = reference$map), VariableFeatures(object = seurat_integrated))
  
  # Find anchors between query and reference
  anchors <- FindTransferAnchors(
    reference = reference$map,
    query = seurat_integrated,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = "refAssay",
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = features,
    dims = 1:50,
    n.trees = 20,
    mapping.score.k = 100
  )
  
  
  refdata <- lapply(X = "annotation.l1", function(x) {
    reference$map[[x, drop = TRUE]]
  })
  names(x = refdata) <- "annotation.l1"
  if (FALSE) {
    refdata[["impADT"]] <- GetAssayData(
      object = reference$map[['ADT']],
      slot = 'data'
    )
  }
  
  
  seurat_integrated <- TransferData(
    reference = reference$map,
    query = seurat_integrated,
    dims = 1:50,
    anchorset = anchors,
    refdata = refdata,
    n.trees = 20,
    store.weights = TRUE
  )
  
  # Calculate the embeddings of the query data on the reference SPCA
  seurat_integrated <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = reference$map,
    query = seurat_integrated,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE
  )
  
  # Calculate the query neighbors in the reference
  # with respect to the vitro embeddings
  seurat_integrated[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(reference$map[["refDR"]]),
    query = Embeddings(seurat_integrated[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  seurat_integrated <- Azimuth:::NNTransform(
    object = seurat_integrated,
    meta.data = reference$map[[]]
  )
  
  # Project the query to the reference UMAP.
  seurat_integrated[["proj.umap"]] <- RunUMAP(
    object = seurat_integrated[["query_ref.nn"]],
    reduction.model = reference$map[["refUMAP"]],
    reduction.key = 'UMAP_'
  )
  
  
  # Calculate mapping score and add to metadata
  seurat_integrated <- AddMetaData(
    object = seurat_integrated,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  
  # First predicted metadata field, change to visualize other predicted metadata
  id <- "annotation.l1"[1]
  predicted.id <- paste0("predicted.", id)
  ```


# Apply the function to each element of the list
lapply(alpha_glucose_celltype, my_function)



## Primary Visualization 
I believe visualization is one of the most demending part for the journal publication, ironically. Anyhow, here is the sample for figures.
markers can be tailor made according to your data set yet I will share it here as well. Many markers are also obtained from Azimuth pancreatic islet reference, which contains datas from very prominent publications.  

```c
alphacell <- c("GCG","TTR","CRYBA2","SCG2")
betacell <- c("INS","IAPP","MAFA","DLK1")
deltacell <-c("SST","LY6H","LEPR","SYT1")
PPcell <- c("PPY","MEIS2","ID2","AQP3")
epsiloncell <-c("PHGR1","SPRN","UGT2B4","BHMT")
cyclingcell <- c("CENPF","TOP2A","MKI67","NUSAP1")
acinarcell <- c("REG1A","REG3A","PRSS1","CTRB2")
accell <- c("COL1A1","COL1A2","TIMP3","IGFBP5")
qscell <- c("RGS5","OLFML2A1","CSPG4","CCDC3")
ductalcell <- c("CFTR","SPP1","MMP7","KRT7")
immune <- c("CD3E","FCER1G","C1QB","CD68")
endothelialcell <- c("PLVAP","ESM1","RGCC","PECAM1")
schwanncell <- c("COL8A11","NGFR","CDH19","RUNX21","PLP1","PTPRZ1","DCT")
mastcell <- c("TPSAB1","CPA3","KIT")



allcellmarkers <- list("Alpha" = alphacell, "Beta" = betacell, "Delta" = deltacell, "PP" = PPcell, "Acinar" = acinarcell,
                       "Ductal"=ductalcell, "Endo"=endothelialcell, "TMM"=immune, "Mast"=mastcell, "Schwann"=schwanncell, "AS"=accell,"QS"=qscell,"Cycling"=cyclingcell)


nuclear_endomarkers <- list("Alpha"=c("PTPRT","FAP","PDK4","LOXL4"),
                            "Beta"=c("ZNF385D","TRPM3","LRFN2","PLUT"),
                            "Delta"=c("LRFN5","ADARB2","ERBB4","KCNT2"),
                            "Gamma"=c("CNTNAP5","CACNA2D3","THSD7A","RBFOX3"))


  levels(seurat_integrated) <- c("alpha","beta","delta","gamma",
                              "acinar","ductal","endothelial","schwann","immune",
                              "activated_stellate","quiescent_stellate")
  
  DotPlot(subset(seurat_integrated, idents=c("alpha","beta","delta","gamma")),scale.min=0, scale.max=100,
          col.min=-2.5, col.max=2.5, features=nuclear_endomarkers) + scale_colour_gradient2(low = "darkblue", mid = "white", high = "red") +  
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + scale_colour_gradient2(low = "darkblue", mid = "white", high = "red") + 
    theme(text = element_text(face = "bold"), 
          axis.text.x=element_text(angle=45, hjust=1, size=10), 
          axis.text.y=element_text(size=10)) + theme(panel.border=element_rect(color="darkgreen", fill=NA, size=2)) + 
    ggtitle("Whatever Data Names", subtitle="Annotation = Azimuth")
 ``` 
 

## Identify Cell Type Markers from DEG 
If you use Azimuth as your primary annotation resource, then cell type info is under meta.data$predicted.annotation.l1 in most of the cases, or you can change this any label that you have assigned. (Check out your seurat@met.data section.) If you run the 

```c

Idents(seurat_integrated) <- "predicted.annotation.l1"
seurat_celltype_list <- levels(seurat_integrated)

my_function <- function(celltype) {
  FindMarkers(subset(seurat_integrated, predicted.celltype == celltype), ident.1=celltype,
              ident.2=NULL) -> celltype_marker
  write.csv(celltype_marker, file=paste0("./",celltype,"_marker.csv"))
} 


lapply(seurat_celltype_list, my_function)
```

## Internal Referencing 
If you have your own reference like this study, you can simply do the reference based anntation & clustering as well 
```c

DefaultAssay(seurat_ingrated) <-"SCT"
DefaultAssay(seurat_vivo) <-"SCT"


anchors <- FindTransferAnchors(reference = seurat_ingrated, query=seurat_vivo, 
                               reference.reduction = "pca", normalization.method = "SCT",
                               reference.assay= "SCT", query.assay="SCT",recompute.residuals =FALSE, 
                               features=rownames(seurat_ingrated[["SCT"]])) 


seurat_vivo <- MapQuery(anchorset = anchors, reference = seurat_ingrated, 
                              query = seurat_vivo, 
                              refdata = list(celltype = "celltype"), reference.reduction = "pca", 
                              reduction.model = "umap") 
```
## Monocle 3 Pseudotime Analysis with Seurat Object 
```c

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)

gene_annotation <- as.data.frame(rownames(betaonly@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(betaonly@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"


cell_metadata <- as.data.frame(betaonly@assays[["RNA"]]@counts@Dimnames[[2]],
                               row.names = betaonly@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"


New_matrix <- betaonly@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(betaonly@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

library(monocle3)

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition



list_cluster <- betaonly@active.ident
names(list_cluster) <- betaonly@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster



cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-betaonly@reductions[["umap"]]@cell.embeddings


cds_from_seurat@preprocess_aux$gene_loadings <- betaonly@reductions[["pca"]]@feature.loadings

cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)

plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=4)
```

## Pathway Enrichment Test 
I believe the escape is an extremely useful package to calculate enrichment scores for each single cell. The beauty of this tool is that you can project the enrichment score in FeaturePlot so you can practically localize enrichment area as well. Anyhow, if you are interested, please refer to the github <https://github.com/ncborcherding/escape> or BioConductor Turorial page <http://bioconductor.org/packages/devel/bioc/vignettes/escape/inst/doc/vignette.html>. By the way, the whole pathway panel enrichment is extremely resource intensive and will take long. So you can select the pathways of interest as well. For our data, we ran the whole panel anyway but if you wish to do it, plan it wisely. We have done this analysis in beta cell subcluster only, so the name of the Seurat object is betaonly, and subcluster information is under betaonly@meta.data$subtype.


```c

library(escape)
library(dittoSeq)
library(SingleCellExperiment)
library(Seurat)
library(GSEABase)

# Getting Gene Set

GS.hallmark <- getGeneSets(library = "C")
GS.C2 <- getGeneSets(species="Homo sapiens",library="C2")
GS.C5 <- getGeneSets(species="Homo sapiens",library="C5")

ES.seurat <- enrichIt(obj = betaonly, gene.sets = GS.C2, groups = 1000, cores = 36)
ES.seurat <- enrichIt(obj = betaonly, gene.sets = GS.C5, groups = 1000, cores = 36)

# Add enrichment score to metadata section 
betaonly <- Seurat::AddMetaData(betaonly, ES.seurat)


colors <- colorRampPalette(c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20"))

#This may take  long
dittoHeatmap(betaonly, genes = NULL, metas = names('your_desired_pathway_name'), 
             annot.by = "subtype", 
             fontsize = 7, 
             cluster_cols = TRUE,
             heatmap.colors = colors(50))


ES2 <- data.frame(betaonly[[]], Idents(betaonly))
colnames(ES2)[ncol(ES2)] <- "cluster"

## plot
ridgeEnrichment(ES2, gene.set = "YOUR_PATHWAY_NAME", group = "sub1", add.rug = TRUE)


# The fitting can be done by different statistical model, refer to the tutorial link. 
output <- getSignificance(ES2, group = "cluster", fit = "linear.model")
```

## RNA Velocity - Loom File Creation  
First you need Loom file for each data, and you can get this by running velocyto algorithm in shell. For detailed tutorail for Velocyto, please refer to the <https://velocyto.org/velocyto.py/tutorial/cli.html>. 

```c
# Run this code in shell. 
velocyto run10x /your_data_folder /your_reference_folder/genes/genes.gtf -m /folder_contains_mask_file/grch38_rmsk.gtf
```

## RNA Velocity - scVelo
Sam Morivato shared excellent analysis pipeline, so check out his original script for doing this <https://smorabit.github.io/tutorials/8_velocyto/>. There is way to doing this from integrated Seurat with greater convenience, but I used following analysis pipeline for the publication . 

```c
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

sample_2c = anndata.read_loom("/media/randy/Expansion/H5/SCvsSN/B2C/velocyto/B2C.loom")
sample_2n = anndata.read_loom("/media/randy/Expansion/H5/SCvsSN/B2N/velocyto/B2N.loom")
sample_3c = anndata.read_loom("/media/randy/Expansion/H5/SCvsSN/B3C/velocyto/B3C.loom")
sample_3n = anndata.read_loom("/media/randy/Expansion/H5/SCvsSN/B3N/velocyto/B3N.loom")
sample_5c = anndata.read_loom("/media/randy/Expansion/H5/SCvsSN/B5C/velocyto/B5C.loom")
sample_5n = anndata.read_loom("/media/randy/Expansion/H5/SCvsSN/B5N/velocyto/B5N.loom")


sample_obs = pd.read_csv("your_cell_id_info.csv")
umap_cord = pd.read_csv("your_umap_coordinate_info.csv")
cell_clusters = pd.read_csv("your_cluster_naming_info.csv")
```

## Barcode Info & Loom 
set up the cell-id (barcode) info and adjust the string names for AnnData. 
```c
### B2C
sample_obs = sample_obs.rename(columns = {'x':'CellID'})


cellID_obs_B2C = sample_obs[sample_obs['CellID'].str.contains("-1_1")]


cellID_obs_B2C.CellID = "B2C:" + cellID_obs_B2C.CellID 
cellID_obs_B2C= cellID_obs_B2C["CellID"].str.replace("-1_1","x")

sample_2c = sample_2c[np.isin(sample_2c.obs.index,cellID_obs_B2C)]


### B3C

cellID_obs_B3C = sample_obs[sample_obs['CellID'].str.contains("-1_2")]


cellID_obs_B3C.CellID = "B3C:" + cellID_obs_B3C.CellID 
cellID_obs_B3C= cellID_obs_B3C["CellID"].str.replace("-1_2","x")

sample_3c = sample_3c[np.isin(sample_3c.obs.index,cellID_obs_B3C)]


### B5C


cellID_obs_B5C = sample_obs[sample_obs['CellID'].str.contains("-1_3")]


cellID_obs_B5C.CellID = "B5C:" + cellID_obs_B5C.CellID 
cellID_obs_B5C= cellID_obs_B5C["CellID"].str.replace("-1_3","x")

sample_5c = sample_5c[np.isin(sample_5c.obs.index,cellID_obs_B5C)]


### B2N 



cellID_obs_B2N = sample_obs[sample_obs['CellID'].str.contains("-1_4")]


cellID_obs_B2N.CellID = "B2N:" + cellID_obs_B2N.CellID 
cellID_obs_B2N= cellID_obs_B2N["CellID"].str.replace("-1_4","x")

sample_2n = sample_2n[np.isin(sample_2n.obs.index,cellID_obs_B2N)]

# B3N



cellID_obs_B3N = sample_obs[sample_obs['CellID'].str.contains("-1_5")]


cellID_obs_B3N.CellID = "B3N:" + cellID_obs_B3N.CellID 
cellID_obs_B3N= cellID_obs_B3N["CellID"].str.replace("-1_5","x")

sample_3n = sample_3n[np.isin(sample_3n.obs.index,cellID_obs_B3N)]


#B5N


cellID_obs_B5N = sample_obs[sample_obs['CellID'].str.contains("-1_6")]


cellID_obs_B5N.CellID = "B5N:" + cellID_obs_B5N.CellID 
cellID_obs_B5N= cellID_obs_B5N["CellID"].str.replace("-1_6","x")

sample_5n = sample_5n[np.isin(sample_5n.obs.index,cellID_obs_B5N)]



sample_2c.var_names_make_unique()
sample_3c.var_names_make_unique()
sample_5c.var_names_make_unique()

sample_2n.var_names_make_unique()
sample_3n.var_names_make_unique()
sample_5n.var_names_make_unique()

seurat_merged = sample_2c.concatenate(sample_3c, sample_5c,sample_2n,sample_3n, sample_5n)

seurat_cell = sample_2c.concatenate(sample_3c, sample_5c)
seurat_nuclear = sample_2n.concatenate(sample_3n, sample_5n)
```

## UMAP Coordinate & Cluster info
```c

umap_coord_B2C = umap_cord[umap_cord['CellID'].str.contains("-1_1")]
umap_coord_B3C = umap_cord[umap_cord['CellID'].str.contains("-1_2")]
umap_coord_B5C = umap_cord[umap_cord['CellID'].str.contains("-1_3")]
umap_coord_B2N = umap_cord[umap_cord['CellID'].str.contains("-1_4")]
umap_coord_B3N = umap_cord[umap_cord['CellID'].str.contains("-1_5")]
umap_coord_B5N = umap_cord[umap_cord['CellID'].str.contains("-1_6")]

### The Numbering disrupted by the python, so don't panic. #################################################

umap_coord_B2C.CellID = "B2C:" + umap_coord_B2C.CellID 
umap_coord_B2C.CellID = umap_coord_B2C["CellID"].str.replace("-1_1","x-0")


umap_coord_B3C.CellID = "B3C:" + umap_coord_B3C.CellID 
umap_coord_B3C.CellID = umap_coord_B3C["CellID"].str.replace("-1_2","x-1")


umap_coord_B5C.CellID = "B5C:" + umap_coord_B5C.CellID 
umap_coord_B5C.CellID = umap_coord_B5C["CellID"].str.replace("-1_3","x-2")





umap_coord_B2N.CellID = "B2N:" + umap_coord_B2N.CellID 
umap_coord_B2N.CellID = umap_coord_B2N["CellID"].str.replace("-1_4","x-3")


umap_coord_B3N.CellID = "B3N:" + umap_coord_B3N.CellID 
umap_coord_B3N.CellID = umap_coord_B3N["CellID"].str.replace("-1_5","x-4")


umap_coord_B5N.CellID = "B5N:" + umap_coord_B5N.CellID 
umap_coord_B5N.CellID = umap_coord_B5N["CellID"].str.replace("-1_6","x-5")

umap_merged = pd.concat([umap_coord_B2C,umap_coord_B2N, umap_coord_B3C, umap_coord_B3N, umap_coord_B5C, umap_coord_B5N])


umap_cell = pd.concat([umap_coord_B2C, umap_coord_B3C,  umap_coord_B5C])
umap_nuclear = pd.concat([umap_coord_B2N, umap_coord_B3N, umap_coord_B5N])




### Extract Index 

merged_seurat_index = pd.DataFrame(seurat_merged.obs.index)

seurat_cell_index = pd.DataFrame(seurat_cell.obs.index)
seurat_nuclear_index = pd.DataFrame(seurat_nuclear.obs.index)




# Merging UMAP data to the Sample 1 Index

umap_ordered = merged_seurat_index.merge(umap_merged, on = "CellID")



umap_cell_ordered = seurat_cell_index.merge(umap_cell, on="CellID")
umap_nuclear_ordered = seurat_nuclear_index.merge(umap_cell, on="CellID")



### rematch index for umap 

umap_nuclear.CellID = umap_nuclear["CellID"].str.replace("x-3","x-0")
umap_nuclear.CellID = umap_nuclear["CellID"].str.replace("x-4","x-1")
umap_nuclear.CellID = umap_nuclear["CellID"].str.replace("x-5","x-2")



umap_nuclear_ordered = seurat_nuclear_index.merge(umap_nuclear, on="CellID")




umap_ordered = umap_ordered.iloc[:,1:]

umap_cell_ordered = umap_cell_ordered.iloc[:,1:]
umap_nuclear_ordered = umap_nuclear_ordered.iloc[:,1:]




seurat_merged.obsm['X_umap'] = umap_ordered.values
seurat_cell.obsm['X_umap'] = umap_cell_ordered.values
seurat_nuclear.obsm['X_umap'] = umap_nuclear_ordered.values





#### Cell Type Annotation  

#### Renaming the cell ID Index


cell_clusters = cell_clusters.rename(columns = {'Unnamed: 0':'CellID'})


cluster_b2c = cell_clusters[cell_clusters['CellID'].str.contains("-1_1")]
cluster_b3c = cell_clusters[cell_clusters['CellID'].str.contains("-1_2")]
cluster_b5c = cell_clusters[cell_clusters['CellID'].str.contains("-1_3")]
cluster_b2n = cell_clusters[cell_clusters['CellID'].str.contains("-1_4")]
cluster_b3n = cell_clusters[cell_clusters['CellID'].str.contains("-1_5")]
cluster_b5n = cell_clusters[cell_clusters['CellID'].str.contains("-1_6")]


cluster_b2c.CellID = "B2C:" + cluster_b2c.CellID 
cluster_b2c.CellID = cluster_b2c["CellID"].str.replace("-1_1","x-0")

cluster_b3c.CellID = "B3C:" + cluster_b3c.CellID 
cluster_b3c.CellID = cluster_b3c["CellID"].str.replace("-1_2","x-1")

cluster_b5c.CellID = "B5C:" + cluster_b5c.CellID 
cluster_b5c.CellID = cluster_b5c["CellID"].str.replace("-1_3","x-2")

cluster_b2n.CellID = "B2N:" + cluster_b2n.CellID 
cluster_b2n.CellID = cluster_b2n["CellID"].str.replace("-1_4","x-3")


cluster_b3n.CellID = "B3N:" + cluster_b3n.CellID 
cluster_b3n.CellID = cluster_b3n["CellID"].str.replace("-1_5","x-4")

cluster_b5n.CellID = "B5N:" + cluster_b5n.CellID 
cluster_b5n.CellID = cluster_b5n["CellID"].str.replace("-1_6","x-5")



cluster_merged = pd.concat([cluster_b2c,cluster_b3c, cluster_b5c, cluster_b2n,cluster_b3n,cluster_b5n])

cluster_cell = pd.concat([cluster_b2c,cluster_b3c, cluster_b5c])
cluster_nuclear = pd.concat([cluster_b2n,cluster_b3n, cluster_b5n])




#### Rename Nuclear Clusters 

cluster_nuclear.CellID = cluster_nuclear["CellID"].str.replace("x-3","x-0")
cluster_nuclear.CellID = cluster_nuclear["CellID"].str.replace("x-4","x-1")
cluster_nuclear.CellID = cluster_nuclear["CellID"].str.replace("x-5","x-2")



## if you just merge by sample one index, you will pick up sample 1 only
cluster_merged = merged_seurat_index.merge(cluster_merged, on = "CellID")
cluster_cell_merged = seurat_cell_index.merge(cluster_cell, on = "CellID")
cluster_nuclear_merged = seurat_nuclear_index.merge(cluster_nuclear, on = "CellID")


cluster_ordered = cluster_merged.iloc[:,1]

cell_ordered = cluster_cell_merged.iloc[:,1]

nuclear_ordered = cluster_nuclear_merged.iloc[:,1]
```

## Now, add values to the clusters 
```c
seurat_merged.obs['cluster'] = cluster_ordered.values

seurat_cell.obs['cluster'] = cell_ordered.values
seurat_nuclear.obs['cluster'] = nuclear_ordered.values




seurat_merged.uns['clusters'] = cluster_ordered.values
seurat_merged.uns['cluster_colors'] = cluster_ordered.values





seurat_cell.uns['clusters'] = cluster_cell.values
seurat_cell.uns['cluster_colors'] = cluster_cell.values



seurat_nuclear.uns['clusters'] = cluster_nuclear.values
seurat_nuclear.uns['cluster_colors'] = cluster_nuclear.values

```
## Modeling & Visualization 
```c
scv.pp.filter_and_normalize(seurat_merged)
scv.pp.moments(seurat_merged)
scv.tl.velocity(seurat_merged, mode = "stochastic")

scv.tl.velocity_graph(seurat_merged)
scv.pl.velocity_embedding(seurat_merged, basis = 'umap')

#scv.pl.scatter(seurat_merged, color="cluster", palette={'beta_1': 'red', 'beta_2': 'blue', 'beta_3': 'green'})



scv.pl.proportions(seurat_merged, groupby="cluster")


scv.pl.proportions(seurat_cell, groupby="cluster") 
scv.pl.proportions(seurat_nuclear, groupby="cluster")


scv.pl.velocity_embedding_stream(seurat_merged, basis='umap', color='cluster')

scv.pl.velocity_embedding_stream(seurat_cell, basis='umap', color='cluster', title='scRNA')
scv.pl.velocity_embedding_stream(seurat_nuclear, basis='umap', color='cluster', title='snRNA')



scv.pl.velocity_embedding(seurat_merged, arrow_length=3, arrow_size=2, dpi=120, color='cluster')
```
