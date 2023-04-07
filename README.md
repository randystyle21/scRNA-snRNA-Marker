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
  nExp_poi <- round(0.15*nrow(filtered_seurat@meta.data))  ## Assuming 15% doublet formation rate - tailor for your dataset
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
  You can also download the file that you need for scoring from the link. Here, 
  
  
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

  







