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
  
  











