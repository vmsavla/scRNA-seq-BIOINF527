#title: "scRNA-seq_qc"
#author: "Varunika Savla"
#date: "2024-12-04"


# loading the required packages
library(Seurat)
library(harmony)
library(dplyr)
library(future)
library(magrittr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(data.table)
library(easybio)


# setting up the resources
num_cores <- 10
plan(multisession, workers = num_cores)
options(future.globals.maxSize=80000 *1024^2)


# set up working directory and local directories
# this is the directory that contains all the samples (preprocessed using CellRanger)
setwd('/nfs/turbo/dcmb-class/bioinf527/groups/group_02/varunika_files/prac_pipeline/')
parent_dir <- '/nfs/turbo/dcmb-class/bioinf527/groups/group_02/varunika_files/prac_pipeline/'
sample_dirs <- list.dirs(path = parent_dir, full.names = TRUE, recursive = FALSE)


## Function to perform Quality Control

# define the function perform_qc
perform_qc <- function (sample_dir) {
  # read in data
  data <- Read10X(data.dir = sample_dir)
  
  # create Seurat object
  data_Seurat_object <- CreateSeuratObject(data, project = basename(sample_dir))
  
  # calculate MT percentage
  data_Seurat_object[['percent_MT']] <- PercentageFeatureSet(data_Seurat_object, pattern = '^MT')
  
  # calculate the median number of genes
  median_genes <- median(data_Seurat_object$nFeature_RNA)
  
  # filter the cells
  data_Seurat_object <- subset(data_Seurat_object, subset = nFeature_RNA > 500 &
                                 nFeature_RNA < 2* median_genes & percent_MT < 10)
  
  # output file
  output_file <- file.path(sample_dir, paste0(basename(sample_dir),"_qc_results.rds"))
  saveRDS(data_Seurat_object, file = output_file)
}

# run quality control for each sample
for (sample_dir in sample_dirs) {
  tryCatch({
  sample_name <- basename(sample_dir)
  
  cat("Processing sample:", sample_name, "\n")
  
  perform_qc(sample_dir)
  
  cat("Processed sample:", sample_name, "\n")
  }, error = function(e) {
  cat("Error processing sample:", sample_name, "\n")
  cat("Error message:", conditionMessage(e), "\n")
  })
}


## Merging all Seurat objects

# merging seurat objects
seurat_files <- list.files(parent_dir, pattern= "_qc_results.rds", full.names = TRUE, recursive = TRUE)
seurat_objects <- lapply(seurat_files,readRDS)
basename(dirname(seurat_files))
merged_seurat <- merge(seurat_objects[[1]],
                       y = seurat_objects[-1],
                       add.cell.ids = basename(dirname(seurat_files)),
                       project = "MergedProject")

## Visualize merged data

# violin plot
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3, pt.size= 0)
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent_MT")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

## Further processing

# NormalizeData and FindVariableFeatures
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(merged_seurat), 10)
plot1 <- VariableFeaturePlot(merged_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# CellCycleScoring and ScaleData
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
merged_seurat <- JoinLayers(merged_seurat)
merged_seurat <- CellCycleScoring(merged_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
merged_seurat <- ScaleData(merged_seurat, vars.to.regress = c("S.Score", "G2M.Score"))

# RunPCA, RunHarmony, RunUMAP
merged_seurat <- RunPCA(merged_seurat, pc.genes = merged_seurat@var.genes, npcs = 30, verbose = FALSE)
options(repr.plot.height = 2.5, repr.plot.width = 6)
merged_seurat <- merged_seurat %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(merged_seurat, 'harmony')
harmony_embeddings[1:5, 1:5]
merged_seurat <- merged_seurat %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

## Visualize the UMAP

# UMAP by cluster and by sample
p1 <- DimPlot(merged_seurat, reduction = "umap", label = TRUE, pt.size = .5) + ggtitle("UMAP: Cells Colored by Clusters")
p2 <- DimPlot(mRCC, reduction = "umap", group.by = "orig.ident", label = FALSE) + ggtitle("UMAP: Cells Colored by Original Sample")
p1 + p2
ggsave("umap_clusters_and_samples_19.png", width = 16, height = 8)


## Identifying significant genes

# FindAllMarkers
merged_seurat.markers <- FindAllMarkers(merged_seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
merged_seurat.markers <- merged_seurat.markers[order(merged_seurat.markers$cluster, -merged_seurat.markers$avg_log2FC),]
write.csv(merged_seurat.markers,file="all_markers_clusters_19_ordered.csv", row.names = FALSE)

top_markers <- merged_seurat.markers %>%
  group_by(cluster) %>%
  top_n(50, wt = avg_log2FC)
write.csv(top_markers,file="top50_markers_clusters_19_ordered.csv", row.names = FALSE)


## Manual Cell Annotation

# bubble plot
# features1 is a list of markers used in the paper
features1 <- c("NDUFA4L2","CA9","KRT18","KRT8","KLRD1","KLRB1","GNLY","NKG7","CD3D","CD3E","CD8A","CD8B","IL7R","CD68","CD163","GPNMB","SLC40A1","MSR1","ACTA2","PDGFRB","COL1A2","PECAM1","KDR","CDH5","S100A8","S100A9","LYZ","KCNQ1OT1","CP","CTLA4","FOXP3","CD1C","CD1E","ACKR1","VWF","MKI67","TOP2A","CD79A","CD79B","IGKC","IGLC2","MS4A1","TPSB2","TPSAB1","KIT","KRT19","WFDC2")
p <- DotPlot(merged_seurat, features = features1, cols = c("blue","red"),dot.scale = 8)+RotatedAxis()
ggsave("dotplot_ref_markers.png", width = 16, height = 8)

# marker gene annotation - from cell marker reference used in the paper
features1_markers <- list(
  "ccRCC" = c("NDUFA4L2","CA9","KRT18","KRT8","KRT19","WFDC2"),
  "NKcells" = c("KLRD1","KLRB1","GNLY","NKG7"),
  "Tcells" =c("CD3D","CD3E","CD8A","CD8B","IL7R"),
  "Macrophages" = c("CD68","CD163","GPNMB","SLC40A1","MSR1"),
  "CAF" =c("ACTA2","PDGFRB","COL1A2"),
  "EndothelialCells"=c("PECAM1","KDR","CDH5","ACKR1","VWF"),
  "Monocytes" = c("S100A8","S100A9","LYZ"),
  "TregCells" = c("CTLA4","FOXP3"),
  "DendriticCells" = c("CD1C","CD1E"),
  "ProCD8+Tcells" = c("MKI67","TOP2A"),
  "Bcells" = c("CD79A","CD79B","MS4A1"),
  "PlasmaCells" = c("IGKC","IGLC2"),
  "MastCells" = c("TPSB2","TPSAB1","KIT")
  )

new.cluster.ids <- c(
  "0" = "ccRCC1",
  "1" = "CD8Tcells",
  "2" = "NKcells",
  "3" = "CD4Tcells",
  "4" = "Macrophages",
  "5" = "CAF",
  "6" = "EndoCells1",
  "7" = "Monocytes",
  "8" = "TregCells",
  "9" = "DendriticCells",
  "10" = "ProCD8+Tcells",
  "11" = "EndoCells2",
  "12" = "TAM",
  "13" = "ccRCC2",
  "14" = "ccRCC3",
  "15" = "PlasmaCells",
  "16" = "Bcells",
  "17" = "MastCells",
  "18" = "ccRCC4"
)

merged_seurat <- RenameIdents(merged_seurat, new.cluster.ids)
# UMAP with manual annotations
p3 <- DimPlot(merged_seurat, reduction = "umap", label = TRUE, pt.size = 0.5)+ ggtitle("UMAP: Cells Labeled by Clusters")
ggsave("umap_labeled_clusters_19.png", width = 8, height = 8)


## Pie chart for the proportion of cell types

# pie chart
b <- table(Idents(merged_seurat))
b<- sapply(levels(Idents(merged_seurat)),function(x){
  length(WhichCells(merged_seurat, idents = x))
})
sum(b)
sum_tcells <- sum(b[c(2,4,9,11)])
sum_macro <- sum(b[c(5,13)])
sum_endo <- sum(b[c(7,12)])
sum_dc <- sum(b[c(10)])
sum_mast <- sum(b[c(18)])
sum_tumor <- sum(b[c(1,14,15,19)])
sum_nk <- sum(b[c(3)])
sum_caf <- sum(b[c(6)])
sum_b <- sum(b[c(17,16)])
b_group <- c(sum_tcells,sum_macro,sum_endo,sum_dc,sum_mast,sum_tumor,sum_nk,sum_caf,sum_b)
column <- c("T cells", "Macrophages", "Endothelial cells", "Dendritic cells", "Mast cells", "Tumor cells", "NK cells", "CAF cells", "B cells")
row <- c("cell")
A <- rbind(b_group)
dimnames(A) = list(row,column)
lals_pct <- paste(column, "", "(",round(A/sum(A)*100,2), "%",")",sep="")
pie(A, labels = lals_pct, col = RColorBrewer::brewer.pal(n=9, name = "Set3"))
title("Proportion of cells in each category", cex.main=2)
dev.copy(png, "pie_prop_celltype.png", width = 800, height = 800)
dev.off()


## Visualization of manual cell annotation distribution of specific markers

# FeaturePlot for tumor cell genes
#FeaturePlot(merged_seurat, features = c("NDUFA4L2","CA9","KRT18","KRT8","KLRD1","CD3D","IL7R","CD68","MSR1","ACTA2","PECAM1","S100A8","KCNQ1OT1","CP","CTLA4","CD1C","MKI67","CD79A","IGLC2","TPSB2","KRT19","WFDC2"), ncol = 4)
FeaturePlot(merged_seurat, features = c("NDUFA4L2","CA9"),ncol=2)
ggsave("umap_feature_cancer.png", width = 16, height = 16)


## easybio to Perform Cell Type Annotation

# easybio
markerTop50Matched <- matchCellMarker2(marker = merged_seurat.markers, n = 50, spc = "Human")
easybio_top3_cluster <- markerTop50Matched[,head(.SD, 3), by = cluster][,1:4] |>knitr::kable()
write.csv(easybio_top3_cluster, file = "easybio_top3_celltype_per_cluster.csv", row.names = FALSE)
plotPossibleCell(markerTop50Matched[,head(.SD), by = cluster], min.uniqueN = 2)
cl2cell <- markerTop50Matched[,head(.SD, 1), by = cluster]
cl2cell<- setNames(cl2cell[["cell_name"]], cl2cell[["cluster"]])
merged_seurat@meta.data[["easy_annotation"]] <- cl2cell[as.character(Idents(merged_seurat))]
# UMAP with easybio annotation
color <- colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(19)
DimPlot(merged_seurat, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "easy_annotation")+ scale_color_manual(values = color)+ ggtitle("UMAP: Cells Labeled by Clusters")
