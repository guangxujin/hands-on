### seurat hands-on

library(dplyr)
library(Seurat)
library(patchwork)
# set working directory
setwd("~/Downloads/Lecture_March 1st and 3rd/Install_packages/seurat")
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#save(pbmc.data,file='pbmc.data.RData')
load("pbmc.data.RData")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# cell seletion
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)




## Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

## Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
library(ggplot2)
library(ggrepel)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2




## Scaling the data

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

## Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

## Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot=ElbowPlot(pbmc)
ElbowPlot

## Cluster the cells

pbmc1 <- FindNeighbors(pbmc, dims = 1:10)
#pbmc1 <- FindClusters(pbmc1, resolution = 1, min=0.2)
pbmc1 <- FindClusters(pbmc1, resolution = 0.5)
head(Idents(pbmc), 5)

pbmc1 <- RunUMAP(pbmc1, dims = 1:10)
DimPlot(pbmc1, reduction = "umap")

pbmc1 <- RunTSNE(pbmc1, dims = 1:10)
DimPlot(pbmc1, reduction = "tsne")


## use the same parameters as tutorial
pbmc1 <- FindNeighbors(pbmc, dims = 1:10)
#pbmc1 <- FindClusters(pbmc1, resolution = 1, min=0.2)
pbmc1 <- FindClusters(pbmc1, resolution = 0.5)
head(Idents(pbmc), 5)

pbmc1 <- RunUMAP(pbmc1, dims = 1:10)
DimPlot(pbmc1, reduction = "umap")
pbmc<-pbmc1
##Finding differentially expressed features (cluster biomarkers)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)


library(dplyr)
library(tidyr)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

VlnPlot(pbmc, features = c("CCR7", "GZMK"))
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

## cell annotation
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) #+ NoLegend()

mkdir<-function(d){
  if (! file.exists(d)){
    dir.create(d)
  }
}
mkdir("../output")
saveRDS(pbmc, file = "../output/pbmc3k_final.rds")



