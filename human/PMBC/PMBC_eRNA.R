library(dplyr)
library(Seurat)
library(patchwork)
pbmc.data <- Read10X(data.dir = "~/eRNA/pbmcs/cell/reads/pos/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc10k", min.cells = 3, min.features = 80)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

plot<-VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf("~/eRNA/pbmcs/read_count_p.pdf")
plot1
plot2
plot
dev.off()



pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("~/eRNA/pbmcs/feature selection.pdf")
plot1
plot2
dev.off()

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
plot3<-ElbowPlot(pbmc)
pdf("~/eRNA/pbmcs/pca.pdf")

plot3
dev.off()
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
umap<-DimPlot(pbmc, reduction = "umap")
pdf("~/eRNA/pbmcs/umap-ADD.pdf")
umap

dev.off()


write.table(Idents(pbmc),"~/eRNA/pbmcs/cell/reads/pos//cluster.txt",quote=F)
