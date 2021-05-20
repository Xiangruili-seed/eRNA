library(dplyr)
library(Seurat)
library(patchwork)
args<-commandArgs(T)

dir<-args[1]
name<-args[2]
pbmc.data <- Read10X(data.dir = paste(dir,sep = ""))
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = paste(name,sep = ""), min.cells = 3, min.features = 300, names.field = 1, names.delim = "_", meta.data = NULL)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

plot<-VlnPlot(pbmc,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(dir,'read_count_p.pdf',sep = ""))
plot1
plot2
plot
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(paste(dir,'feature selection.pdf',sep = ""))
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
pdf(paste(dir,'pca.pdf',sep = ""))

plot3
dev.off()
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
umap<-DimPlot(pbmc, reduction = "umap")
pdf(paste(dir,'umap.pdf',sep = ""))
umap

dev.off()



pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




write.table(Idents(pbmc),paste(dir,'cluster.txt',sep = ""),quote=F)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
q10<-top10[which(top10$p_val_adj<0.05),]

top100 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top1000 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 1000, wt = avg_log2FC)
q<-top1000[which(top1000$p_val_adj<0.05),]
save(pbmc,file=paste(dir,'/.Rdata',sep = ""))

write.table(q,paste(dir,"top_0.05.txt",sep = ""),quote=F)
