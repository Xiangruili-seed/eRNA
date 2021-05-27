library(dplyr)
library(Seurat)
library(patchwork)
pbmc.data <- Read10X(data.dir = "~/eRNA/pbmcs/cell/reads/pos_add")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc10k", min.cells = 3, min.features = 300)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

plot<-VlnPlot(pbmc,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf("~/eRNA/pbmcs/pos/read_count_p.pdf")
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
pdf("~/eRNA/pbmcs/pos_add/feature selection.pdf")
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
pdf("~/eRNA/pbmcs/pos_add/pca.pdf")

plot3
dev.off()
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
umap<-DimPlot(pbmc, reduction = "umap")
pdf("~/eRNA/pbmcs/pos_add/umap.pdf")
umap

dev.off()






pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(Idents(pbmc),"~/eRNA/pbmcs/pos_add/cluster.txt",quote=F)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
q10<-top10[which(top10$p_val_adj<0.05),]
write.table(top10,"~/eRNA/pbmcs/pos_add/top_10.txt",quote=F)

top100 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top1000 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 1000, wt = avg_log2FC)
q<-top1000[which(top1000$p_val_adj<0.05),]

write.table(top100,"~/eRNA/pbmcs/pos_add/top_100.txt",quote=F)

write.table(q,"~/eRNA/pbmcs/pos_add/top_0.05.txt",quote=F)

load("/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/filtered_feature_bc_matrix/MAESTRO/data/human.immune.CIBERSORT.RData")

RNAAnnotateCelltypeCluster <- function(genes, signatures = "human.immune.CIBERSORT", cluster = 0){
    if(class(signatures) == "character"){
        data(list = signatures)
        signatures = get(signatures)
    }
    celltypes <- as.character(unique(signatures[,1]))
    signature_list <- sapply(1:length(celltypes),function(x){
                      return(toupper(as.character(signatures[which(signatures[,1]==celltypes[x]),2])))})
    names(signature_list) <- celltypes
    idx = genes$cluster==cluster
    avglogFC = genes$avg_log2FC[idx]
    names(avglogFC) = toupper(genes$gene[idx])
    score_cluster = sapply(signature_list, function(x){
                    score = sum(avglogFC[x], na.rm = TRUE) / log2(length(x))
                    return(score)
    })
    return(sort(score_cluster, decreasing=T)) 
}

celtype.score_0 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 0)
celtype.score_1 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 1)
celtype.score_2 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 2)
celtype.score_3 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 3)
celtype.score_4 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 4)
celtype.score_5 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 5)
celtype.score_6 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 6)
celtype.score_7 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 7)
celtype.score_8 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 8)
celtype.score_9 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 9)
celtype.score_10 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 10)
celtype.score_11 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 11)
celtype.score_12 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 12)
celtype.score_13 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 13)
celtype.score_14 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 14)
celtype.score_15 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 15)
celtype.score_16 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 16)
celtype.score_17 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 17)
celtype.score_18 <- RNAAnnotateCelltypeCluster(q, human.immune.CIBERSORT, cluster = 18)

names<-c(names(celtype.score_0[1]),names(celtype.score_1[1]),names(celtype.score_2[1]),names(celtype.score_3[1]),names(celtype.score_4[1]),names(celtype.score_5[1]),names(celtype.score_6[1]),names(celtype.score_7[1]),names(celtype.score_8[1]),names(celtype.score_9[1]),names(celtype.score_10[1]),names(celtype.score_11[1]),names(celtype.score_12[1]),names(celtype.score_13[1]),names(celtype.score_14[1]),names(celtype.score_15[1]),names(celtype.score_16[1]),names(celtype.score_17[1]),names(celtype.score_18[1]))

cluster<-cbind(names,c(0:18))
write.table(cluster,"~/eRNA/pbmcs/pos_add/cell_type.txt",quote=F)

new.cluster.ids <- cluster[,1]
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
plot4<-DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
pdf("~/eRNA/pbmcs/pos_add/umap selection.pdf")
plot4

dev.off()






erna<-read.table("~/eRNA/pbmcs/pos_add/top_0.05_erna.txt",header=F)

cell_typeA_marker_gene_list <-as.character(erna[,8])

 plot25<-FeaturePlot(pbmc, combine=F,keep.scale = "all",features =cell_typeA_marker_gene_list)


#~/eRNA/pbmcs/cell/reads/.RData


pdf("~/eRNA/pbmcs/pos_add/cell_typeA_score1.pdf")
plot25
dev.off()

erna<-read.table("~/eRNA/pbmcs/pos_add/top_10_erna.txt",header=F)

cell_typeA_marker_gene_list <-as.character(erna[,8])

 plot26<-FeaturePlot(pbmc, combine=F,keep.scale = "all",features =cell_typeA_marker_gene_list)

pdf("~/eRNA/pbmcs/pos_add/cell_typeA_score.pdf")
plot26
dev.off()


erna<-read.table("~/eRNA/pbmcs/pos_add/top_100_erna.txt",header=F)

cell_typeA_marker_gene_list <-as.character(erna[,8])

 plot27<-FeaturePlot(pbmc, combine=F,keep.scale = "all",features =cell_typeA_marker_gene_list)

pdf("~/eRNA/pbmcs/pos_add/cell_typeA_score100.pdf")
plot27
dev.off()
save(pbmc,file=paste("~/eRNA/pbmcs/pos_add//.Rdata",sep = ""))
