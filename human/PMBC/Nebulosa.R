#~/eRNA/pbmcs/pos/top_0.05_cor.txt

load("~/eRNA/pbmcs/cell/reads/.RData")
#A Genome-wide Framework for Mapping Gene Regulation via Cellular Genetic Screens
library("Nebulosa")
library("Seurat")
library("BiocFileCache")
data<-read.table("~/eRNA/pbmcs/pos/top_0.05_cor_cut.txt",header=F,stringsAsFactors = F)

pdf(file = paste('~/eRNA/pbmcs/pos/Nebulosa_cut.pdf',sep=""))

for (i in 1:nrow(data)) {
	plot<-plot_density(pbmc, as.character(data[i,1:2]), joint = TRUE,method = "wkde",pal = "viridis")
	print(plot+ plot_layout(ncol = 1))

}

dev.off()

