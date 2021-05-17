library(ggplot2)
args<-commandArgs(T)
annotation<-args[1]
GM12878_GROcap_over_one <- read.table(paste('/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/fpos_plus_',annotation,'.txt',sep = ""), header = F)
GM12878_GROcap_below_one <- read.table(paste('/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/fpos_minus_',annotation,'.txt',sep = ""), header = F)
GM12878_GROcap <- read.table(paste('/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/fpos_',annotation,'.txt',sep = ""), header = F)

merged_score <- data.frame(cbind(rbind(GM12878_GROcap_over_one,GM12878_GROcap_below_one,GM12878_GROcap),c(rep('plus', 600),rep('minus', 600),rep('all', 600))))
colnames(merged_score) <- c('distance', 'Reads', 'group')

pdf(file = paste('/data/tusers/lixiangr/eRNA/single-cell/hg38/PBMCs/data/fpos_',annotation,'.pdf',sep=""), width = 8, height = 8)
ggplot(data = merged_score, aes(x = distance, y = Reads, color = group)) + geom_line(size = 1) + theme_classic() + theme(legend.position = 'top') + guides(color = guide_legend(title = NULL, override.aes = list(size = 1))) + ggtitle(paste('Reads for',annotation,' within cCREs-dELS',sep = ""))+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
