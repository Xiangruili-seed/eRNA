library(ggplot2)
args<-commandArgs(T)
annotation<-args[1]
type<-args[2]
dELS<- args[3]
GM12878_GROcap_over_one <- read.table(paste(annotation,'.txt',sep = ""), header = F)
GM12878_GROcap_below_one <- read.table(paste(annotation,'.',type,'.txt',sep = ""), header = F)

merged_score <- data.frame(cbind(rbind(GM12878_GROcap_over_one,GM12878_GROcap_below_one),c(rep('before_cut', 1000),rep('after_cut', 1000))))
colnames(merged_score) <- c('distance', 'Reads', 'group')

pdf(file = paste(annotation,'.',type,dELS,'.pdf',sep=""), width = 8, height = 8)
ggplot(data = merged_score, aes(x = distance, y = Reads, color = group)) + geom_line(size = 1) + theme_classic() + theme(legend.position = 'top') + guides(color = guide_legend(title = NULL, override.aes = list(size = 1))) + ggtitle(paste(annotation,sep = ""))+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
