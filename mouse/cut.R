library(ggplot2)
library (plotrix)
library(grid)
library(ggrepel)
args<-commandArgs(T)
annotation<-args[1]

GM12878_GROcap_over_one <- read.table(paste(annotation,sep = ""), header = F)
 data<-data.frame(cbind(c(1:nrow(GM12878_GROcap_over_one)),sort(GM12878_GROcap_over_one[,7])))
 #This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
 calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
 	inputVector <- sort(inputVector)
	inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
	slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
	xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
	y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.

	if(drawPlot){  #if TRUE, draw the plot
		plot(1:length(inputVector), inputVector,type="l",...)
		b <- y_cutoff-(slope* xPt)
		abline(v= xPt,h= y_cutoff,lty=2,col=8)
		points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
		abline(coef=c(b,slope),col=2)
		title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
		axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
	}
	return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
	yPt <- myVector[x]
	b <- yPt-(slope*x)
	xPts <- 1:length(myVector)
	return(sum(myVector<=(xPts*slope+b)))
}
 calculate_cutoff(data[,2], drawPlot=FALSE)$absolute

data[which(data[,2]== calculate_cutoff(data[,2], drawPlot=FALSE)$absolute),][1,1]
colnames(data)<-c("Rank","Reads_counts")

p1<-ggplot(data=data,mapping=aes(x=Rank,y=Reads_counts))+ geom_line(size = 1) + theme_classic() + theme(legend.position = 'top') + guides(color = guide_legend(title = NULL, override.aes = list(size = 1))) + ggtitle(paste('reads within dELS_reads_',annotation))+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +geom_point(aes(
data[which(data[,2]== calculate_cutoff(data[,2], drawPlot=FALSE)$absolute),][1,1], calculate_cutoff(data[,2], drawPlot=FALSE)$absolute
), col="red", size=3)+annotate("text",x = c(0,50000,100000), y =  c(9000,9000,9000), label = c('cut_dot:',data[which(data[,2]== calculate_cutoff(data[,2],drawPlot=FALSE)$absolute),][1,1], calculate_cutoff(data[,2], drawPlot=FALSE)$absolute), size = 3)

pdf(paste(annotation,'.pdf',sep = ""))
print(p1)
dev.off()
write.table(c(data[which(data[,2]== calculate_cutoff(data[,2], drawPlot=FALSE)$absolute),][1,1], calculate_cutoff(data[,2], drawPlot=FALSE)$absolute
),paste(annotation,'_cut_pot.txt',sep = ""))
write.table(calculate_cutoff(data[,2], drawPlot=FALSE)$absolute,paste(annotation,'_dot.txt',sep = ""), col.names= F, row.names= F)
