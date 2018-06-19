#Code credit to Alex
#input: 1)SampleInfo file 2)output dirctory to store pdf/table 3)pathway to .clus file(no need to specify batch,run)
Args <- commandArgs(TRUE)

SampleInfo<-as.matrix(read.table(Args[1],colClasses="character",sep="\t"))

tmp<-tail(sapply(strsplit(Args[1],'/'),'[',-1),1)
batch<-sapply(strsplit(tmp,'_'),'[',2)
run<-sapply(strsplit(tmp,'_'),'[',3)

Signal2Noise_TSS <- matrix(numeric(),dim(SampleInfo)[1],2)
colnames(Signal2Noise_TSS) <- c("Unfiltered","Filtered")

pdf(paste(Args[2],"/ATAC_Signal2Noise_TSS_",batch,"_",run,"_unfiltered_filtered.pdf",sep=""))
for (i in 1:dim(SampleInfo)[1]){
  print(i)
  ### unfiltered data
  tempunfiltered <- as.matrix(read.table(paste(Args[3],
  "/HOMER_hg19_TSS_nochr.bed_",SampleInfo[i,1],".recal.bam.noMT.bed.clus",sep=""),colClasses="character",skip=1))
  tempunfiltered <- tempunfiltered[,3:dim(tempunfiltered)[2]]
  mode(tempunfiltered) <- "numeric"
  
  ### remove the 1% outlier
  tempcutoff <- round(dim(tempunfiltered)[1]*0.99,0)
  tempunfiltered <- tempunfiltered[order(rowSums(tempunfiltered)),]         
  tempunfiltered <- tempunfiltered[1:tempcutoff,]
  tempmeans <- colMeans(tempunfiltered)
  
  ### using the first 5 and the last 5 bins as background
  tempnormalized <- tempmeans/mean(c(tempmeans[1:5],tempmeans[(length(tempmeans)-4):length(tempmeans)]))

  ### using TSS as signal
  Signal2Noise_TSS[i,1] <- mean(tempnormalized[200:201])
  
  ### filtered data
  tempfiltered <- as.matrix(read.table(paste(Args[3],"/HOMER_hg19_TSS_nochr.bed_",
  SampleInfo[i,1],".recal.bam.dedup.mapped.noMT.clipped.MQ.noGL.bam.bed.clus",sep=""),colClasses="character",skip=1))
  tempfiltered <- tempfiltered[,3:dim(tempfiltered)[2]]
  mode(tempfiltered) <- "numeric"
  
  ### remove the 1% outlier
  tempcutoff <- round(dim(tempfiltered)[1]*0.99,0)
  tempfiltered <- tempfiltered[order(rowSums(tempfiltered)),]
  tempfiltered <- tempfiltered[1:tempcutoff,]
  tempmeans <- colMeans(tempfiltered)
  
  ### using the first 5 and the last 5 bins as background
  tempnormalized2 <- tempmeans/mean(c(tempmeans[1:5],tempmeans[(length(tempmeans)-4):length(tempmeans)]))
  
  ### using TSS as signal
  Signal2Noise_TSS[i,2] <- mean(tempnormalized2[200:201])
  
  ### plot
  tempymax = 1.1*max(c(tempnormalized,tempnormalized2))
  names(tempnormalized) <- seq(-2000,1999,10)
  plot(tempnormalized,type="l",cex=0.2,col="black",main=paste("ATAC TSS enrichment",SampleInfo[i,1]),ylim=c(0,tempymax),xaxt="n",xlab="Relative position to TSS")
  axis(1,at=seq(0,400,50),labels=seq(-2000,2000,500))
  lines(tempnormalized2,type="l",cex=0.2,col="blue")
}
dev.off()

write.table(cbind(apply(SampleInfo,1,function(x){paste(x[2],x[1],x[3],x[4],sep=":")}),Signal2Noise_TSS[,2]),file=paste(Args[2],"/",batch,"_",run,".signal2noise",sep=""),row.names=F,col.names=F,sep="\t",quote=F)
