#! Rscript --vanilla
### Input: 1) gwasinfo
### 2) intmarksinfo
### 3) number of samples

### Final output file: number of samples x markers matrix

Args <- commandArgs(TRUE)
if (length(Args)!=4){
        print("Usage: i) gwasinfo; ii) intmarksinfo; iii) numsamples; iv) output")
        quit()
}

### numsamples <- 10000;gwasinfo <- as.matrix(read.table("../data/common_v3g1keurimputed_r27_Allchr_ldblock_maf_numgenes",colClasses="character",sep="\t")); intmarksinfo <- as.matrix(read.table("../data/indep66_Allchr_ldblock_maf_numgenes",colClasses="character",sep="\t"))

numsamples <- as.numeric(Args[3])
gwasinfo <- as.matrix(read.table(Args[1],colClasses="character",sep="\t"))
intmarksinfo <- as.matrix(read.table(Args[2],colClasses="character",sep="\t"))
results <- matrix(character(),numsamples,dim(intmarksinfo)[1])

for (j in 1:dim(intmarksinfo)[1]){
        tempm1 <- intmarksinfo[j,1]
        templ1 <- as.numeric(intmarksinfo[j,5])
        tempf1 <- as.numeric(intmarksinfo[j,6])
        tempg1 <- as.numeric(intmarksinfo[j,7])
        ### cannot sample itself or markers in same locus
        tempc <- intmarksinfo[j,2]
        temps <- as.numeric(intmarksinfo[j,3])-100000
        tempe <- as.numeric(intmarksinfo[j,4])+100000
        tempgwas <- gwasinfo[(gwasinfo[,1]!=tempm1) & !((as.numeric(gwasinfo[,2])==tempc) & (as.numeric(gwasinfo[,3])>=temps) & (as.numeric(gwasinfo[,4])<=tempe)),]


        ### give higher probability to sample the marker if it has maf/ld length/number genes similar to the interested marker
        tempr1 <- 1/rank(abs(as.numeric(tempgwas[,5])-templ1))
        tempr2 <- 1/rank(abs(as.numeric(tempgwas[,6])-tempf1))
        tempr3 <- 1/rank(abs(as.numeric(tempgwas[,7])-tempg1))
        ### rank product
        temprp <- tempr1*tempr2*tempr3
        ### if apply stringent matching criteria, we will not end up with higher diversity of sampled markers
        ###tempw <- (temprp^2)/(sum(temprp^2))
        tempw <- (temprp)/(sum(temprp))
        results[,j] <- sample(tempgwas[,1],numsamples,prob=tempw,replace=T)
}

write.table(results,file=Args[4],row.names=F,col.names=F,quote=F,sep="\t")
