#! Rscript --vanilla
### Input: 1) metafile for the genomeic features
### 2) intmarkerGFoverlapped
### 3) randommarkerGFoverlapped
### 4) number of samples

### Final output files: i) feature specific (p-val_emp; p-val_norm (z-score); # of obs overlapped; expected overlap); ii) cell specific (min p-val_emp; feat_minp_emp; min p-val_norm (z-score); feat_minp_norm; cell p-val_norm_Q; cell p-val_norm_sample_Q; obs_overlapped_average; exp_overlapped_average; median FC ;numberoffeatures)

Args <- commandArgs(TRUE)
if (length(Args)!=5){
        print("Usage: i) metafile; ii) intmarkerGFoverlapped; iii) randommarkerGFoverlapped; iv) number of samples; v) output prefix")
        quit()
}

library(MASS)

metafile <- as.matrix(read.table(Args[1],colClasses="character"));
intmarkerGFoverlapped  <- as.matrix(read.table(Args[2],colClasses="numeric",header=T,check.names=F));
randommarkerGFoverlapped  <- as.matrix(read.table(Args[3],colClasses="numeric",sep="\t",header=T,check.names=F));
numsamples <- as.numeric(Args[4])

### metafile <- as.matrix(read.table("/net/assembly/alextsoi/Researches/Psoriasis_Immunochip_Kiel_GAIN_WTCCC2_Finemapping/MetaAnalysis/r2_common_known_novel5/data/GFmetafile",colClasses="character"));intmarkerGFoverlapped  <- as.matrix(read.table("/net/assembly/alextsoi/Researches/Psoriasis_Immunochip_Kiel_GAIN_WTCCC2_Finemapping/MetaAnalysis/r2_common_known_novel5/function/temp.intmarkerGFoverlapped",colClasses="numeric",header=T,check.names=F));randommarkerGFoverlapped  <- as.matrix(read.table("/net/assembly/alextsoi/Researches/Psoriasis_Immunochip_Kiel_GAIN_WTCCC2_Finemapping/MetaAnalysis/r2_common_known_novel5/function/temp.randommarkerGFoverlapped",colClasses="numeric",sep="\t",header=T,check.names=F));numsamples <- 2000

features <- colnames(intmarkerGFoverlapped)
cells <- unique(metafile[match(metafile[,1],features,nomatch=0)!=0,3])

######## A) feature specific result
print("Computing p-value for each faeture...")
result_f <- matrix(numeric(),length(features),4)
rownames(result_f) <- features
colnames(result_f) <- c("Emp_p","Norm_p","Obs","Exp")
for (f in 1:length(features)){
        tempobs <- intmarkerGFoverlapped[,features[f]]
        tempsample <- randommarkerGFoverlapped[,features[f]]
        ### !!! use adjusted p-value (give half score when samplecount==obscount)
        temp1 <- (sum(tempsample>tempobs)+0.5*sum(tempsample==tempobs))/length(tempsample)
        temp2 <- pnorm((tempobs-mean(tempsample))/sqrt(var(tempsample)),lower.tail=F)
        result_f[f,] <- c(temp1,temp2,tempobs,mean(tempsample))
}

####### B) cell type specific result
print("Computing p-value for each cell type...")
result_c <- matrix(numeric(),length(cells),10)
rownames(result_c) <- cells
colnames(result_c) <- c("minp_emp","minp_emp_feat","minp_norm","minp_norm_feat","Qpval_norm","Qpval_norm_sample","Obs_overlapped_average","Exp_overlapped_average","Median_FC","Number_of_features")

##### !! cannot use the feature if it does not have any variation in the tempsample
availablefeatures <- colnames(randommarkerGFoverlapped)[colSums(randommarkerGFoverlapped)>0]
#####
for (c in 1:length(cells)){
        tempfeatures <- intersect(availablefeatures,metafile[metafile[,3]==cells[c],1])
        tempresult_f <- result_f[tempfeatures,]

        if (length(tempfeatures)>1){
                temp1 <- min(tempresult_f[,1])
                temp2 <- paste(rownames(tempresult_f)[tempresult_f[,1]==temp1],collapse="; ")
                temp3 <- min(tempresult_f[,2])
                temp4 <- paste(rownames(tempresult_f)[tempresult_f[,2]==temp3],collapse="; ")
        } else {
                temp1 <- tempresult_f[1]
                temp2 <- tempfeatures
                temp3 <- tempresult_f[2]
                temp4 <- tempfeatures
        }
        tempobs <- intmarkerGFoverlapped[,tempfeatures]
        tempsample <- randommarkerGFoverlapped[,tempfeatures]

        ### compute observed z-score
        if (class(tempsample)=="matrix"){
                tempobs_z <- (tempobs-colMeans(tempsample))/apply(tempsample,2,function(x){sqrt(var(x))})
        } else {
                tempobs_z <- (tempobs-mean(tempsample))/sqrt(var(tempsample))
        }

        ### compute random z-score using already existing random samples
        if (class(tempsample)=="matrix"){
                tempsample_z <- (tempsample-t(matrix(colMeans(tempsample),dim(tempsample)[2],dim(tempsample)[1])))/t(matrix(apply(tempsample,2,function(x){sqrt(var(x))}),dim(tempsample)[2],dim(tempsample)[1]))
                ###temp5 <- (sum(rowSums(tempsample_z^2) > sum(tempobs_z^2))+0.5*sum(rowSums(tempsample_z^2) == sum(tempobs_z^2)))/dim(tempsample_z)[1]
                temp5 <- (sum(rowSums(tempsample_z) > sum(tempobs_z))+0.5*sum(rowSums(tempsample_z) == sum(tempobs_z)))/dim(tempsample_z)[1]
        } else {
                tempsample_z <- (tempsample-mean(tempsample))/sqrt(var(tempsample))
                ###temp5 <- (sum((tempsample_z^2) > (tempobs_z^2))+0.5*sum((tempsample_z^2) > (tempobs_z^2)))/length(tempsample_z)
                temp5 <- (sum((tempsample_z) > (tempobs_z))+0.5*sum((tempsample_z) > (tempobs_z)))/length(tempsample_z)

        }
        ### compute random z-score by sampling from multivariate normal in which the covariance learned from existing random samples
        if (class(tempsample)=="matrix"){
                tempcov <- cov(tempsample_z)
                tempsample_mvrnormz <- mvrnorm(numsamples,rep(0,length(tempobs_z)),tempcov)
                ###temp6 <- (sum(rowSums(tempsample_mvrnormz^2) > sum(tempobs_z^2))+0.5*sum(rowSums(tempsample_mvrnormz^2) == sum(tempobs_z^2)))/dim(tempsample_mvrnormz)[1]
                temp6 <- (sum(rowSums(tempsample_mvrnormz) > sum(tempobs_z))+0.5*sum(rowSums(tempsample_mvrnormz) == sum(tempobs_z)))/dim(tempsample_mvrnormz)[1]
        } else{
                tempcov <- var(tempsample)
                tempsample_mvrnormz <- rnorm(numsamples)
                ###temp6 <- (sum((tempsample_mvrnormz^2) > (tempobs_z^2))+0.5*sum((tempsample_mvrnormz^2) > (tempobs_z^2)))/length(tempsample_mvrnormz)
                temp6 <- (sum((tempsample_mvrnormz) > (tempobs_z))+0.5*sum((tempsample_mvrnormz) > (tempobs_z)))/length(tempsample_mvrnormz)
        }
        temp7 <- mean(tempobs)
        if (class(tempsample)=="matrix"){
                temp.8 <- colMeans(tempsample)
                temp8 <- mean(colMeans(tempsample))
        } else{
                temp.8 <- tempsample
                temp8 <- mean(tempsample)
        }
        temp9 <- median(tempobs/temp.8)
        temp10 <- length(tempobs)
        result_c[c,] <- c(temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10)
}

write.table(result_f,file=paste(Args[5],"result.feature",sep="."),row.names=T,col.names=T,sep="\t",quote=F)
write.table(result_c,file=paste(Args[5],"result.cell",sep="."),row.names=T,col.names=T,sep="\t",quote=F)
