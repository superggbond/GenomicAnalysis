#! /usr/bin/Rscript --vanilla 
### Input: 1) meta file of BAM; 2) 
### Final output pdf: 

#Code credit to Alex Tsoi

Args <- commandArgs(TRUE)
if (length(Args)!=2){
        print("Usage: i) BAM files ; ii) output pdf")
        quit()
}

#Args <- c("/net/assembly/alextsoi/Researches/Psoriasis_ATAC_seq_preliminary2/QCs/metabams","temp.pdf")

tempbams <- as.matrix(read.table(Args[1],colClasses="character"))

pdf(Args[2])

for (i in 1:dim(tempbams)[1]){
	print(tempbams[i,1])
	### output insert size
	tempcommand <- paste("samtools view -f66 ",tempbams[i,1]," | cut -f9 | sed 's/^-//' > tempinsert")
	system(tempcommand)
	### only plot insert size < 1,000 bp
	tempname <- strsplit(tempbams[i,1],split="/")[[1]]
	tempname <- tempname[length(tempname)]
	tempinsert <- as.matrix(read.table("tempinsert",colClasses="numeric"))
	hist(tempinsert[tempinsert<1000],breaks=1000,xlab="Insert Size",main=tempname)
}
system(paste("rm -f tempinsert"))
dev.off()
