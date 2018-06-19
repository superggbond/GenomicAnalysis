"""
### processed gotCloud bam files for ATAC-seq and report # of reads in each step

Input: i) list of files for gotCloud generated bam 

Output: i) output directory for processed bams; ii) output table
"""

#Code credit to Alex Tsoi, modified by Yuhua Zhang to allow multiprocessing

from sys import *
import os
import subprocess as sub
import pandas as pd
import itertools
import multiprocessing
import re

'''if len(argv)!=4:
	print("Usage:\tmetagotCloudbam\toutdir\toutputtable\tjob\n")
	exit()'''

def proc_sum_gotCloud(param):
	tempstr=param[0]
	outdir=param[1]
	intermediate_file=param[2]
	tempstr=tempstr.rstrip()
	tempsample=tempstr.rstrip().split("/")[-1]
	print("Processing " +tempsample+" ...")	
	readcounts=pd.DataFrame(columns=[tempsample])
	### # of raw reads (there could be multiple alignments for each read, so only need the primary ones, and -F 256 would include non-mapped one as well)
	tempcount=os.popen("samtools view -F 256 -c "+tempstr).readline().rstrip()
	readcounts=readcounts.append({tempsample:tempcount},ignore_index=True)
	### # of mapped reads
	tempcount=os.popen("samtools view -c -F 260 "+tempstr).readline().rstrip()
	readcounts=readcounts.append({tempsample:tempcount},ignore_index=True)
	### # of self- and mate- mapped reads
	#tempcount=os.popen("samtools view -c -F 268 "+tempstr).readline().rstrip()
	#readcounts=readcounts.append({tempsample:tempcount},ignore_index=True)
	### # of self- and mate- mapped reads which are on the same chromosome
	#tempcommand="samtools view -F 268 "+tempstr+" | awk '($7==\"=\")' | wc -l"
	#tempcount=sub.Popen(tempcommand,stdout=sub.PIPE,stderr=sub.PIPE,shell=True).communicate()[0].rstrip()
	#readcounts=readcounts.append({tempsample:tempcount},ignore_index=True)

	### # of self- and mate- mapped reads which are on the same chromosome, oriented towards each other, and with a sensible insert size
	tempcommand="samtools view -F 268 -f 2 "+tempstr+" | awk '($7==\"=\")' | wc -l"	
	tempcount=sub.Popen(tempcommand,stdout=sub.PIPE,stderr=sub.PIPE,shell=True).communicate()[0].rstrip()
	readcounts=readcounts.append({tempsample:tempcount},ignore_index=True)	
	### # of reads after removing mitochondria reads
	os.system("samtools view -h -F 268 -f 2 "+tempstr+" | awk '($3!=\"MT\")' | samtools view -bS - > "+outdir+"/"+tempsample+".mapped.noMT.bam ")
	tempcount=os.popen("samtools view -c "+outdir+"/"+tempsample+".mapped.noMT.bam ").readline().rstrip()
	readcounts=readcounts.append({tempsample:tempcount},ignore_index=True)
	### # of reads after removing duplicated reads
	os.system("bam squeeze --in "+outdir+"/"+tempsample+".mapped.noMT.bam --out "+outdir+"/"+tempsample+".dedup.mapped.noMT.bam --keepOQ")
	tempcount=os.popen("samtools view -c -F 268 -f 2 "+outdir+"/"+tempsample+".dedup.mapped.noMT.bam ").readline().rstrip()
	os.system("rm -f "+outdir+"/temp."+tempsample+".dedup.mapped.noMT.bam")
	readcounts=readcounts.append({tempsample:tempcount},ignore_index=True)
	### clipOverlap 
	os.system("bam clipOverlap --in "+outdir+"/"+tempsample+".dedup.mapped.noMT.bam --out "+outdir+"/"+tempsample+".dedup.mapped.noMT.clipped.bam --storeOrig")
	### Mapping quality 30
	os.system("samtools view -h -F 268 -f 2 -q 30 "+outdir+"/"+tempsample+".dedup.mapped.noMT.clipped.bam | samtools view -bS - > "+outdir+"/"+tempsample+".dedup.mapped.noMT.clipped.MQ.bam")
	tempcount=os.popen("samtools view -c "+outdir+"/"+tempsample+".dedup.mapped.noMT.clipped.MQ.bam").readline().rstrip()
	readcounts=readcounts.append({tempsample:tempcount},ignore_index=True)
	### Remove reads mapped to GL chromosome
	os.system("samtools view -h "+outdir+"/"+tempsample+".dedup.mapped.noMT.clipped.MQ.bam | awk ' ($3!~ \"GL\")' | samtools view -bS - > "+outdir+"/temp."+tempsample+".dedup.mapped.noMT.clipped.MQ.noGL.bam")
	os.system("samtools sort -n "+outdir+"/temp."+tempsample+".dedup.mapped.noMT.clipped.MQ.noGL.bam | samtools fixmate -O bam - - | samtools view -h -f 1 - | samtools sort - -o "+outdir+"/"+tempsample+".dedup.mapped.noMT.clipped.MQ.noGL.bam ")
	tempcount=os.popen("samtools view -c "+outdir+"/"+tempsample+".dedup.mapped.noMT.clipped.MQ.noGL.bam").readline().rstrip()
	os.system("rm -f "+outdir+"/temp."+tempsample+".dedup.mapped.noMT.clipped.MQ.noGL.bam")
	readcounts=readcounts.append({tempsample:tempcount},ignore_index=True)
	### Index
	os.system("samtools index "+outdir+"/"+tempsample+".dedup.mapped.noMT.clipped.MQ.noGL.bam")
	if intermediate_file==False:
		os.system('rm '+outdir+'/'+tempsample+'.mapped.noMT.bam')
		os.system('rm '+outdir+'/'+tempsample+'.dedup.mapped.noMT.bam')
		os.system('rm '+outdir+'/'+tempsample+'.dedup.mapped.noMT.clipped.bam')
		os.system('rm '+outdir+'/'+tempsample+'.dedup.mapped.noMT.clipped.MQ.bam')
	readcounts.to_csv(tempsample,index=None,sep='\t')

def generate_meta_file(outdir,filename,batch,run):
	file=open(outdir+'/metagotCloudbamfiles_filtered_Batch'+batch+'_Run'+run,'a')
	for ele in filename:
		file.write(outdir+'/'+ele+'.dedup.mapped.noMT.clipped.MQ.noGL.bam\n')
	file.close()

def concat_func(command):
	print('Processing filtering...')
	# gererate the output except the output table
	outdir=command[1]
	intermediate=command[3]
	batch=command[4]
	run=command[5]
	tempstring=[]
	filename=[]
	for ele in open(command[0]):
		tempstring.append(ele)
		temp=ele.rstrip()
		filename.append(temp.rstrip().split("/")[-1])
	paramlist=list(itertools.product(tempstring,[outdir],[intermediate]))
	job=int(command[2])
	pool=multiprocessing.Pool(job)
	pool.map(proc_sum_gotCloud,paramlist)

	# generate the output table
	outtable=pd.read_table(filename[0])
	os.system('rm '+filename[0])
	for ele in filename[1:]:
		tmp=pd.read_table(ele)
		outtable=pd.concat([outtable,tmp],axis=1)
		os.system('rm '+ele)
	tmp={'ReadCount':['Raw','Mapped','PairProperlymapped','RemovedMito','Duplicate','MQ30','noGL']}
	tmp=pd.DataFrame(data=tmp)
	outtable=pd.concat([tmp,outtable],axis=1)
	outtable.to_csv(command[1]+'/readcount_Batch'+batch+'_Run'+run,index=None,sep='\t')

	#generate meta file foor plotting
	generate_meta_file(outdir,filename,batch,run)

	print('Filtering done')

if __name__=="__main__":
	concat_func(command)