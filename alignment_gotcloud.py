#!/usr/bin/python

#Code contributed by Yuhua Zhang

import pandas as pd
import os
import re
import sys
import itertools
import multiprocessing

'''def get_batch_run(pathway):
	#get the batch and run from the output directory
	batch=re.search('(Batch\d*\d)',pathway).group(1)
	run=re.search('(Run\d*\d)',pathway).group(1)
	return(batch,run)'''

def get_sampleID(pathway):
	#get the sampleID from the pathway specified by user
	SampleID=[]
	for file in os.listdir(pathway+'/elder'):
		if (re.search('Sample*',file)):
			SampleID.append(re.sub('Sample_','',file))
	SampleID.sort()
	return(SampleID)

def fastq_index(pathway,index_pathway,batch,run):
	#generate the input list of fastq files
	#input: pathway of trimmed data specified by user
	#output: fastq_index file
	SampleID=get_sampleID(pathway)
	#(batch,run)=get_batch_run(pathway)
	lanes=[]
	file1=[]
	file2=[]
	for ele in SampleID:
		lane=0
		for file in os.listdir(pathway+'/elder/Sample_'+ele):
			if(re.search('R1_\d*\d.trimmed.fastq.gz',file)):
				lane=lane+1
				file1.append(file)
			if(re.search('R2_\d*\d.trimmed.fastq.gz',file)):
				file2.append(file)
		lanes.append(lane)
	
	for i in range(len(SampleID)):
		MERGE_NAME=pd.DataFrame(columns=['MERGE_NAME'])
		FASTQ1=pd.DataFrame(columns=['FASTQ1'])
		FASTQ2=pd.DataFrame(columns=['FASTQ2'])
		for lane in range(lanes[i]):
			MERGE_NAME=MERGE_NAME.append({'MERGE_NAME':SampleID[i]},ignore_index=True)
			File1=file1[sum(lanes[0:i])+lane]
			File2=file2[sum(lanes[0:i])+lane]
			FASTQ1=FASTQ1.append({'FASTQ1':pathway+'/elder/Sample_'+SampleID[i]+'/'+File1},ignore_index=True)
			FASTQ2=FASTQ2.append({'FASTQ2':pathway+'/elder/Sample_'+SampleID[i]+'/'+File2},ignore_index=True)
		fastq_index=pd.concat([MERGE_NAME,FASTQ1,FASTQ2],axis=1)
		fastq_index_name=index_pathway+'/fastq_Batch'+batch+'_Run'+run+'_index_'+SampleID[i]
		fastq_index.to_csv(fastq_index_name,index=None,sep='\t')

'''def generate_config_file(outdir,index_pathway,pathway,batch,run):
	#generate the configuration file
	#input: output directory, pathway to index file
	#output: configuration file
	#(batch,run)=get_batch_run(pathway)
	SampleID=get_sampleID(pathway)
	for ele in SampleID:
		filename='gotcloud_Batch'+batch+'_Run'+run+'_config_file_'+ele
		file=open(index_pathway+'/'+filename,'a')
		file.write('FASTQ_LIST = '+index_pathway+'/fastq_Batch'+batch+'_Run'+run+'_index_'+ele+'\n')
		file.write('OUT_DIR = '+outdir+'/Sample_'+ele+'\n')
		file.write('#BATCH_TYPE = mosix\n')
		file.write('#BATCH_OPTS = -j21,22,23,24,30\n')
		file.write('############\n')
		file.write('# References\n')
		file.write('REF_DIR = /data/local/ref/gotcloud.ref\n')
		file.write('AS = NCBI37\n')
		file.write('#REF = $(REF_DIR)/human.g1k.v37.fa\n')
		file.write('REF = /net/1000g/mktrost/seqshop/gotcloud/gotcloud.ref/human.g1k.v37.fa\n')
		file.write('DBSNP_VCF = $(REF_DIR)/dbsnp_135.b37.vcf.gz\n')
		file.write('HM3_VCF = $(REF_DIR)/hapmap_3.3.b37.sites.vcf.gz\n')
		file.write('#\n')
		file.write('MAP_TYPE = BWA_MEM\n')
		file.write('BAMUTIL_THINNING = --phoneHomeThinning 0\n')
		file.close()'''

def run_gotcloud(param):
	outdir=param[0]
	index_pathway=param[1]
	ele=param[2]
	batch=param[3]
	run=param[4]
	conf=param[5]
	print('Processing alignment of Sample_'+ele)
	#(batch,run)=get_batch_run(pathway)
	os.system('mkdir '+outdir+'/Sample_'+ele)
	os.system('gotcloud align --conf '+conf+' --outdir '+outdir+'/Sample_'+ele
		+' --list '+index_pathway+'/fastq_Batch'+batch+'_Run'+run+'_index_'+ele)
	print('Alignment of Sample_'+ele+'  done')

def generate_meta_file(outdir,index_pathway,pathway,batch,run):
	SampleID=get_sampleID(pathway)
	#(batch,run)=get_batch_run(pathway)
	Flag=False
	for bam_file in os.listdir(outdir):
		if (re.search('Sample_',bam_file)):
			Flag=True
			break
	if Flag:
		file=open(index_pathway+'/metagotCloudbamfiles_Batch'+batch+'_Run'+run,'a')
		for ele in SampleID:
			file.write(outdir+'/Sample_'+ele+'/bams/'+ele+'.recal.bam\n')
		file.close()
	else:
		file=open(index_pathway+'/metagotCloudbamfiles_Batch'+batch+'_Run'+run,'a')
		for ele in SampleID:
			file.write(outdir+'/bams/'+ele+'.recal.bam\n')
		file.close()

def concat_func(command):
	'''outdir=''
	pathway=''
	index_pathway=os.path.dirname(os.path.realpath(__file__))
	process_input=3
	try:
		opts, args=getopt.getopt(argv,'o:p:d:j:h',['output=','pathway=','directory=','job=','help'])
	except getopt.GetoptError:
		print('python alignment_gotcloud.py -p <pathway_to_trimmed_files> -o <output_directory>')
		sys.exit(2)

	for opt, arg in opts:
		if opt in ('-h', '--help'):
			print('python alignment_gotcloud.py -p <pathway_to_trimmed_files> -o <output_directory> -d <directory_to_write_index&config_file> -j <job>')
			print('Please do not include / at the end of the pathway')
			print('The default directory to write index and config file is current directory')
			print('Please pay attention to the memory when specifying more than 1 job')
			sys.exit(2)
		elif opt in ('-p','--pathway'):
			pathway=arg
		elif opt in ('-o','--output'):
			outdir=arg
		elif opt in ('-d','--directory'):
			index_pathway=arg
		elif opt in ('-j','--job'):
			process_input=arg
		else:
			print('python alignment_gotcloud.py -p <pathway_to_trimmed_files> -o <output_directory>')
			sys.exit(2)'''
	print('Processing alignment through gotcloud...')
	outdir=command[0]
	pathway=command[1]
	index_pathway=command[2]
	process_input=command[3]
	batch=command[4]
	run=command[5]
	conf=command[6]
	
	#write the index files
	fastq_index(pathway,index_pathway,batch,run)

	#write the config files
	#generate_config_file(outdir,index_pathway,pathway,batch,run)

	#run gotcloud separately for each sample
	SampleID=get_sampleID(pathway)
	paramlist=list(itertools.product([outdir],[index_pathway],SampleID,[batch],[run],[conf]))
	process_input=int(process_input)
	pool=multiprocessing.Pool(processes=process_input)
	pool.map(run_gotcloud,paramlist)

	#generate metacloud file
	generate_meta_file(outdir,index_pathway,pathway,batch,run)

	print('Alignment done')

if __name__=="__main__":
	concat_func(command)