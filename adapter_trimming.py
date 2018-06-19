#!/usr/bin/python

#Code contributed by Yuhua Zhang

import sys
import pandas as pd
import re
import glob
import os
import itertools	
import multiprocessing

def parse_file(filename):
	#Generate patient ID and description
	#input: core_info_file
	#output: patient ID and description(cell type etc)
	core_info=pd.read_table(filename)
	tmp=core_info[core_info['Lane']==1]
	tmp=tmp.drop(tmp.index[len(tmp)-1])

	#patient ID
	patient_ID=tmp['Description']
	patient_ID=pd.Series(patient_ID).str.extract('(\A\d*\d)',expand=True)
	
	cell_type=tmp['Description']
	cell_type=pd.Series(cell_type).str.extract('(CD4[a-zA-Z]+|CD8[a-zA-Z]+|mDC)',expand=True)
	cell_type=cell_type.apply(lambda x:x.astype(str).str.upper())
	cell_type=cell_type.apply(lambda x:x.astype(str).str.replace('MDC','mDC',case=False))

	hour=tmp['Description']
	hour=pd.Series(hour).str.extract('(\d*h)',expand=True)
	
	operation=tmp['Description']
	operation=pd.Series(operation).str.extract('(s[a-zA-Z]*|c[a-zA-Z]*)',expand=True)

	#Generate the description
	descp=pd.concat([cell_type,hour,operation],axis=1)
	descp=descp.apply(lambda x:x.astype(str).str.cat(sep='_'),axis=1)
	descp=descp.apply(lambda x:re.sub('_nan','',x))
	return(patient_ID,descp)

def generate_sample_info(filename,pathway):
	#combine with the information from parse_file and generate sample_info
	#input: core_info_file, pathway provided bu user
	#output: sample info
	(patient_ID,Description)=parse_file(filename)
	(lanes,samples,sampleID)=get_lanes_samples(filename)
	tmp=re.search('(Run_\d*\d)',filename).group(1)
	tmp=re.sub('Run_','',tmp)
	batch=re.search('(Batch\d*\d)',filename).group(1)
	sample_info_name=batch+'_Run'+tmp
	tmp=sample_info_name
	batch_run=pd.DataFrame(columns=['A'])
	for i in range(samples):
		batch_run=batch_run.append({'A':tmp},ignore_index=True)
	sample_info=pd.concat([sampleID,batch_run,Description,patient_ID],axis=1)
	sample_info_name='SampleInfo_'+sample_info_name
	return(sample_info,sample_info_name)

def get_lanes_samples(filename):
	#get sample number and lanes
	#input: core_info_file
	#output: lanes, sample numbers and sampleID
	core_info=pd.read_table(filename)
	lanes=len(core_info.groupby(['Lane']).count())
	samples=len(core_info.groupby(['SampleID']).count())-1
	tmp=core_info[core_info['Lane']==1]
	tmp=tmp.drop(tmp.index[len(tmp)-1])
	sampleID=tmp['SampleID']
	return(lanes,samples,sampleID)

def check_sampleID(filename,pathway):
	#check whether samples in the folder is consistant with the core_info
	#input: core_info_file, pathway provoded by user
	#output: a bool variable; if passed, true, vice versa
	(lanes,samples,sampleID)=get_lanes_samples(filename)
	Flag=True
	Samples=len(glob.glob(pathway+'/elder/Sample*'))
	#First compare number of samples
	if samples!=Samples:
		print("Numder of Samples is not correct")
		Flag=False

	else:
		print("Number of Samples is correct")
		#Next, compare sampleID
		SampleID=[]
		for file in os.listdir(pathway+'/elder'):
			if (re.search('Sample*',file)):
				SampleID.append(re.sub('Sample_','',file))
		SampleID.sort()
		if set(sampleID)!=set(SampleID):
			print("The SampleID doesn't match")
			Flag=False

		else:
			print("The SampleID matches")
			for ele in SampleID:
				count_=0
				for file in os.listdir(pathway+'/elder/Sample_'+str(ele)+'/'):
					if (re.search('\d.fastq.gz',file)):
						count_+=1
				#Compare Lanes
				if count_!=2*lanes:
					print("Number of lanes in Sample"+str(ele)+" doesn't match")
					Flag=False
					break
				else: 
					print('Number of lanes in Sample'+str(ele)+" matches")
	return(Flag)

def generate_trimmed_file_para(param):
	pathway=param[0]
	i=param[1]
	ele=param[2]
	os.system('python /net/fantasia/home/alextsoi/Software/atactk/scripts/trim_adapters '
		+pathway+'/elder/Sample_'+ele+'/*L00'+str(i+1)+'*R1*.fastq.gz '
		+pathway+'/elder/Sample_'+ele+'/*L00'+str(i+1)+'*R2*.fastq.gz')

def concat_func(command):
	#read in the file specified by user
	#first check if the files match with core_info
	#then generate sample_info and adapter-trimmed file
	
	'''pathway=''
	filename=''
	process_input=''
	try:
		opts, args=getopt.getopt(argv,'p:f:j:h',['pathway=','filename=','job=','help'])
	except getopt.GetoptError:
		print('python parse_core_info.py -p <pathway> -f <core_info_file> -j <job_in_parallel>')
		sys.exit(2)

	for opt, arg in opts:
		if opt in ('-h','--help'):
			print('python parse_core_info.py -p <pathway> -f <core_info_file> -j <job_in_parallel>')
			print('Please do not include / at the end of the pathway')
			sys.exit(2)
		elif opt in ('-p','--pathway'):
			pathway=arg
		elif opt in ('-f','--filename'):
			filename=arg
		elif opt in ('-j','--job'):
			process_input=arg
		else:
			print('python parse_core_info.py -p <pathway> -f <core_info_file> -j <job_in_parallel>')
			sys.exit(2)
	'''
	print('Processing adapter trimming...')
	pathway=command[0]
	filename=command[1]
	process_input=command[2]
	entire_output=command[3]
	Flag=check_sampleID(filename,pathway)
	if Flag:
		(sample_info,sample_info_name)=generate_sample_info(filename,pathway)
		sample_info.to_csv(entire_output+'/'+sample_info_name,header=None,index=None,sep='\t')
		#call trim_adapter in parallel
		(lanes,sample,SampleID)=get_lanes_samples(filename)
		lane=range(lanes)
		path=[pathway]
		paramlist=list(itertools.product(path,lane,SampleID))
		process_input=int(process_input)
		pool=multiprocessing.Pool(processes=process_input)
		pool.map(generate_trimmed_file_para,paramlist)
		print('Adapter trimming done')

if __name__=="__main__":
	concat_func(command)