#!/usr/bin/python

#Code created by Alex Tsoi in bash script, modified by Yuhua Zhang to python

import pandas as pd
import os
import re
import sys

def generate_metabam(pathway):
	#input: pathway to store the processed .bam files
	#output: the matebam files
	col1=pd.DataFrame(columns=['1'])
	col2=pd.DataFrame(columns=['2'])
	for file in os.listdir(pathway):
		if(re.search('.noGL.bam\Z',file)):
			col1=col1.append({'1':pathway+'/'+file},ignore_index=True)
			tmp_file=re.sub('.bam\Z','',file)
			col2=col2.append({'2':tmp_file},ignore_index=True)
	#col1=col1.sort(['1'])
	#col2=col2.sort(['2'])
	file=pd.concat([col1,col2],axis=1)
	file=file.sort_values('1')
	return(file)

def generate_metatable(metabam,output_table):
	os.system('python depth_of_coverage.py '+metabam+' '+output_table)

def concat_func(command):
	pathway=command[0]
	outdir=command[1]
	batch=command[2]
	run=command[3]
	#outdir=outdir+'/Batch'+batch+'_Run'+run
	#if not os.path.exists(outdir):
	#	os.system('mkdir '+outdir)
	tmp_file=generate_metabam(pathway)
	filename='metabams_Batch'+batch+'_Run'+run
	tmp_file.to_csv(outdir+'/'+filename,header=None,index=None,sep='\t')
	generate_metatable((outdir+'/'+filename),(outdir+'/'+filename+'.table'))
