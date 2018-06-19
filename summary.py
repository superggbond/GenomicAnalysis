#!/usr/bin/python

#Code contributed by Yuhua Zhang

import pandas as pd
import os
import re
import sys

def get_sampleinfo(out_sampleinfo,batch,run):
	#sample_info=pd.DataFrame(columns=['SampleID','Batch','Condition','Patient'])
	for file in os.listdir(out_sampleinfo):
		if(re.search(('Batch'+batch+'_Run'+run),file)):
			tmp_file=pd.read_table(out_sampleinfo+'/'+file,header=None)
			tmp_file.columns=['SampleID','Batch','Condition','Patient']
			#print(tmp_file)
			#sample_info=sample_info.append(tmp_file,ignore_index=True)
	return(tmp_file)

def get_readcounts(out_proc_bam,batch,run):
	for file in os.listdir(out_proc_bam):
		if(re.search(('readcount_Batch'+batch+'_Run'+run),file)):
			tmp_file=pd.read_table(out_proc_bam+'/'+file)
			tmp_file=tmp_file.set_index('ReadCount').transpose()
			#tmp_file.sort_index(inplace=True)
			#tmp_file=tmp_file[tmp_file.columns[1:]]
			tmp_file.reset_index(drop = True, inplace = True)
	return(tmp_file)

def get_coverage(out_coverage,batch,run):
	for file in os.listdir(out_coverage):
		if (re.search(('_Batch'+batch+'_Run'+run+'.table'),file)):
			tmp_file=pd.read_table(out_coverage+'/'+file)
			#tmp_file.sort_index(inplace=True)
			tmp_file=tmp_file[tmp_file.columns[1:]]
	return(tmp_file)

def get_signaltonoise(out_signaltonoise,batch,run):
	#signaltonoise=pd.DataFrame(columns=['SignalToNoise'])
	for file in os.listdir(out_signaltonoise):
		if (re.search(('Batch'+batch+'_Run'+run+'.signal2noise'),file)):
			#print('label')
			tmp_file=pd.read_table(out_signaltonoise+'/'+file,header=None)
			#tmp_file.sort_index(inplace=True)
			tmp_file=tmp_file[tmp_file.columns[1:]]
			tmp_file.columns=['SignalToNoise']
	return(tmp_file)

def concat_func(command):
	out_sampleinfo=command[0]
	out_proc_bam=command[1]
	out_coverage=command[2]
	out_signaltonoise=command[3]
	batch=command[4]
	run=command[5]
	out_summary=command[6]
	sample_info=get_sampleinfo(out_sampleinfo,batch,run)
	print(sample_info)
	readcount=get_readcounts(out_proc_bam,batch,run)
	print(readcount)
	coverage=get_coverage(out_coverage,batch,run)
	print(coverage)
	signaltonoise=get_signaltonoise(out_signaltonoise,batch,run)
	print(signaltonoise)
	dat=pd.concat([sample_info,readcount,coverage,signaltonoise],axis=1)
	dat.to_csv(out_summary+'/Batch'+batch+'_Run'+run+'.dat',index=None,sep='\t')

if __name__ == '__main__':
	concat_func(command)