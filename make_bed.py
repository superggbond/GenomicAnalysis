#!/usr/bin/python

#Code originally credit to Alex Tsoi in bash script, rewrite to python by Yuhua Zhang

import re
import os
import sys

def generate_bed_unfiltered(pathway,outdir):
	#input: pathway to .bam files
	Flag=False
	for bam_file in os.listdir(pathway):
		if (re.search('Sample_',bam_file)):
			Flag=True
			break
	if Flag:
		for file in os.listdir(pathway):
			for sub_dir in os.listdir(pathway+'/'+file+'/bams'):
				if(re.search('.recal.bam\Z',sub_dir)):
					filename=sub_dir
					os.system('bedtools bamtobed -i '+pathway+'/'+file+'/bams/'+sub_dir+' > '+outdir+'/'+filename+'.bed')
	else:
		for file in os.listdir(pathway+'/bams'):
			if (re.search('.recal.bam\Z',file)):
				filename=file
				os.system('bedtools bamtobed -i '+pathway+'/bams/'+file+' > '+outdir+'/'+filename+'.bed')

def generate_bed_filtered(pathway,outdir):
	#input: pathway to filtered .bam files
	for file in os.listdir(pathway):
		if (re.search('.noGL.bam\Z',file)):
			filename=file
			os.system('bedtools bamtobed -i '+pathway+'/'+file+' > '+outdir+'/'+filename+'.bed')

def remove_MT(outdir):
	#input: dir to store the bed files
	for file in os.listdir(outdir):
		if (re.search('recal.bam.bed',file)):
			tmp_file=re.sub('.bed','',file)
			os.system('grep -v ^[MGH] '+outdir+'/'+file+' > '+outdir+'/'+tmp_file+'.noMT.bed')

def concat_func(command):
	pathway1=command[0]
	pathway2=command[1]
	outdir=command[2]
	generate_bed_unfiltered(pathway1,outdir)
	generate_bed_filtered(pathway2,outdir)
	remove_MT(outdir)
