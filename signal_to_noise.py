#!/usr/bin/python

#Code contributed by Yuhua Zhang

import re
import os
import sys

def generate_clus_parameter(bedprofile):
	#input: file contains the parameter to run bedprofilecounts
	file=open(bedprofile)
	para={}
	for ele in file:
		tmp=ele.rstrip().split('=')
		para[tmp[0]]=tmp[1:]
	file.close()
	return(para)

def generate_clus(bedprofile,dir_bed,dir_clus):
	#input: file contains the parameter to run bedprofilecounts
	#input: pathway to the processed .bed files
	#input: output directory to store the output
	para=generate_clus_parameter(bedprofile)
	for file in os.listdir(dir_bed):
		if (re.search('(.recal.bam.noMT.bed)',file)) or (re.search('(noGL.bam.bed)',file)):
			win_size=para['win_size'][0]
			step_size=para['step_size'][0]
			regions=para['regions'][0]
			norm_flag=para['norm_flag'][0]
			read_type=para['read_type'][0]
			avg_frag_len=para['avg_frag_len'][0]
			count_type=para['count_type'][0]
			chr_sizes=para['chr_sizes'][0]
			os.system('/net/fantasia/home/alextsoi/code/ATAC_Seq/ChangLab/atacseq_tools_changrila/bin/bedProfileCount '
				+win_size+' '+step_size+' '+regions+' '+norm_flag+' '+read_type+' '+avg_frag_len+' '+count_type+' '+chr_sizes+' '+dir_bed+'/'+file)
	os.system('mv ./*.clus '+dir_clus)

def generate_signal_to_noise(core_info,outdir,dir_clus):
	#input: core_info file generated in the previous step
	#input: output directory to store the generated pdf file
	#input: pathway to dir_clus file
	os.system('Rscript signal_to_noise.R '+core_info+' '+outdir+' '+dir_clus)

def concat_func(command):
	bedprofile=command[0]
	dir_bed=command[1]
	dir_clus=command[2]
	core_info=command[3]
	dir_pdf=command[4]
	#generate .clus file
	generate_clus(bedprofile,dir_bed,dir_clus)

	#generate s2n plot
	generate_signal_to_noise(core_info,dir_pdf,dir_clus)