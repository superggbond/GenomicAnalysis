#! /usr/bin/python2.7
# This script calculates the depth of coverage and breadth of coverage for a given bam. 
# Outputs a dictionary containing the contig/chromosome names and the depth and breadth of coverage for each
# and for the entire genome.
#
# If you optionally specify the name of the mitochondrial chromosome (e.g. mtDNA, chrM, chrMT)
# The script will also generate breadth and depth of coverage for the nuclear genome AND the ratio
# of mtDNA:nuclearDNA; which can act as a proxy in some cases for mitochondrial count within an individual.
# 
# Author: Daniel E. Cook
# Website: Danielecook.com
#
# script modified by Alex Lam C Tsoi to compute genomewide breadth under different coverage cutoff

from sys import *
import os
import re
from subprocess import Popen, PIPE

def get_contigs(bam):
    header, err = Popen(["samtools","view","-H",bam], stdout=PIPE, stderr=PIPE).communicate()
    if err != "":
        raise Exception(err)
    # Extract contigs from header and convert contigs to integers
    contigs = {}
    for x in re.findall("@SQ\WSN:(?P<chrom>[A-Za-z0-9_]*)\WLN:(?P<length>[0-9]+)", header):
        if x[0]=="MT":
		continue
	contigs[x[0]] = int(x[1])
    return contigs

def to_int(x):
    if x=='':
        x=0
    else:
        x=int(x)
    return(x)

def coverage(bam, mtchr = None):
    # Check to see if file exists
    if os.path.isfile(bam) == False:
        raise Exception("Bam file does not exist")
    contigs = get_contigs(bam)

    # Guess mitochondrial chromosome
    mtchr = [x for x in contigs if x.lower().find("m") == 0]
    if len(mtchr) != 1:
        mtchr = None
    else:
        mtchr = mtchr[0]

    coverage_dict = {}
    for c in contigs.keys():
        ### print(c)
        coverage_dict[c] = {}
	### >= 1X	
	command = "samtools depth -r %s %s | awk '{sum+=$3;cnt++}END{print cnt \"\t\" sum}'" % (c, bam)
        coverage_dict[c]["Bases Mapped"], coverage_dict[c]["Sum of Depths"] = map(to_int,Popen(command, stdout=PIPE, shell = True).communicate()[0].strip().split("\t"))

	### >= 2X
	command = "samtools depth -r %s %s | perl -lane 'if ($F[2]>=2){print $_;}' | awk '{cnt++}END{print cnt}'" % (c, bam)
        coverage_dict[c]["Bases Mapped 2X"] = map(to_int,Popen(command, stdout=PIPE, shell = True).communicate()[0].strip().split("\t"))[0]

	### >= 3X
	'''command = "samtools depth -r %s %s | perl -lane 'if ($F[2]>=3){print $_;}' | awk '{cnt++}END{print cnt}'" % (c, bam)
        coverage_dict[c]["Bases Mapped 3X"] = map(int,Popen(command, stdout=PIPE, shell = True).communicate()[0].strip().split("\t"))[0]'''

	### >= 5X
	command = "samtools depth -r %s %s | perl -lane 'if ($F[2]>=5){print $_;}' | awk '{cnt++}END{print cnt}'" % (c, bam)
        coverage_dict[c]["Bases Mapped 5X"] = map(to_int,Popen(command, stdout=PIPE, shell = True).communicate()[0].strip().split("\t"))[0]
	

	"""### >= 10X
	command = "samtools depth -r %s %s | perl -lane 'if ($F[2]>=10){print $_;}' | awk '{cnt++}END{print cnt}'" % (c, bam)
        coverage_dict[c]["Bases Mapped 10X"] = map(int,Popen(command, stdout=PIPE, shell = True).communicate()[0].strip().split("\t"))[0]

        coverage_dict[c]["Breadth of Coverage"] = coverage_dict[c]["Bases Mapped"] / float(contigs[c])
        coverage_dict[c]["Depth of Coverage"] = coverage_dict[c]["Sum of Depths"] / float(contigs[c])
        coverage_dict[c]["Length"] = int(contigs[c])
	"""

    # Calculate Genome Wide Breadth of Coverage and Depth of Coverage
    genome_length = float(sum(contigs.values()))
    coverage_dict["genome"] = {}
    coverage_dict["genome"]["Length"] = int(genome_length)
    coverage_dict["genome"]["Bases Mapped"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k != "genome"])
    coverage_dict["genome"]["Sum of Depths"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k != "genome"])
    coverage_dict["genome"]["Breadth of Coverage"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)
    coverage_dict["genome"]["Breadth of Coverage 2X"] = sum([x["Bases Mapped 2X"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)
    ###coverage_dict["genome"]["Breadth of Coverage 3X"] = sum([x["Bases Mapped 3X"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)
    coverage_dict["genome"]["Breadth of Coverage 5X"] = sum([x["Bases Mapped 5X"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)
    ###coverage_dict["genome"]["Breadth of Coverage 10X"] = sum([x["Bases Mapped 10X"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)

    coverage_dict["genome"]["Depth of Coverage"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)

    if mtchr != None:
        # Calculate nuclear breadth of coverage and depth of coverage
        ignore_contigs = [mtchr, "genome", "nuclear"]
        coverage_dict["nuclear"] = {}
        coverage_dict["nuclear"]["Length"] = sum([x["Length"] for k,x in coverage_dict.iteritems() if k not in ignore_contigs ])
        coverage_dict["nuclear"]["Bases Mapped"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs])
        coverage_dict["nuclear"]["Sum of Depths"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs])
        coverage_dict["nuclear"]["Breadth of Coverage"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs]) / float(coverage_dict["nuclear"]["Length"])
        coverage_dict["nuclear"]["Depth of Coverage"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs]) / float(coverage_dict["nuclear"]["Length"])

        # Calculate the ratio of mtDNA depth to nuclear depth
        coverage_dict["genome"]["mt_ratio"] = coverage_dict[mtchr]["Depth of Coverage"] / float(coverage_dict["nuclear"]["Depth of Coverage"])

    # Flatten Dictionary 
    """coverage = []
    for k,v in coverage_dict.items():
        for x in v.items():
            coverage += [(k,x[0], x[1])]
    return coverage"""
    return coverage_dict["genome"]

if len(argv)!=3:
	print("Usage:\tmetaBams\tbame2breadth.table")
	exit()

### metaBams a file with the location and sampleName
breadth=[]
for tempstr in open(argv[1]):
	tempstr=tempstr.rstrip().split("\t")
	print(tempstr[1])
	### if key: breadth[tempstr[1]]=coverage(tempstr[0])
	breadth.append([tempstr[1],coverage(tempstr[0])])	

outf=open(argv[2],'w')
outf.write("SampleName\t1X\t2X\t5X\n")
###outf.write("SampleName\t1X\t2X\t3X\t5X\n")
###outf.write("SampleName\t1X\t2X\t3X\t5X\t10X")
for k in range(len(breadth)):
	outf.write(breadth[k][0]+"\t"+str(breadth[k][1]["Breadth of Coverage"])+"\t"+str(breadth[k][1]["Breadth of Coverage 2X"])+"\t"+str(breadth[k][1]["Breadth of Coverage 5X"])+"\n")
outf.close()
