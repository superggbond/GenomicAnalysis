#! /usr/bin/python
"""
Given the metafiles (ID\tfilelocation\tcell\tFeature) of the genomic features (e.g. from ENCODE) and the bed file contains the LD blocks of interested markers, create, for each feature, a file annotating which markers overlap with the feature
"""
from optparse import OptionParser
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(sys.argv[0])))

def main():
        usage="usage: %prog -m metafile -g gwas.bed -o outprefix"
        parser=OptionParser(usage)
        parser.add_option("-m", "--metafile",dest="metafile",help="metafile storing the files of the genomic features (ID\tfilelocation\tcell\tFeature)")
        parser.add_option("-g", "--gwas.bed",dest="gwasbed",help="bed file contains the LD blocks of markers")
        parser.add_option("-o", "--outDir",dest="outDir",help="output directory and prefix")

        (options,args)=parser.parse_args()
        if options.metafile==None or options.gwasbed==None or options.outDir==None:
                parser.error("Input arguments are not fully specified")

        metafile=options.metafile
        gwasbed=options.gwasbed
        outDir=options.outDir

        ############################################################
        tempout=outDir+".temp.bed"
        for tempstr in open(metafile):
                tempstr=tempstr.rstrip().split()
                tempfile=tempstr[1]
                tempname=tempstr[0]
                print tempfile
                tempscript="bedtools intersect -wa -u -a "+gwasbed+" -b "+tempfile+" > "+ tempout
                os.system(tempscript)

                ### only need the identifiers
                tempscript="cut -f4 "+tempout+" | sort | uniq > "+outDir+"."+tempname
                os.system(tempscript)
        os.system("rm -f "+tempout)

###################################################################
if __name__ == "__main__":
        main()
