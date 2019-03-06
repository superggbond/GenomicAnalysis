#! /usr/bin/python
"""
Given the metafiles (ID\tfilelocation\tcell\tFeature) of the genomic features (e.g. from ENCODE), the GWAS info file (marke,chr,start,end,ld-length,maf,numgenes), the interested marker info file, the prefix of the folder containing files of GF overlapping GWAS marker, and number of samplings: perform genomic feature enrichment analysis for individual feature and cell-type based while considering correlations between features within same cell type
"""

from optparse import OptionParser
import sys
import os
from multiprocessing import *
from math import *
sys.path.append(os.path.dirname(os.path.realpath(sys.argv[0])))

def main():
        usage="usage: %prog -m metafile -g gwasinfo -i intmarksinfo -f GF_folder_fileprefix [-s numsamples -t numProcesses] -o outprefx"
        parser=OptionParser(usage)
        parser.add_option("-m", "--metafile",dest="metafile",help="metafile storing the files of the genomic features (ID\tfilelocation\tcell\tFeature)")
        parser.add_option("-g", "--gwasinfo",dest="gwasinfo",help="gwas info file (marke\tchr\tstart\tend\tld-length\tmaf\tnumgenes)")
        parser.add_option("-i", "--intmarksinfo",dest="intmarksinfo",help="interested markers info file (marke\tchr\tstart\tend\tld-length\tmaf\tnumgenes)")
        parser.add_option("-f", "--featurepre",dest="featurepre",help="folder and prefix of the files storing the features overlapping GWAS markers")
        parser.add_option("-s", "--numsamples",dest="numsamples",type="int",default=10000)
        parser.add_option("-t","--numProcess",dest="numProcess",type="int",default=1)
        parser.add_option("-o", "--outDir",dest="outDir",help="output directory and prefix")
        parser.add_option("-p", "--intmarks_ldsnps",dest="intmarks_ldsnps",help="Indicating using ld-markers instead of ld-regions for the interested markers: need to provide a bed file showing each ld-marker per row (4th column indicates the interested marker)")


        (options,args)=parser.parse_args()
        if options.metafile==None or options.gwasinfo==None or options.intmarksinfo==None or options.featurepre==None or options.numsamples==None or options.outDir==None:
                parser.error("Input arguments are not fully specified")

        metafile=options.metafile
        gwasinfo=options.gwasinfo
        intmarksinfo=options.intmarksinfo
        featurepre=options.featurepre
        numsamples=options.numsamples
        numProcess=options.numProcess
        outDir=options.outDir

        if options.intmarks_ldsnps !=None:
                print "Using ld-markers instead of ld-regions for genomic feature overlap.."
                intmarks_ldsnps=options.intmarks_ldsnps

        ###########################################################
        ### A) sample markers
        print "Sampling markers from the GWAS list..."
        tempout=outDir+".sampledmarkers"
        tempscript="Rscript --vanilla "+os.path.dirname(os.path.realpath(sys.argv[0]))+"/samplemarkers.R "+gwasinfo+" "+intmarksinfo+" "+str(numsamples)+" "+tempout
        os.system(tempscript)

        ###########################################################
        ### B) overlap with genomic features
        features=[]
        for tempstr in open(metafile):
                features.append(tempstr.rstrip().split()[0])
        ### i) int markers
        print "Counting numbers of queried markers overlapping with each feature..."
        intmarks_feat_X=[]
        tempout=outDir+".tempout"
        if options.intmarks_ldsnps ==None:
                ### create temporary bed file for int markers
                tempbed=outDir+".intmarkerbed"
                tempscript="cat "+intmarksinfo+" | perl -lane 'print \"chr$F[1]\t$F[2]\t$F[3]\t$F[0]\";' > "+tempbed
                os.system(tempscript)
                for tempstr in open(metafile):
                        tempstr=tempstr.rstrip().split()
                        tempfile=tempstr[1]
                        tempname=tempstr[0]
                        tempscript="bedtools intersect -wa -u -a "+tempbed+" -b "+tempfile+" > "+tempout
                        os.system(tempscript)
                        count=0
                        for tempstr2 in open(tempout):
                                count+=1
                        ### print tempname+": "+str(count)
                        intmarks_feat_X.append(str(count))
        else:
                ### compute overlap using only ld-marks instead of ld-block
                for tempstr in open(metafile):
                        tempstr=tempstr.rstrip().split()
                        tempfile=tempstr[1]
                        tempname=tempstr[0]
                        tempscript="bedtools intersect -wa -u -a "+intmarks_ldsnps+" -b "+tempfile+" | cut -f4 | uniq > "+tempout
                        os.system(tempscript)
                        count=0
                        for tempstr2 in open(tempout):
                                count+=1
                        ### print tempname+": "+str(count)
                        intmarks_feat_X.append(str(count))

        os.system("rm -f "+tempout)
        tempout=open(outDir+".intmarkerGFoverlapped",'w')
        tempout.write('\t'.join(features)+"\n")
        tempout.write('\t'.join(intmarks_feat_X))
        tempout.close()
        
        #############################################################################
        ### ii) sampled markers
        print "Counting numbers of sampled markers overlapping with each feature..."
        tempout=open(outDir+".randommarkerGFoverlapped",'w')
        tempout.write('\t'.join(features))
        ### dict storing the overlapped GWAS markers per feature
        gfmarks={}
        for f in features:
                gfmarks[f]=[]
                for tempstr2 in open(featurepre+"."+f):
                        gfmarks[f].append(tempstr2.rstrip())
        ### random genes
        randmarkslist=[]
        for tempstr in open(outDir+".sampledmarkers"):
                randmarkslist.append(tempstr.rstrip().split())
        ### Multi processing
        if numProcess>1:
                print "Start multi-processing...."
        randgenes_gf_count=[]
        processes=[]
        queues=Queue()
        avgnumrandmarks=int(ceil(numsamples/float(numProcess)))
        culnumrandmarks=0
        for i in range(numProcess):
                if i==(numProcess-1):
                        temprandmarks=randmarkslist[culnumrandmarks:]
                else:
                        temprandmarks=randmarkslist[culnumrandmarks:(culnumrandmarks+avgnumrandmarks)]
                culnumrandmarks+=avgnumrandmarks
                processes.append(Process(target=process_randmarks_gf_overlap,args=(temprandmarks,gfmarks,features,queues)))

        ### start new processes
        for t in processes:
                t.start()

        ### gather the counts
        for t in range(numProcess):
                randgenes_gf_count+=queues.get()

        for r in randgenes_gf_count:
                tempout.write("\n"+'\t'.join(r))
        tempout.close()
        ###################################################################################
        ### C) computed p-values

        tempscript="Rscript --vanilla "+os.path.dirname(os.path.realpath(sys.argv[0]))+"/computep.R "+metafile+" "+outDir+".intmarkerGFoverlapped "+outDir+".randommarkerGFoverlapped "+str(numsamples)+" "+outDir
        os.system(tempscript)

###################################################################
def process_randmarks_gf_overlap(temprandmarks,gfmarks,features,queues):
        tempresult=[]
        count=0
        for r in temprandmarks:
                count+=1
                temptempresult=[]
                print "Processed "+str(count)+" samples in one of the processes..."
                for f in features:
                        temptempresult.append(str(len(list(set(r) & set(gfmarks[f])))))
                tempresult.append(temptempresult)
        queues.put(tempresult)

###################################################################
if __name__ == "__main__":
        main()
