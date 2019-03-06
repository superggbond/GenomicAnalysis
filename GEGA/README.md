# GEGA

GEGA is an approach which uses a sampling technique to perform regulatory element enrichment. It can calculate the cell-type specific score by combining the chromatin marks for that cell type. It also allows sampling a small amount and uses z-statistics to estimate the p-value. An example command to run GEGA is as follows:

./GEGA.py -m ../data/GFmetafile -g ../data/common_v3g1keurimputed_r27_Allchr_ldblock_maf_numgenes -i ../data/indep66_Allchr_ldblock_maf_numgenes -f ../data/GEGA_ENCODE_gwas_ldsnp/common_v3g1keurimputed_r27_Allchr_tag0.8 -p ../data/indep66_Allchr_tag0.8snp.bed -s 20000 -t 10 -o ../results/GEGA_20000_snp
