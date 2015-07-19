#!/bin/sh
# the purpose of this script is to append the SIFT predictions for each mutation in the 'auxotrophy_snps.csv' file
# to do this, I will iterate through each mutation in the 'auxotrophy_snps.csv' file, use the information to look for what population it is in
# and then parse out the SNPs in the 'population.freebayes_SIFTpredictions.vcf' file to get the matching SIFT prediction
# i would like to get the preidiction, probability of toleratance (want <0.05), and coverage (how many matching genomes were used for MSA)

#documentation of SIFT:http://sift-db.bii.a-star.edu.sg/Commandline.html
#http://sift-db.bii.a-star.edu.sg/AboutSIFT4G.html


#usage
 #python Population_Level_Seq_Analysis/parse_SIFT_to_auxo.py 
 #-strainfile /fh/fast/shou_w/NextGenSeq/PopSeq/Parameter+Strain_Files/AaronMiller_strains.txt
 #-auxofile /fh/fast/shou_w/NextGenSeq/PopSeq/MDunham/Miller_SulfateGlucosePhosphate_2015/auxotrophy_snps.csv
 #-dir /fh/fast/shou_w/NextGenSeq/PopSeq/MDunham/Miller_SulfateGlucosePhosphate_2015/

#import libraries needed
import sys
import os
import argparse

#define parsing options
parser = argparse.ArgumentParser()
parser.add_argument('-auxofile', type=argparse.FileType('r'))
parser.add_argument('-dir', dest='dir')

args = parser.parse_args()

# get the path where the data are stored. will be used with other strings later
path=args.dir


ll=args.auxofile.readlines() # read auxotrophy file, parse through for each SNP

p=open('auxotrophy_snps_SIFT.csv','w')

print('Merging SIFT info the auxofile....')
for line in ll:

	line=line.rstrip('\n')
	line_copy=line # will be used later for writing to file a new csv with SIFT predictions
	line=line.split(',')
	#print line
	if(line[0]!='sample'):
		pop=line[0]
		os.chdir(pop)
		#os.system('pwd')
		# find snp being analyzed in SIFTresults.vcf file
		chr=line[1]
		pos=line[2]
		file=pop+'.freebayes_SIFTpredictions.vcf'
		g=open(file,'r')
		gg= g.readlines
		for gline in g:
			gline=gline.rstrip('\n')
			gline=gline.split('\t')
			if(len(gline)>4): 
				g_chr=gline[0]
				g_pos=gline[1]
				if(g_pos==pos):
					if(g_chr==chr):
						#print gline
						#print line # at this point I have a match for the SNP in my auxotrophy file and the vcf file with SIFT predictions
						# PARSE OUT SIFT predictions here
						#print '***********************'
						info=gline[7].split(';')
						SIFT_info=info[-1].split('|')
						#print SIFT_info[3,5,6,7,8,9,10,12]
						#print SIFT_info
						if(len(SIFT_info)>1): # i'm missing SIFT preidcitions for some mutations. Need to compare to figure out what % my new file has compared to auxotrophy_snps.csv file
							# roughly 20 or so mutations with no SIFT annotation. Will ignore for now.
							SIFT_add=SIFT_info[3]+','+SIFT_info[5]+','+SIFT_info[6]+','+SIFT_info[7]+','+SIFT_info[8]+','+SIFT_info[9]+','+SIFT_info[10]+','+SIFT_info[12]+'\n'
				
							line_copy=line_copy+','+SIFT_add
							p.write(line_copy)
							



		#print chr,pos
		
		
		g.close()
		os.chdir('..')
	
	else:
		header=line_copy
		#print header
		#HSV2,NONSYNONYMOUS,D/A,198,0.014,2.73,15,DELETERIOUS
		header=header+',SIFT_Gene,SIFT_MutationType,SIFT_AAChange,SIFT_AAPos,SIFT_Prob,SIFT_MedianSeqConserve,SIFT_Coverage,SIFT_Prediction'+'\n'
		#print header
		p.write(header)


p.close()
