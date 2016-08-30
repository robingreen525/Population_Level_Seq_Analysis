#!/bin/sh

#import libraries needed
import sys
import os
import argparse

#define parsing options
parser = argparse.ArgumentParser()
parser.add_argument('-strainfile', type=argparse.FileType('r'))
parser.add_argument('-dir', dest='dir')

args = parser.parse_args()

# get the path where the data are stored. will be used with other strings later
path=args.dir


strains=[]
### use the strain file flag to get the list of strains to analyze
ll=args.strainfile.readlines()
for line in ll:
	line=line.split('\t')
	strains.append(line[0])
	

for strain in strains:
	full_path=path+strain
	print full_path
	vcf=full_path+'/'+strain+'.freebayes.vcf'
	
	#java -jar <Path to SIFT4G_Annotator> -c -i <Path to input vcf file> -d <Path to SIFT4G database directory> -r <Path to your results folder> -t
	w='java -jar /fh/fast/shou_w/bin/SIFT4G_Annotator_v2.2.jar -r'+full_path+ \
	' -d /fh/fast/shou_w/NextGenSeq/reference/EF4.74/ -c -i '+vcf
	print(w)
	os.system(w)
