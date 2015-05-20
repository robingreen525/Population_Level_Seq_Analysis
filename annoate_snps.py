# the purpose of this script is to annotate each SNP with the gene identifier (standard and systematic name)
# along with what codon and position in the codon. This information can be used to determine which frame 
# the mutation is in and predict amino acid substitutions  

#import libraries needed
import sys
import os
import argparse

#define parsing options
parser = argparse.ArgumentParser()
parser.add_argument('-genefile',type=argparse.FileType('r'))
parser.add_argument('-strainfile', type=argparse.FileType('r'))
parser.add_argument('-dir', dest='dir')

args = parser.parse_args()


strains=[]
### use the strain file flag to get the list of strains to analyze
ll=args.strainfile.readlines()
for line in ll:
	line=line.split('\t')
	strains.append(line[0])
	
# build a lookup table of all genes in cerevisea genome	
ll=args.genefile.readlines()	
lookup_table=[] # will be a list of lists where each sublist is the gene information


for line in ll:
	f=line[0:3] #first three letters of line, want to get 'chr' starting lines in .gff file
	if(f =='chr'):
		line=line.split()
		if(line[2]=='gene'):
			chr=line[0] #chromosome
			start=line[3] # start site (relative to plus strand)
			end=line[4] # end site (relative to plus strand)
			strand=line[6] # what strand, if not plus, start and end will be flipped later
			info=line[8].split(';') # where the ID and gene name information is
			for thing in info:
				thing=thing.split('=')
				if(thing[0]=='ID'): ID=thing[1] # systematic name
				elif(thing[0]=='gene'):  # common/standard name
					gene=thing[1]
					
			if 'gene' in globals():	## see if gene name was assigned along with systematic name
				lookup_table.append([chr,start,end,strand,ID,gene])
			else:
				lookup_table.append([chr,start,end,strand,ID,'unknown'])