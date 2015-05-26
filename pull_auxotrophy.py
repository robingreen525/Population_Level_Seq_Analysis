# the purpose of this script is to extract the auxotrophy mutations in the 'master.csv'
#file that I have generated so I restrict my analysis only to SNPs that could potentially 
# result in auxotrophy.

#import libraries needed
import sys
import os
import argparse

#define parsing options
parser = argparse.ArgumentParser()
parser.add_argument('-auxofile',type=argparse.FileType('r'))
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
	
	
auxotrophy={}	
### use the auxofile flag to build a table of auxotrophic genees
ll=args.auxofile.readlines()
for line in ll:
	line=line.split('\t')
	if len(line)==10: # if non-header information line
		#print line
		#print line[1],line[7],line[-1].rstrip('\n')
		auxotrophy[line[1]]=[line[7],line[-1].rstrip('\n')]
		
file=path+'/master.csv'
f=open(file,'r')
ll=f.readlines()

newfile=path+'/auxotrophy_snps.csv'
a=open(newfile,'w')

a.write('sample,chr,position,ref_allele,ref_reads,alt_allele,alt_reads,strand,SysName,gene,freq,Mutation,ref_codon,alt_codon,frame,auxotrophy,ref\n')


#S1,chrXIV,196580,T,38,G,5,-,YNL241C,ZWF1,0.116279069767,E455A,gag,gcg,2,methionine,Slekar KH, et al. (1996) PMID:8910528


for line in ll:
	line=line.rstrip('\n')
	line=line.split(',')
	gene=line[8]
	if gene in auxotrophy.keys():
		data=''
		i=1
		for thing in line:
			data+=thing
			if(i<len(line)):
				data+=','
				i+=1
		
		i=1
		data+=','
		for thing in auxotrophy[gene]:
			data+=thing
			if(i<len(auxotrophy[gene])):
				data+=','
				i+=1
		data=data+'\n'
		a.write(data)
		 
f.close()
a.close()