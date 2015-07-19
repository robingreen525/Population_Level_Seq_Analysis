
#usage
#python Population_Level_Seq_Analysis/assign_RefSeq_to_auxotrophfile.py 
#-RefSeqfile S.cerevisea_RefSeq_header.fasta 
#-dir /fh/fast/shou_w/NextGenSeq/PopSeq/MDunham/Miller_SulfateGlucosePhosphate_2015/


#import libraries needed
import sys
import os
import argparse

#define parsing options
parser = argparse.ArgumentParser()
parser.add_argument('-RefSeqfile', type=argparse.FileType('r'))
parser.add_argument('-dir', dest='dir')


args = parser.parse_args()

# get the path where the data are stored. will be used with other strings later
path=args.dir

auxofile=path+'/auxotrophy_snps.csv'
auxofile=open(auxofile,'r')
auxofile=auxofile.readlines()

ll=args.RefSeqfile.readlines()


refseq={}

for line in ll:
	line=line.rstrip('\n')
	line=line.split('|')
	gene=line[-1].split()
	#print gene
	for thing in gene:
		letter=thing[0]
		#print thing
		if letter.isupper() and thing !='S288c]':
			print line[1],thing
			if line[1] not in refseq.keys():
				if(len(thing)==4): ## try to ensure that protein names are chosen
					refseq[thing]=line[1]
		



ll=auxofile
f=open('auxotrophy_snps_giNumber.csv','w')

for line in ll:
	
	line=line.rstrip('\n')
	line_copy=line
	line = line.split(',')

	if line[9]!='gene':
		if(line[9] in refseq.keys()):
			print line[9],refseq[line[9]],line[15]
			line_copy=line_copy+','+refseq[line[9]]+'\n'
			f.write(line_copy)
			
		else:
			#print line[9],'########',line[15]
			line_copy=line_copy+','+'##########'+'\n'
			f.write(line_copy)
	else: # writes header information to new file
		header=line_copy+',gi_Number'+'\n'
		f.write(header)
			
f.close()

### want to build file to submit to SIFT
unique={}
f=open('auxotrophy_snps_giNumber.csv','r')
s=open('MillerSNPsAuxoSIFT_Batch1.csv','w')
ll=f.readlines()
for line in ll:
	line=line.rstrip('\n')
	line=line.split(',')
	if(line[-1][0]!='#'):
		#print line[11],line[-1]
		key=line[-1]+'_'+line[11]
		if key not in unique.keys():
			unique[key]=1
			add=line[-1]+','+line[11]+'\n'
			s.write(add)
			

f.close()
s.close()
			



