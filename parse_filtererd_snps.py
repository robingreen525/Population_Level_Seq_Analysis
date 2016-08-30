# the purpose of this script is to parse out useful information from the filtered vcf files
# I want the following information


#usage:
#python Population_Level_Seq_Analysis/parse_filtererd_snps.py -strainfile ../../Parameter+Strain_Files/AaronMiller_strains.txt -dir /fh/fast/shou_w/NextGenSeq/PopSeq/MDunham/Miller_SulfateGlucosePhosphate_2015/
#note: I'm using a symbolic link to my scripts directory (Github repo)


#For each snp in each sample, I want the following information:

#1) Chromosome
#2) Position on chromosome
#3) Reference allele
#4) Alternative allele
#5) Reads supporting Alternative allele
#6) Reads supporting reference allele
		
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

#for each strain, generate vcf.gz and index file for vcf-isec command in next loop
for strain in strains:
	
	full_path=path+strain
	
	vcf = full_path+'/'+strain+'.filtered.vcf'
	parsed_snps=full_path+'/'+strain+'.parsed.txt'
	
	f=open(vcf)
	d=open(parsed_snps,'w')
	
	ll=f.readlines()
	for line in ll:
		if(line[0])!='#': #only care about lines with snp information.
			line=line.split()
			if(len(line))>1:
				info=line[7] ## has the AO and RO information I want
				chr=line[0] #chromosome
				pos=line[1] # position
				ref_allele=line[3] # ref allele
				alt_allele=line[4]# alt allele
				
				info=info.split(';') # split the info string by ';' where I can access each piece of information seperatly
				for thing in info:
					thing=thing.split('=')
					if(thing[0] == 'RO' or thing[0] == 'AO'):
						if(thing[0]=='RO'):ref_obv = (thing[1])
						elif(thing[0]=='AO'):alt_obv = (thing[1])
				
				str= chr+','+pos+','+ref_allele+','+ref_obv+','+alt_allele+','+alt_obv+'\n'
				d.write(str) # writes information to 'sample.parsed.txt'
	f.close()
	d.close()			
	