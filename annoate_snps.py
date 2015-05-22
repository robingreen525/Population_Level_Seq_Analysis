# the purpose of this script is to annotate each SNP with the gene identifier (standard and systematic name)
# along with what codon and position in the codon. This information can be used to determine which frame 
# the mutation is in and predict amino acid substitutions  


#usage:
# python Population_Level_Seq_Analysis/annoate_snps.py -genefile /fh/fast/shou_w/NextGenSeq/reference/S288C.gff -strainfile /fh/fast/shou_w/NextGenSeq/PopSeq/Parameter+Strain_Files/AaronMiller_strains.txt -d /fh/fast/shou_w/NextGenSeq/PopSeq/MDunham/Miller_SulfateGlucosePhosphate_2015/
# works if in /fh/fast/shou_w/NextGenSeq/PopSeq/MDunham/Miller_SulfateGlucosePhosphate_2015



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

# get the path where the data are stored. will be used with other strings later
path=args.dir


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
				lookup_table.append([chr,start,end,strand,ID,'unknown']) # if no standard name, just add 'unknown'

# for each strain.parsed.txt file, annotate snp calls (gene, mutation, aa change?) and write to
# strain.parsed.annotated.txt file



for strain in strains:
	
	full_path=path+strain
	file=full_path+'/'+strain+'.parsed.txt'
	newfile=full_path+'/'+strain+'.annotated.parsed.txt'
	print(file)
	if(os.path.exists(file)):
		file=open(file,'r')
		newfile=open(newfile,'w')
		ll=file.readlines()
		for line in ll:  #for each snp, find where it is in the lookup_table
			line=line.split(',')
			if(len(line)>1):
				chr=line[0] #chromosome
				pos=int(line[1]) # position
				for gene in lookup_table:
					chr_ref=gene[0] #chromosome of each snp in lookup table
					start_pos_ref=int(gene[1]) # starting position (from the Waston strand (+))
					end_pos_ref=int(gene[2])  # ending position (from the Watson strand (+))
					if(chr_ref==chr):
						if(start_pos_ref<=pos):
							if(end_pos_ref>=pos): ## if snp found within gene
								ref_reads=float(line[3])
								alt_reads=float(line[len(line)-1].rstrip('\n'))
								freq=(float (alt_reads/(alt_reads+ref_reads)))
								string=line[0]+','+line[1]+','+line[2]+','+line[3]+','+line[4]+','+line[5].rstrip('\n')+','+gene[3]+','+gene[4]+','+gene[5]+','+str(freq)+'\n'
								newfile.write(string)
								#newfile.write(line[0],line[1],line[2],line[3],line[4],line[5].rstrip('\n'), gene[3],gene[4],gene[5],freq,'\n')
	
	file.close()
	newfile.close()							