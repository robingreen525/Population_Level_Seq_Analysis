# the purpose of this script is to determine if the observed SNPS will result in an amino acid substitution.
# i have a fasta file that gives the sequence for each gene in the cerviesea genome (from SGD). I will build a dictionary
# of the genes and use that information to assess frame and mutations for each SNP called.



# Dictionary holding the code to life: Codon to Amino Acid map
code = {     'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
             'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
             'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
             'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
             'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
             'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
             'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
             'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
             'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
             'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
             'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
             'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
             'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
             'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
             'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
             'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
        }

#import libraries needed
import sys
import os
import argparse

#define parsing options
parser = argparse.ArgumentParser()
parser.add_argument('-genefasta',type=argparse.FileType('r'))
parser.add_argument('-strainfile', type=argparse.FileType('r'))
parser.add_argument('-dir', dest='dir')
parser.add_argument('-genefile',type=argparse.FileType('r'))

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





gene_seq={}

### build a dictionary of genes and their sequences
ll=args.genefasta.readlines()
end=len(ll)
i=0

while i<end-1:
	line=ll[i]
	if line[0]=='>': # for each header, get the corresponding seqeunce for each gene
		line=line.split()
		sysname=line[0][1:] # get the systematic name
		geneseq=''
		i+=1
		newline=ll[i]
		while newline[0]!='>' and i<end-1: # keep going until you reach another header line
			geneseq=geneseq+newline.rstrip('\n')
			i+=1
			newline=ll[i]
		gene_seq[sysname]=geneseq # assign seq to key (systematic name)
	
	
all_mutations=path+'/master.csv'	
a=open(all_mutations,'w')	
	
for strain in strains:
	full_path=path+strain
	file=full_path+'/'+strain+'.annotated.parsed.txt'
	
	print(file)
	if os.path.exists(file):
		f=open(file)
		ll=f.readlines() # open sample.annotated.parsed.txt file
		for line in ll:
			line=line.split(',')
			sys_gene=line[7] # get systamtic name for gene where snp is
			seq=gene_seq[sys_gene] # get the sequence of the gene
			for gene in lookup_table: # find the start,end, and strand  for the gene being analyzed
				sys_gene_ref=gene[4]
				if sys_gene_ref==sys_gene: # if the gene being analyzed is found, get positional information
					start=int(gene[1])
					end=int(gene[2])
					strand=gene[3]
					pos=int(line[1])
					
					#print sys_gene,sys_gene_ref,start, end, strand
					#print(strand)
					if(strand=='+'):
						frame=((pos-start)%3) # 1=first base, 2=second base, 0= third base
						
						if(frame==1):
							codon=seq[(start-start):(pos-start)+3]
							codon=codon[-3:]
							newcodon=list(codon)
							newcodon[0]=line[4]
							newcodon=''.join(newcodon)
							#print(line)
							#print codon
							#print newcodon
							codon=codon.lower()
							newcodon=newcodon.lower()
							
							if(len(codon)==3 and len(newcodon)==3): #would miss indels, fix later
								aa_ref=code[codon]
								aa_new=code[newcodon]
							
								if(aa_ref!=aa_new):
									#print 'Nonsynamous mutation'
									#print aa_ref
									#print aa_new
									data=strain+','
									
									for thing in line:
										data=data+thing.strip('\n')+','
									
									data=data+','+aa_ref+str(len(seq)/3)+aa_new+'\n'
									#print(data)
									a.write(data)
									
								
						elif(frame==2):
							#codon=seq[(start-start):(pos-start)+3]
							codon=seq[(start-start):(pos-start)+2]
							codon=codon[-3:]
							newcodon=list(codon)
							newcodon[0]=line[4]
							newcodon=''.join(newcodon)
							#print(line)
							#print codon
							#print newcodon
							codon=codon.lower()
							newcodon=newcodon.lower()
							
							if(len(codon)==3 and len(newcodon)==3): #would miss indels, fix later
								aa_ref=code[codon]
								aa_new=code[newcodon]
							
								if(aa_ref!=aa_new):
									#print 'Nonsynamous mutation'
									#print aa_ref
									#print aa_new
									data=strain+','
									
									for thing in line:
										data=data+','+thing.strip('\n')
									
									data=data+','+aa_ref+str(len(seq)/3)+aa_new+'\n'
									#print(data)
									a.write(data)
						elif(frame==0):
							#codon=seq[(start-start):(pos-start)+3]
							codon=seq[(start-start):(pos-start)+1]
							codon=codon[-3:]
							newcodon=list(codon)
							newcodon[0]=line[4]
							newcodon=''.join(newcodon)
							#print(line)
							#print codon
							#print newcodon
							codon=codon.lower()
							newcodon=newcodon.lower()
							
							if(len(codon)==3 and len(newcodon)==3): #would miss indels, fix later
								aa_ref=code[codon]
								aa_new=code[newcodon]
							
								if(aa_ref!=aa_new):
									#print 'Nonsynamous mutation'
									#print aa_ref
									#print aa_new
									data=strain+','
									
									for thing in line:
										data=data+thing.strip('\n')+','
									
									data=data+','+aa_ref+str(len(seq)/3)+aa_new+'\n'
									#print(data)
									a.write(data)	

a.close()							
		