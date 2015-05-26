# the purpose of this script is to determine if the observed SNPS will result in an amino acid substitution.
# i have a fasta file that gives the sequence for each gene in the cerviesea genome (from SGD). I will build a dictionary
# of the genes and use that information to assess frame and mutations for each SNP called.

#usage
#python Population_Level_Seq_Analysis/record_aminoacid_changes.py -genefasta /fh/fast/shou_w/NextGenSeq/reference/S288C_orf_genomic_all.fasta -strainfile /fh/fast/shou_w/NextGenSeq/PopSeq/Parameter+Strain_Files/AaronMiller_strains.txt -dir /fh/fast/shou_w/NextGenSeq/PopSeq/MDunham/Miller_SulfateGlucosePhosphate_2015/ -genefile /fh/fast/shou_w/NextGenSeq/reference/S288C.gff



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



rev_comp={'a':'t','t':'a','c':'g','g':'c'}

def rev_comp_codon(codon,dict):
	#new=''
	#for base in codon:
		#new=new+dict[base]
	
	#return(new[::-1])
	return(codon)

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
	
	
# create master.csv file where all mutation information for every strain will be stored	
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

					if(strand=='+'): # annotation for Waston strand
						frame=((pos-start+1)%3) # 1=first base, 2=second base, 0= third base
						
						if(frame==1): # first base in codon
							codon1=seq[(start-start):((pos-start)+3)] # get 2 bases after  base 
							codon=codon1[-3:] # last 3 bases in strong
							newcodon=list(codon)
							newcodon[0]=line[4] #change first base in codon to SNP
							newcodon=''.join(newcodon)
							codon=codon.lower()
							newcodon=newcodon.lower()
							
							if(len(codon)==3 and len(newcodon)==3): #would miss indels, fix later
								aa_ref=code[codon] # get AA of ref
								aa_new=code[newcodon] # get AA of alt
															
								if(aa_ref!=aa_new): #if not a silent mutation, record 
									data=strain+','	
									i=1
									for thing in line:
										data=data+thing.strip('\n')
										if i<len(line):
											data=data+','
											i+=1								
									data=data+','+aa_ref+str(len(codon1)/3)+aa_new+','+codon+','+newcodon+','+str(frame)+'\n'
									a.write(data)
									
								
						elif(frame==2): # if mutation in second position
							codon1=seq[(start-start):((pos-start)+2)] # get one base after snp base
							codon=codon1[-3:]
							newcodon=list(codon)
							newcodon[1]=line[4]
							newcodon=''.join(newcodon)
							codon=codon.lower()
							newcodon=newcodon.lower()
							
							if(len(codon)==3 and len(newcodon)==3): #would miss indels, fix later
								aa_ref=code[codon]
								aa_new=code[newcodon]
							
								if(aa_ref!=aa_new):
									data=strain+','
									
									i=1
									for thing in line:
										data=data+thing.strip('\n')
										if i<len(line):
											data=data+','
											i+=1
									
									data=data+','+aa_ref+str(len(codon1)/3)+aa_new+','+codon+','+newcodon+','+str(frame)+'\n'
									a.write(data)
						elif(frame==0): # if snp in 3ed base
							codon1=seq[(start-start):(pos-start)+1] # last base is 3rd position in codon
							codon=codon1[-3:]
							newcodon=list(codon)
							
							newcodon[-1]=line[4]
							newcodon=''.join(newcodon)
							
							codon=codon.lower()
							newcodon=newcodon.lower()
							
							if(len(codon)==3 and len(newcodon)==3): #would miss indels, fix later
								codon=str(rev_comp_codon(codon,rev_comp)) 
								newcodon=str(rev_comp_codon(codon,rev_comp))
								aa_ref=code[codon]
								aa_new=code[newcodon]
							
								if(aa_ref!=aa_new):
									data=strain+','
									
									i=1
									for thing in line:
										data=data+thing.strip('\n')
										if i<len(line):
											data=data+','
											i+=1
									
									
									data=data+','+aa_ref+str(len(codon1)/3)+aa_new+','+codon+','+newcodon+','+str(frame)+'\n'
									a.write(data)	

					elif(strand=='-'):
						# if on crick strand, must swap end and start positions. Also do rev comp of bases
						temp=end
						end=start
						start=temp
						
						frame=((start-pos)+1)%3# 1=first base, 2=second base, 0= third base
						
						if(frame==1):
							codon1=seq[(start-start):((start-pos)+3)]
							codon=codon1[-3:]

							if(len(line[4])==1 and len(codon)==3): #would miss indels, fix later
								codon=codon.lower()
								codon=str(rev_comp_codon(codon,rev_comp))
								newcodon=list(codon)
								
								newcodon[0]=rev_comp[line[4].lower()]
								newcodon=''.join(newcodon)
								newcodon=newcodon.lower()		
								aa_ref=code[codon]
								newcodon=str(rev_comp_codon(newcodon,rev_comp))
								aa_new=code[newcodon]
								
								if(aa_ref!=aa_new):
									data=strain+','
									
									i=1
									for thing in line:
										data=data+thing.strip('\n')
										if i<len(line):
											data=data+','
											i+=1
									
									data=data+','+aa_ref+str(len(codon1)/3)+aa_new+','+codon+','+newcodon+','+str(frame)+'\n'	
									a.write(data)

						elif(frame==2):
							codon1=seq[(start-start):((start-pos)+2)]
							codon=codon1[-3:]

							if(len(line[4])==1 and len(codon)==3): #would miss indels, fix later
								codon=codon.lower()
								codon=str(rev_comp_codon(codon,rev_comp))
								newcodon=list(codon)
								newcodon[1]=rev_comp[line[4].lower()]
								newcodon=''.join(newcodon)
								newcodon=newcodon.lower()		
								aa_ref=code[codon]
								newcodon=str(rev_comp_codon(newcodon,rev_comp))
								aa_new=code[newcodon]
								
								if(aa_ref!=aa_new):
									data=strain+','
									
									i=1
									for thing in line:
										data=data+thing.strip('\n')
										if i<len(line):
											data=data+','
											i+=1
									
									data=data+','+aa_ref+str(len(codon1)/3)+aa_new+','+codon+','+newcodon+','+str(frame)+'\n'
									a.write(data)
						elif(frame==0):
							codon1=seq[(start-start):((start-pos))+1]
							codon=codon1[-3:]
							newcodon=list(codon)

							if(len(line[4])==1 and len(codon)==3): #would miss indels, fix later
								codon=codon.lower()
								codon=str(rev_comp_codon(codon,rev_comp))
								newcodon=list(codon)
								newcodon[-1]=rev_comp[line[4].lower()]
								newcodon=''.join(newcodon)
								newcodon=newcodon.lower()		
								aa_ref=code[codon]
								newcodon=str(rev_comp_codon(newcodon,rev_comp))
								aa_new=code[newcodon]
							
								if(aa_ref!=aa_new):
									
									data=strain+','
									
									i=1
									for thing in line:
										data=data+thing.strip('\n')
										if i<len(line):
											data=data+','
											i+=1
									
									
									data=data+','+aa_ref+str(len(codon1)/3)+aa_new+','+codon+','+newcodon+','+str(frame)+'\n'
									a.write(data)	
a.close()							
		