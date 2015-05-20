
# the purpose of this script is to run vcf-isec of all Aaron Miller population level data against the FY4
# ancestor strain. This will generate SNPs that are unique to the population.

# usage:
#python Population_Level_Seq_Analysis/vcf_isec_all.py -strainfile ../../Parameter+Strain_Files/AaronMiller_strains.txt -dir /fh/fast/shou_w/NextGenSeq/PopSeq/MDunham/Miller_SulfateGlucosePhosphate_2015/ -anc FY4


#import libraries needed
import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-strainfile', type=argparse.FileType('r'))
parser.add_argument('-dir', dest='dir')
parser.add_argument('-anc',dest='anc')
                      
args = parser.parse_args()

# get the path where the data are stored. will be used with other strings later
#get the ancestral sample that all SNPs will be filtered against
path=args.dir
anc=args.anc




strains=[]
### use the strain file flag to get the list of strains to analyze
ll=args.strainfile.readlines()
for line in ll:
	line=line.split('\t')
	strains.append(line[0])
	



#for each strain, generate vcf.gz and index file for vcf-isec command in next loop
for strain in strains:

	full_path=path+strain
	anc_path=path+anc
	
	vcf_file=full_path+'/'+strain+'.freebayes.vcf'
	
	if(os.path.exists(vcf_file)): # if vcf file exist not bgzipped, bgzip it
		cmd='bgzip '+vcf_file
		print(cmd)
		os.system(cmd)
		
	bg_file=full_path+'/'+strain+'.freebayes.vcf.gz'
	if(os.path.exists(bg_file)): # if bgzip file exist, index
		cmd='tabix -f -p vcf '+bg_file
		print(cmd)
		os.system(cmd)
		
		
#vcf-isec -c A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz	
# for each indexed vcf.gz file, filter against anc.vcf.gz	
for strain in strains:

	full_path=path+strain
	anc_path=path+anc
	
	bg_file=full_path+'/'+strain+'.freebayes.vcf.gz'
	anc_bg_file=anc_path+'/'+anc+'.freebayes.vcf.gz'

	
	cmd='vcf-isec -c '+bg_file+' '+anc_bg_file+' > '+full_path+'/'+strain+'.filtered.vcf'
	print(cmd)
	os.system(cmd)
	
	
	# i get a broken pipe write error when i run this script to completion, but the filtered vcf file is generated.
	# I will ignore this issue for now
