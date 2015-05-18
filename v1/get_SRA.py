#!/bin/sh


import sys
import os

f=open(sys.argv[1],'r')

ll=f.readlines()

for line in ll:
	line=line.rstrip('\n')
	line=line.rstrip('^M')
	dir_exist=os.path.isdir(line) ## for each SRA sample, check if directory already exists, if not, then make one
	print line
	line1=line
	if(dir_exist==False):
		
		line='mkdir '+line
		os.system(line)

	
	os.chdir(line1)
	os.system('pwd')

	##import reads from SRA archive

	main=line1[:6]
	#print(main)
	#test='ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR129/SRR1299103/SRR1299103.sra'
	q='ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/'+main+'/'+line+'/'+line1+'.sra'
	w='wget '+q

	os.system(w)
	






	os.chdir('../')
