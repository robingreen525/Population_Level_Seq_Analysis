#!/bin/sh


import sys
import os

f=open(sys.argv[1],'r')
d=open('Sherlock_strains.txt','w')


ll=f.readlines()

for line in ll:
	line=line.rstrip('\n')
	print line
	dir_exist=os.path.isdir(line) ## for each SRA sample, check if directory already exists, if not, then make one
	line1=line
	if(dir_exist==False):
		line='mkdir '+line
		os.system(line)

	
	os.chdir(line1)
	os.system('pwd')

	w='cp -s /fh/fast/shou_w/NextGenSeq/AllReads/Outside_Reads/GSherlock_Stanford/'+line1+'/*.fastq.gz ./'
	#w='cp -s /fh/fast/shou_w/NextGenSeq/AllReads/Outside_Reads/JBarrick_UTexas/'+line1+'/*.fastq.gz ./'
	os.system(w)

	w=line1+'\t'+'S288C\t'+'WY1335\n'
	d.write(w)
	
	
	os.chdir('../')


d.close()
f.close()
