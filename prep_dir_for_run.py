#!/bin/sh


import sys
import os

f=open(sys.argv[1],'r')

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

	#os.system(w)

	os.chdir('../')