
#import libraries needed
import sys
import os
import argparse

#define parsing options
parser = argparse.ArgumentParser()
parser.add_argument('-addfile',type=argparse.FileType('r'))
parser.add_argument('-genefile',type=argparse.FileType('r'))
parser.add_argument('-dir', dest='dir')


args = parser.parse_args()

# get the path where the data are stored. will be used with other strings later
path=args.dir

ll=args.genefile.readlines()

for line in ll:
	line=line.rstrip('\n')
	line=line.split('\t')
	if len(line)>8:
		line=line[8]
		line = line.split(';')
		if len(line)>3:
			for thing in line:
				thing=thing.split('=')
				if thing[0]=='Alias':
					print thing