


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