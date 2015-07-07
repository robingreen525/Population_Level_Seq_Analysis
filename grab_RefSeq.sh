#!/bin/sh

python get_refseq_yeast.py
grep '>' S.cerevisea_RefSeq.fasta  > S.cerevisea_RefSeqOnly.fasta 
rm S.cerevisea_RefSeq.fasta 
