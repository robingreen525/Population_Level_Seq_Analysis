#!/bin/sh

/fh/fast/shou_w/bin/Robin/freebayes/freebayes -b ../alignment/aligned_reads.sorted.bam -f ../reference/sacCer3.fa -P 0.999 -p 20  -K >sample.freebyase.pop.pval99.popploidy20.vcf