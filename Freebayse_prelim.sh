#!/bin/sh

/fh/fast/shou_w/bin/Robin/freebayes/freebayes -b ../alignment/aligned_reads.sorted.bam -f ../reference/saccharomyces_cerevisiae_rm11-1a_1_supercontigs.fasta -P 0.999 -p 20  -K >sample.freebyase.pop.pval99.popploidy20.vcf

python /fh/fast/shou_w/bin/Robin/user_scripts/population_level/parse_population_vcf.py sample.freebyase.pop.pval99.popploidy20.vcf >filtered_calls.txt



python /fh/fast/shou_w/bin/Robin/user_scripts/population_level/annotate_population_data.py /fh/fast/shou_w/bin/Robin/user_scripts/population_level/saccharomyces_cerevisiae_rm11-1a_1_genome_summary_per_gene.txt ~/shougroup/_DATABASE/_Sequences_AND_Strain_Keys/reference_genomes/RM11-1a/saccharomyces_cerevisiae_rm11-1a_1_transcripts.gtf filtered_calls.txt /fh/fast/shou_w/bin/Robin/user_scripts/population_level/auxotrophy_annotations.txt >annotated_call.txt