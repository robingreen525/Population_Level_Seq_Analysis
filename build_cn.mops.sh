#!/bin/sh

echo "Cleaning files"
rm -rf cn.mops/bam_files/*
ls |while read i;
	do
		if [ -d $i ];
			then
			cd $i;
			dirc="alignment";
			b=$(pwd);
			if [ -d "$dirc" ]; #requires space before/after each bracket because shell is stupid
				then
					cd $dirc;
					pwd;
					#ls;
					cp aligned_reads.sorted.bam ../../cn.mops/bam_files/$i
					cd ../
					
				fi
			cd ../
			fi
	done;
	

cd cn.mops;
sbatch -n 1 -c 1 run_cn.mops.sh

