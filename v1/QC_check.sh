#the purpose of this script is to run FastQC on all read files being analyzed in the 'compare' diretory

ls |while read i; # iterate in all files
	do
		if [ -d $i ];
			then
			cd $i;
			dirc="reads"; 
			b=$(pwd);
			# if the 'reads' directory is present, go into directory and do FastQC check. If not, move on.
			if [ -d "$dirc" ]; #requires space before/after each bracket because shell is stupid
				then
					cd $dirc;
					pwd;
					rm -rf QC_report/
					mkdir QC_report/
					fastqc *.fastq
					mv *fastqc* QC_report/
					#ls;
					cd ../
					
				fi
			cd ../
			fi
	done;
	

