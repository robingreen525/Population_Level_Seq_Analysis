#!/usr/bin/perl
# Requirements:
# Have X11 installed
# log into lynx with: ssh -X -u username@lynx.fhcrc.org
# Required programs: samtools, vcftools, GATK, picard, freebayes
# Currently installed versions
# samtools: 1.0
# VCFtools (v0.1.12b)
# GATK 3.2-2-gec30cee (also needs Java v1.7)
# Picard 1.102
 
# freebayes
#Robin's command
#	freebayes -b aligned_reads.sorted.bam -f saccharomyces_cerevisiae_rm11-1a_1_supercontigs.fasta -P 0.95 -K  -F 0.05 --standard-filters  >sample.freebyase.pop.pval99.vcf
# OPTIONAL Improvements: 
#change so read files don't have to be there?cd

###RJG EDIT LOG 2015-05-014
#Picard mark duplicate seems to be failing with large bam file from popseq data. Added -Xmx4g, which seems to use smaller allocations and memory. Suggest at
# http://sourceforge.net/p/samtools/mailman/samtools-help/thread/4DE11062.6000606@ukmuenster.de/

#### RJG EDIT LOG 2015-05-15


use strict; use warnings;  use File::Copy; use Getopt::Std;

# check that the -s and -p flags were set
our($opt_s, $opt_p, $opt_n); getopts('ns:p:');
unless($opt_p){die "\n\nUSAGE: PopSeqPipeline.pl -s Strains.txt -p params.txt\n\n"}
unless($opt_s){die "\n\nUSAGE: PopSeqPipeline.pl -s Strains.txt -p params.txt\n\n"}

# Setting up timestamp and runtime start.
my $totalstart=time; my $strainstart; my $functionstart;
my@time=localtime();
my$stamp= ($time[5]+1900).sprintf("%02d",$time[4]+1).sprintf("%02d",$time[3]).sprintf("%02d",$time[2]).sprintf("%02d",$time[1]);

# saving original STDOUT & STDERR
open(SAVEDSTDOUT, ">&STDOUT") || warn "Can't save STDOUT\n"; #saves STDOUT so it can later be restored
open(SAVEDSTDERR, ">&STDERR") || warn "Can't save STDERR\n"; #saves STDERR so it can later be restored

##########################################			
#  Read/Check Parameter file from $opt_p #
##########################################
open PARAMS, "<$opt_p" or die $!;
my$straindir; # directory path containing the folders for each strain
my@straindirs;  #subdirectories found in #straindir
my$refdir; # directory path where reference sequences are found
my$jardir; # directory path where java jars are found
my$bindir; # directory path where scripts arefound
my$logdir; #directory path where logs are saved
my%params; # program specific parameter settings
my%overwrite; #overwrite existing files settings 
my %rmint; # remove intermediate files settings
my $paramprint=""; #parameter print messsage
$paramprint.="\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\nRunning the following programs:\n\n";

while(<PARAMS>){
	chomp;
	s/""/ /g; s/"//g; 	#remove quotes around fields in param file, replace empty quotes with single space
	if(/^#/){next}elsif(/^\s*$/){next}		#skip commented and empty lines
	elsif(/^strainDir\t(\S+)$/){$straindir=$1; $paramprint.="\tstrains are in folder $straindir\n"}	#saving strain directory
	elsif(/^refDir\t(\S+)$/){$refdir=$1; $paramprint.="\treference genomes are in folder $refdir\n"}	#saving strain directory
	elsif(/^jarDir\t(\S+)$/){$jardir=$1; $paramprint.="\tjava jar files are in folder $jardir\n"}		#saving reference genome directory
	elsif(/^binDir\t(\S+)$/){$bindir=$1; $paramprint.="\texecutable files are in folder $bindir\n\n"}		#saving reference genome directory
	elsif(/^(\S+)\t([01])\t([^\t]+)\t([01]+)\t([01])$/){		# finding programs, parameters, and settings
		if($2==1){	# if set to run (1)
			if(defined$3){
				$params{$1}=$3; #saving parameters
				$overwrite{$1}=$4; #saving overwrite setting
				$rmint{$1}=$5; # saving remove intermediate files setting
			}
			$paramprint.= "\t$1 has been set to RUN w/ options: $params{$1}\toverwrite:$overwrite{$1}\trmint:$rmint{$1}\n";
		}
		else{
			$params{$1}="dontrun";
			$paramprint.= "\t$1 has been set to NOT RUN \n";
		}
	}
	else{$paramprint.= "\nMalformated line in $opt_p:\n$_\n\n"}; # line in parameter file didn't match any of the above formats
}
print "$paramprint\n";
close PARAMS or die $!;
die "\nERROR:\t\t Reference directory $refdir doesn't exist\n" unless -e $refdir; #stop if reference directory can't be found
die "\nERROR:\t\t BIN directory $refdir doesn't exist\n" unless -e $bindir; #stop if binary directory can't be found
die "\nERROR:\t\t JAR directory $jardir doesn't exist\n" unless -e $jardir; #stop if reference directory can't be found
die "\nERROR:\t\t Strain directory $straindir doesn't exist\n" unless -e $straindir;  #stop if directory of strains can't be found
##########################################			
#  Prepare the Reference Directory	 #
##########################################
print"\nIndexing Reference Files\n";
opendir REFS, $refdir or die "Cannot open $!";
foreach(readdir REFS){	
	if(/\.fasta$/||/\.fa$/){ #for every .fasta or .fa file do the following 4 commands (indexing & Dictionary building)
		system "$bindir"."samtools faidx $refdir$_" unless -s "$refdir$_.fai";
		system "$bindir"."bwa index $refdir$_" unless -s "$refdir$_.ann";
		system "$bindir"."java -jar $jardir"."CreateSequenceDictionary.jar R=$refdir$_ O=$refdir$_.dict" unless -s "$refdir$_.dict";
		system "mrsfast --index $refdir$_" unless -s "$refdir$_.index";
	}
}	
##########################################			
### Read/Check Strain file from $opt_s ###
##########################################

opendir STRAINDIR, $straindir or die "Cannot open $!"; #open $straindir	
foreach (readdir STRAINDIR){push @straindirs, $_} 
closedir STRAINDIR;

open STRAINS, "<$opt_s" or die $!; #read in the file specified after -s
my@strains; my%refs;
while(<STRAINS>){
	chomp;
	if(/^#/){next}elsif(/^\s*$/){next}		#skip commented and empty lines
	elsif(/^(\S+)\t(\S+)\t?\S+?$/){		#get the strain and reference
		unless(-e $straindir.$1){warn "\nERROR:\t\tCouldn't find any directory named $1 in $straindir. Not running line $_ \n";next;} #check strain exists
		unless(-s $refdir.$2.".fasta"){print "\nERROR:\t\tCouldn't find reference genome $2.fasta in $refdir. Not running line $_\n"; next;} # check reference exists
	
		unless(grep(/^$1$/,@strains)){push @strains,$1} # if not yet in strains list, add it;
		$refs{$1}=$2;  #save its reference genome
	}
	else{print "\nSkipping malformated line in $opt_s:\n$_\n\n"} # report lines that aren't in the right format
}
close STRAINS or die $!;

print "\nRunning Analysis on the following strains:\n\n";
foreach my$strain (@strains){
	print "\tStrain $strain with reference $refs{$strain}\n";
}
print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";

##########################################################			
## Analysis of strains start with this loop ###
##########################################################
my$straincount=0;
foreach my$strain (@strains){ #for each strain in the list

	$strainstart=time;  $straincount++; #for progress reporting 
	print "\n#######	Starting analysis of strain $strain with reference $refs{$strain} #######\n\n";
	my$reffile = $refdir.$refs{$strain}.".fasta";	#make path to reference genome 
	my$subdir=$straindir.$strain."\/";	#path to strain folder

	unless( -e $subdir."VS_Ref"){ system "mkdir $subdir"."VS_Ref"} #summary folder for comparison to reference
	unless( -e $subdir."logs"){ system "mkdir $subdir"."logs"} #log folder for this strain
	
	unless($opt_n){
		unless( -e $subdir."logs"){ system "mkdir $subdir"."logs"} #log folder for this strain
		open STDERR, ">>", $subdir."logs/$strain"."_$stamp.log" or die "Can't open logfile.txt: $!\n";  #all STDERR is redirected to the log file, only things printed to SAVEDSTDERR are printed to the screen
		open STDOUT, ">>", $subdir."logs/$strain"."_$stamp.log" or die "Can't open logfile.txt: $!\n";  #all STDOUT is redirected to the log file, only things printed to SAVEDSTDOUT are printed to the screen
		print "$paramprint\n";
	}

##########################################			
## GZIP & CONCATENATION OF FASTQ FILES ###
##########################################

	opendir STRAINDIR, $subdir or die "Cannot open $!"; #open $straindir	
	
	my @readfiles;
	unless($params{gzip} eq "dontrun"){
		$functionstart=time;							#used to determine how long it takes for the function to run
		foreach(sort (readdir STRAINDIR)){				#for each file in $straindir
			if(/^\Q$strain\E_\w+?_L\d+(_R[12])_\d+(\.fastq)(\.gz)?$/){	# check if it starts with $strain and ends in .fastq or .fastq.gz 
				if(-s "$subdir$strain$1$2.gz"){ #if it has already been concatenated and zipped
						print  SAVEDSTDOUT  "GZIP		Looks like $subdir$strain$1$2.gz already exists. Skipping gzip & concatenation of read files\n";
						next;
				}
				else{	#if it hasn't
					push @readfiles, $subdir.$_; 
					#print SAVEDSTDOUT "\nRunning GZIP on $subdir.$_\n";
				}
			}
		}
		foreach my$readfile (@readfiles){ #for each file added to @readfiles combines the fastq files for the forward (R1) and reverse (R2) reads
			if($readfile =~ /^\Q$subdir$strain\E_\S+?(_R[12])_\d+(\.fastq(\.gz)?)$/){
				if(defined$3){										# if already a .gz file
					print  "GZIP		Already zipped, Executing: cat $readfile >> $subdir$strain$1$2\n";
					system("cat $readfile >>$subdir$strain$1$2");			#append / initiate .gz file
				} 	
				else{					# otherwise gzip then concatenate
					print  "GZIP		Not yet zipped, Executing: gzip $subdir$readfile then cat $readfile.gz >> $subdir$strain$1$2.gz\n";
					system("gzip $readfile"); #zip it 
					system("cat $readfile.gz >> $subdir$strain$1$2.gz"); #append / initiate .gz file
				}	
			}			
		}
		print SAVEDSTDOUT "GZIP on $strain COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n GZIP runtime: ".((time-$functionstart)/60)." min\n\n";
	}
	closedir STRAINDIR;			

### finding names of the R*.fastq.gz files
	my@strainfiles;
	opendir STRAINDIR, $subdir or die "Cannot open $!"; #open $straindir	
	foreach (readdir STRAINDIR){
		if(/^(\Q$strain\E_R[12])\.fastq\.gz$/){push @strainfiles, $1}
	}	
	closedir STRAINDIR;
	
	foreach my $file (@strainfiles){			#start of foreach @strainfiles
##########################################			
####			FASTQC				  ####
##########################################	
		unless($params{fastqc} eq "dontrun"){
			$functionstart=time;
			my $infile="$subdir$file.fastq.gz"; #specify the file to be read
		
			if((-s "$subdir$file"."_fastqc.html")&&($overwrite{fastqc} == 0)){print SAVEDSTDOUT "FASTQC		Looks like fastQC was already run for $file.fastq.gz Skipping...\n"} #check if it's already been run and isn't set to overwrite
			else{
				system "rm $subdir$file"."_fastqc.html" if -e "$subdir$file"."_fastqc.html"; #remove any preexisting files used for later success check
				print SAVEDSTDOUT "\nRunning FASTQC on $strain\n";
				my @exec =("$bindir"."fastqc", $infile); #build command: fastq
				if($overwrite{fastqc}==1){print "FASTQC		Overwriting any existing files \n"}
				print"FASTQC		Executing: @exec\n";
				system(@exec); #run command
				@exec =("cp","$subdir$file"."_fastqc.html",$subdir."VS_Ref"); #build command: cp to summary folder
				system(@exec); #run command
			}
			print SAVEDSTDERR "\nERROR		FASTQC on $file FAILED!!!!\n\n" unless -s "$subdir$file"."_fastqc.html";
			print SAVEDSTDOUT "FASTQC on $file COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n FASTQC runtime: ".((time-$functionstart)/60)." min\n\n";
		}
		
##########################################			
## bwa | samtools view | samtools sort  ##
##########################################	
# ../../../bin/bwa mem -M -R '@RG\tID:Sample_WS\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' ../../reference/RM11.fasta ACY45_R1.fastq.gz ACY45_R2.fastq.gz | ../../../bin/samtools view  -Su - | samtools sort - ACY45_sorted

		if($file eq $strainfiles[-1]){ #if $file is the last file in $strainfiles
			$file =~ s/_R[12]$//; # drop R1/R2 from end of file name
		
			unless($params{align} eq "dontrun"){
				$functionstart=time;
				my $infile="$subdir$file"."_R1.fastq.gz"; #file to be read
				if(-s "$subdir$file"."_R2.fastq.gz"){$infile.=" $subdir$file"."_R2.fastq.gz"}			#if paired end run include read2 file
				my $outfile="$subdir$file.sorted"; #file to be generated
				my $params = $params{align}; # get parameters for ALIGN
		
				if((-s "$outfile.mdup.bam")&&($overwrite{align} == 0)){print SAVEDSTDOUT "ALIGN		Looks $outfile,mdup.bam already exists. Skipping....\n"}	
				else{
					system "rm $outfile.mdup.bam" if -e "$outfile.mdup.bam"; #remove any preexisting files used for later success check
					print SAVEDSTDOUT "\nRunning ALIGN on $strain\n";
					
					# Build command bwa | samtools view | samtools sort
					my @exec= ("$bindir"."bwa", $params{align}, $reffile, $infile, "|");
					push @exec, ("$bindir"."samtools view", "-Su","-","|");
					push @exec, ("$bindir"."samtools sort", "-", $outfile);
		
					if($overwrite{align}==1){print "ALIGN		Overwriting any existing files\n"}
					print "\nALIGN		Executing: @exec \n\n";	
					system ("@exec"); #run bwa | samtools view | samtools sort
				
					# Build picard markduplicates
					#@exec = ("$bindir"."java -Xmx4g -jar $jardir"."MarkDuplicates.jar","I=$outfile.bam","O=$outfile.mdup.bam","M=$outfile.mdup.metric");
				
					#print "\nALIGN		Executing: @exec \n\n";
					#system ("@exec"); #run picard markduplicates
					
					
					##############
					#use samtools rmdup instead of picard
					
					# Build samtools rmdup
					@exec = ("$bindir"."samtools "."rmdup "."$outfile.bam "."$outfile.mdup.bam");
				
					print "\nALIGN		Executing: @exec \n\n";
					system ("@exec"); #run samtools rmdup
					
					
					
				
					#samtools index
					@exec = ("$bindir"."samtools index", "$outfile.mdup.bam");
				
					print "\nALIGN		Executing: @exec \n\n";
					system ("@exec"); # run samtools index
				
					if($rmint{align}==1){
						print "\nALIGN		Removing intermediate file: $outfile.bam\n\n";
						system "rm $outfile.bam";
					}
				}
				print SAVEDSTDOUT "\nERROR		ALIGN on $strain FAILED!!!!\n\n" unless -s $outfile.".mdup.bam"; #success check
				print SAVEDSTDOUT "ALIGN on $strain COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n ALIGN runtime: ".((time-$functionstart)/60)." min\n\n";		
			} #end of ALIGN
####################################################		
######	freebayes							 #######
####################################################	
#freebayes -C 5 --pooled-discrete -p1 --standard-filters -b bamfile -f reference.fasta > output.vcf
			unless($params{freebayes} eq "dontrun"){
				$functionstart=time;
			
				my $infile="$subdir$file.sorted.mdup.bam";
				my $outfile="$subdir$file.freebayes.vcf";
				my $params = $params{freebayes};

				if((-s $outfile)&&($overwrite{freebayes} == 0)){print SAVEDSTDOUT "FREEBAYES		Looks $outfile already exists. Skipping...\n"}
				else{
					system "rm $outfile" if -e $outfile; #remove any preexisting files used for later success check
					print SAVEDSTDOUT "\nRunning FREEBAYES on $strain\n";
					if($overwrite{freebayes}==1){print "FREEBAYES  	Overwriting any existing files \n"}
					
					my @exec = ($bindir."freebayes",$params,"-b $infile","-f $reffile", ">$outfile"); 
					print"\nFREEBAYES 		Executing: @exec \n\n";
					system("@exec");		
						
				}
				print SAVEDSTDOUT "\nERROR		FREEBAYES on $strain FAILED!!!!\n\n" unless -s $outfile; #success check
				print SAVEDSTDOUT "FREEBAYES on $strain COMPLETE\n Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n FREEBAYES runtime: ".((time-$functionstart)/60)." min\n\n";
			}# end of freebayes	
#######################################			
######	End of tools		  #####
######################################			
									
		}#end of if last strainfile
	}#end of foreach @strainfiles
	open(STDOUT, ">&SAVEDSTDOUT") || warn "Can't restore STDOUT\n";
	open(STDOUT, ">&SAVEDSTDERR") || warn "Can't restore STDERR\n";
	print "\n$strain COMPLETE!! ($straincount of ".(scalar@strains).")\n  Total runtime: ".((time-$totalstart)/60)." min\n Strain runtime: ".((time-$strainstart)/60)." min\n";
}#foreach @strains
print "\nRUN COMPLETE!!\n  Total runtime: ".((time-$totalstart)/60)." min\n\n";
__END__
