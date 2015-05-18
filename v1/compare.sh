#!/bin/sh
export PERL5LIB=/fh/fast/shou_w/bin/Robin/vcftools/vcftools_0.1.11/perl/
#bash QC_check.sh 


user_dir=/fh/fast/shou_w/Robin_Code/full_pipeline/populations

cd $user_dir
cd ../
cd $user_dir
pwd

mv cn.mops/ ../temp/
ls| while read i; do
 if [ -d $i ];
 then
   cd $i;
   pwd;
   rm slurm*;
   sbatch -n 1 -c 1 run.sh;
   cd ..;
 fi
done  
mv ../temp cn.mops/

ls
cd $user_dir
echo 'Delay.....'
sleep 8h

cd  $user_dir


echo 'Calling SNV and small indels'
cd ancestor/GATK/;
/fh/fast/shou_w/bin/Robin/tabix/bgzip round1.vcf;
/fh/fast/shou_w/bin/Robin/tabix/tabix -p round1.vcf.gz;
cd ../../

anc_snp=$user_dir/ancestor/GATK/round1.vcf.gz

#head $anc_snp

echo $anc_snp
mv cn.mops/ ../temp
ls| while read i; do
 if [ -d $i ];
 then
     cd $user_dir/$i;
     cd GATK;
     pwd;
     gunzip *;
     /fh/fast/shou_w/bin/Robin/tabix/bgzip round1.vcf;
    # tabix -p round1.vcf.gz;
     /fh/fast/shou_w/bin/Robin/tabix/tabix -p vcf  round1.vcf.gz; 
     /fh/fast/shou_w/bin/Robin/vcftools/vcftools_0.1.11/bin/vcf-isec -c round1.vcf.gz $anc_snp > evo_unqiue.vcf
     /fh/fast/shou_w/bin/Robin/vcftools/vcftools_0.1.11/bin/vcf-isec -c $anc_snp round1.vcf.gz  > reversion_to_WT.vcf  
     
     cd ../../;
  fi
done

mv ../temp cn.mops/
echo 'Calling CNV'
bash build_cn.mops.sh
echo 'Done'
