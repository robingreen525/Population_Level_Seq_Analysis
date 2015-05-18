#!/bin/sh

#cp -rf  ancestor/ CC109RS1/
#cp -rf  ancestor/ CC74RL2/
#cp -rf  ancestor/ CT61/

for line in $(cat samples.txt);
   do
    echo $line;
    #a=${line:5)
    p=$(echo $line | cut -d '_' -f2)
    echo $p
    #echo $s
    rm -rf $p/
    cp -rf ancestor/ $p/
    cp ../../../20140620_Seq/$line/* $p/reads/
    gunzip $p/reads/*

  done  


cp ../../../20140620_Seq/Sample_WY1335/* ancestor/reads/
gunzip ancestor/reads/*