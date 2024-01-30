#!/bin/bash
#$ -cwd
#$ -S /bin/bash


# basic argument
Aligned_Path=$1
INPUT=$2
ID=${INPUT%%.bam}

samtools view -b -h -q 30 $Aligned_Path$INPUT > $Aligned_Path${ID}.MAPQ_30.bam
samtools sort $Aligned_Path${ID}.MAPQ_30.bam
samtools index $Aligned_Path${ID}.MAPQ_30.bam 
