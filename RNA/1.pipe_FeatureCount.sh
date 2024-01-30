#!/bin/bash
#$ -cwd
#$ -S /bin/bash


# basic argument
DIR=$1
InputPath=$2
Input=$3
OutputPath=$4
Sample=$5
GTF=$DIR/BED/gencode.v38lift37.annotation.gtf

/home/jen96/tools/subread_featurecount/bin/featureCounts \
-T 2 \
-a $GTF \
-s 2 \
-Q 10 \
-p \
-g gene_id \
-t exon \
-o $OutputPath/${Sample}.RGadded.marked.FeatureCounts_reverse_pairend_min10.txt \
$InputPath$Input


