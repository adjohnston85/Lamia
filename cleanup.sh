#!/bin/bash

rm *.fq.gz
rm $1_s.b*
rm $1.bam
mkdir stats
mkdir ./stats/trim_galore
mv *trimming_report* ./stats/trim_galore/
mv $1*txt ./stats/
mkdir ./stats/methyldackyl
mv $1*svg ./stats/methyldackyl/
mkdir ./stats/picard
mv $1*metrics.test ./stats/picard/
mkdir QC
mv *fastqc* ./QC/
echo 'Cleanup finished'
