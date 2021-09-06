#!/bin/bash

if [ -z $1 ]
then
    SAMPLE_NAME=$(basename $PWD)
else
    SAMPLE_NAME=$1
fi

rm -r ./data
rm *.fastq.gz
rm *.fq.gz
rm ${SAMPLE_NAME}_s.b*
rm ${SAMPLE_NAME}.bam
mkdir stats
mkdir ./stats/trim_galore
mv *trimming_report* ./stats/trim_galore/
mv ${SAMPLE_NAME}*txt ./stats/
mkdir ./stats/methyldackyl
mv ${SAMPLE_NAME}*svg ./stats/methyldackyl/
mkdir ./stats/picard
mv ${SAMPLE_NAME}*metrics.test ./stats/picard/
mkdir QC
mv *fastqc* ./QC/
rm -r ./tmp
mkdir browser_tracks
mv ${SAMPLE_NAME}_sd_CpG.tdf ./browser_tracks/
mkdir MethylSeekR
if grep -q "MethylSeekR not run, coverage < 10x" ${SAMPLE_NAME}_PMD.bed
then
  echo "MethylSeekR not run, inadequate coverage"
  touch ./MethylSeekR/InadequateCoverage.txt
  rm ${SAMPLE_NAME}_PMD.bed,_UMRLMR.bed,_wPMD_UMRLMR.bed
else
  mv *UMR* ./MethylSeekR/
  mv *PMD* ./MethylSeekR/
  mv *CalculateFDR* ./MethylSeekR/
  mv *Segmentation* ./MethylSeekR/
  mv ${SAMPLE_NAME}_AlphaDistribution.pdf ./MethylSeekR/
fi
echo 'Compressing all text files'
find . \( -name "*.txt" -o -name "*.bedGraph" -o -name "*.bed" \) -exec gzip {} \;
echo 'Cleanup finished'
