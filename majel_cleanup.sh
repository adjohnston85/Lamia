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
rm -r ./tmp
mkdir browser_tracks
mv $1_sd_CpG.tdf ./browser_tracks/
mkdir MethylSeekR
if grep -q "MethylSeekR not run, coverage < 10x" $1_PMD.bed
then
  echo "MethylSeekR not run, inadequate coverage"
  touch ./MethylSeekR/InadequateCoverage.txt
  rm $1{_PMD.bed,_UMRLMR.bed,_wPMD_UMRLMR.bed}
else
  mv *UMR* ./MethylSeekR/
  mv *PMD* ./MethylSeekR/
  mv *CalculateFDR* ./MethylSeekR/
  mv *Segmentation* ./MethylSeekR/
  mv $1_AlphaDistribution.pdf ./MethylSeekR/
fi
echo 'Compressing all text files'
find . \( -name "*.txt" -o -name "*.bedGraph" -o -name "*.bed" \) -exec gzip {} \;
echo 'Cleanup finished'
