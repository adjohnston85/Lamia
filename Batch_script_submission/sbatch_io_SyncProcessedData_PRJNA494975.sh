#!/bin/bash
#SBATCH --partition=io
#SBATCH --mem=4G
#SBATCH --time 08:00:00
#SBATCH --output=MappedData_rsync.out

rsync -avzh --whole-file --remove-source-files /scratch1/loc100/$1/ /datasets/work/hb-meth-atlas/work/Data/level_2/hg38/public/PRJNA494975/$1 &> SyncProcessedData_stdout.log
