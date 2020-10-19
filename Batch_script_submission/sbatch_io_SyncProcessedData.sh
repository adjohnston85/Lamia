#!/bin/bash
#SBATCH --partition=io
#SBATCH --mem=4G
#SBATCH --time 08:00:00
#SBATCH --output=MappedData_rsync.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=warwick.locke@csiro.au

rsync -avzh --whole-file --remove-source-files /scratch1/loc100/wgbs/$1/ /datasets/work/hb-meth-atlas/work/Data/level_2/public/cao2020/$1 &> SyncProcessedData_stdout.log
