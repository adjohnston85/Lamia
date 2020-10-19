#!/bin/bash
#SBATCH --partition=io
#SBATCH --mem=4G
#SBATCH --time 08:00:00
#SBATCH --output=MappedData_rsync.out

rsync -avzh --no-g --no-o --whole-file --remove-source-files $2/$1/ $3/$1 &> SyncProcessedData_stdout.log
