#!/bin/bash
#SBATCH --partition=io
#SBATCH --mem=4G
#SBATCH --time 24:00:00
#SBATCH --output=io_batch.out

SCRIPT_DIR=$(dirname $(readlink -f $0))
MAJEL_DIR=$(dirname $SCRIPT_DIR)
SAMPLENAME=$(basename $1)
mkdir -p $2/$SAMPLENAME/sra
dmget -a --recurse $1
RC=1
while [[ $RC -ne 0 ]]
do
  rsync -avzh --no-g $1/ $2/$SAMPLENAME/sra
  RC=$? && ((RC!=0)) && sleep 60
  ((c++)) && ((c>10)) && break && echo "Max Iterations"
done


  
