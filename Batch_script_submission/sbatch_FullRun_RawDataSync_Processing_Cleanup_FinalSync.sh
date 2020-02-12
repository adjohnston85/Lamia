#!/bin/bash
#SBATCH --partition=io
#SBATCH --mem=4G
#SBATCH --time 24:00:00
#SBATCH --output=io_batch.out

SCRIPT_DIR="/datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/develop/majel_wgbspipline/Batch_script_submission/"
MAJEL_DIR="/datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/develop/majel_wgbspipline/"
SAMPLENAME=$(basename $1)
mkdir -p $2/$SAMPLENAME/data
cd $2/$SAMPLENAME/

dmget -a --recurse $1
RC=1
while [[ $RC -ne 0 ]]
do
  rsync -avzh --no-g $1/ $2/$SAMPLENAME/data
  RC=$? && ((RC!=0)) && sleep 60
  ((c++)) && ((c>10)) && break && echo "Max Iterations"
done

sbatch $SCRIPT_DIR/sbatch_majel_submission.sh $SAMPLENAME $2
  
