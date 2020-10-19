#!/bin/bash
#SBATCH --time=03-00
#SBATCH --mem=60gb
#SBATCH --ntasks-per-node=20
#SBATCH --output=sbatch.out

module load bowtie/2.2.9
module load fastqc/0.11.5
module load bismark/0.18.1
module load trimgalore/0.4.4
module load sratoolkit/2.8.2
module load samtools/1.5
module load picard/2.18.11
module load igvtools/2.4.14
module load R/3.6.1
module load bedtools/2.26.0
module load methyldackel/0.4.0
module load python/3.7.2

SCRIPT_DIR="/datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main/Batch_script_submission/"
MAJEL_DIR="/datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main/"
BATCH_SCRIPT_DIR="/datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main/Batch_script_submission/"

python3 $MAJEL_DIR/Majel.py --data_dir $2/$1/data/ --genome hg38 --sample_name $1 --genome_path /datasets/work/hb-meth-atlas/work/pipeline_data/Genomes/ -v 3 -L $2/$1/$1_majel.log --aligner_threads 6 &> slurm_majel_stdout.log
if grep -q "Completed Task = 'methylseekrAndTDF'" slurm_majel_stdout.log
then 
  $MAJEL_DIR/majel_cleanup.sh $1 &> slurm_cleanup_stdout.log
  rm -r ./data
  rm *fastq.gz
  sbatch $SCRIPT_DIR/sbatch_io_SyncProcessedData_flex.sh $1 $2 $3
fi
