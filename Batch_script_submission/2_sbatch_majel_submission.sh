#!/bin/bash
#SBATCH --time=04-00
#SBATCH --mem=60gb
#SBATCH --ntasks-per-node=20
#SBATCH --output=sbatch.out
#SBATCH --mail-type=ALL 

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

HELP='false'

RSYNC_TIME='08:00:00'
RSYNC_MEM='4gb'
SYNC_TO='/datasets/work/hb-meth-atlas/work/Data/level_2/public'

GENOME='hg38'
GENOME_PATH='/datasets/work/hb-meth-atlas/work/pipeline_data/Genomes/'
ALIGNER_THREADS='6'

#fetches the directory from which this script is located and run
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
#the Majel.py script is located one directory up from to bash scripts directory
MAJEL_DIR="$(dirname "$SCRIPT_DIR")"


#function checks if mandatory arguments have been set
check_argument() {
    if [ -z "$2" ];then
        printf '%s\n\n' "Error: --${1}= argument not set"
        exit 1
   fi
   printf '%s\n' "--${1}=$2"
}


while [ $# -gt 0 ]; do
  case "$1" in
    --project-dir=*)
      PROJECT_DIR="${1#*=}"
      ;;
    --sample-name=*)
      SAMPLE_NAME="${1#*=}"
      ;;
    --mail-user=*)
      EMAIL="${1#*=}"
      ;;
    --rsync-time=*)
      RSYNC_TIME="${1#*=}"
      ;;
    --rsync-mem=*)
      RSYNC_MEM="${1#*=}"
      ;;
    --sync-to=*)
      SYNC_TO="${1#*=}"
      ;;
    --genome=*)
      GENOME="${1#*=}"
      ;;
    --genome-path=*)
      GENOME_PATH="${1#*=}"
      ;;
    --aligner-threads=*)
      ALIGNER_THREADS="${1#*=}"
      ;;
    --script-dir=*)
      SCRIPT_DIR="${1#*=}"
      ;;
    --majel-dir=*)
      MAJEL_DIR="${1#*=}"
      ;;
    --majel-args=*)
      if [[ $1 != "--majel-args=" ]]; then
          MAJEL_ARGS="${1#*=} "
      fi
      ;;
    --help*)
      HELP="${1#*p}"
      ;;
    *)
  esac
  shift
done

if [ -z $HELP ]
then
    printf '\n'
    printf '%s\n' 'usage: sbatch_majel_submission_AJ.sh'
    printf '\n'
    printf '%s\n' '  --project-dir=         sets the path to the project directory containing the sample directory (e.g. --project-dir=/scratch1/usr001/PRJNA123456)'
    printf '%s\n' '  --sample-name=         sets name of the sample to run through Majel.py pipeline (e.g. --sample-name=Tissue_Subtissue_CancerType_SampleInfo_SAMN12345678)'
    printf '%s\n' '                         the sample name must correspond to a directory in the project directory and conform to the Majel.py naming conventions'
    printf '%s\n' '                         the sample directory must contain a data/ directory containing either SRA (.sra) or FASTQ (.fq or .fastq) files'
    printf '%s\n' '                         i.e. /path/to/PROJECT_NAME/SAMPLE_NAME/data/file.sra'
    printf '%s\n' '  --mail-user=           sets email for sbatch notifications from 3_sbatch_io_SyncProcessedData.sh'
    printf '\n'
    printf '%s\n' 'Optional arguments:'
    printf '%s\n' '  --rsync-time=          sets time allocated to sbatch_io_SyncProcessedData_AJ.sh (default: --rsync-time=08:00:00)'
    printf '%s\n' '  --rsync-mem=           sets memory allocated to sbatch_io_SyncProcessedData_AJ.sh (default: --rsync-mem=4gb'
    printf '%s\n' '  --sync-to=             sets path to directory being synced to (Default: --sync-to=/datasets/work/hb-meth-atlas/work/Data/level_2/public)'
    printf '%s\n' '  --genome=              used to alter --genome argument for Majel.py (e.g. --majel-genome=hg19)'
    printf '%s\n' '                         default: hg38'
    printf '%s\n' '  --genome-path=         used to alter --genome_path argument for Majel.py (e.g. --majel-genome-path=/path/to/Genomes/)'
    printf '%s\n' '                         default: /datasets/work/hb-meth-atlas/work/pipeline_data/Genomes/'
    printf '%s\n' '  --aligner-threads=     used to alter --aligner_threads for Majel.py (e.g. --majel-threads=4)'
    printf '%s\n' '                         default: 6'
    printf '%s\n' '  --majel-dir=           used to alter path to Majel.py'
    printf '%s\n' '                         default: /datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main'
    printf '%s\n' '  --majel-args=          used to add additional arguments to Majel.py (e.g. --majel-args="--pbat --is_paired_end False"'
    printf '\n'
    exit 1
fi

check_argument "project-dir" $PROJECT_DIR
check_argument "sample-name" $SAMPLE_NAME

if hash slurm 2> /dev/null; then
    check_argument "mail-user" $EMAIL
    RSYNC_PREFIX="sbatch --time=$RSYNC_TIME --mem=$RSYNC_MEM --mail-user=$EMAIL --job-name=RSYNC:$SAMPLE_NAME "
fi

LOG_FILE=$PROJECT_DIR/$SAMPLE_NAME/2_sbatch_majel_submission.log
> $LOG_FILE

cd $PROJECT_DIR/$SAMPLE_NAME

SCRIPT_DIR="$MAJEL_DIR/Batch_script_submission"

SUBMISSION="python3 $MAJEL_DIR/Majel.py --data_dir $PROJECT_DIR/$SAMPLE_NAME/data/ --sample_name $SAMPLE_NAME \
--genome $GENOME --genome_path $GENOME_PATH --aligner_threads $ALIGNER_THREADS ${MAJEL_ARGS}\
-v 3 -L $PROJECT_DIR/$SAMPLE_NAME/${SAMPLE_NAME}_majel.log &> slurm_majel_stdout.log"
TIME=$(date '+%B %d %T %Z %Y')
printf '%s\n\n' "$TIME> $SUBMISSION" | tee -a $LOG_FILE
eval $SUBMISSION

if grep -q "Completed Task = 'methylseekrAndTDF'" slurm_majel_stdout.log
then
  SUBMISSION="$MAJEL_DIR/majel_cleanup.sh $SAMPLE_NAME &>> slurm_majel_stdout.log"
  TIME=$(date '+%B %d %T %Z %Y')
  printf '%s\n\n' "$TIME> $SUBMISSION" | tee -a $LOG_FILE
  eval $SUBMISSION

  SUBMISSION="${RSYNC_PREFIX}$SCRIPT_DIR/3_sbatch_io_SyncProcessedData.sh --sync-to=$SYNC_TO --sync-from=$PROJECT_DIR/$SAMPLE_NAME &>> slurm_majel_stdout.log"
  TIME=$(date '+%B %d %T %Z %Y')
  printf '%s\n\n' "$TIME> $SUBMISSION" | tee -a $LOG_FILE
  eval $SUBMISSION

else
  printf '%s\n\n' "methylseekrAndTDF did not complete - Check for failed tasks" | tee -a $LOG_FILE
fi
