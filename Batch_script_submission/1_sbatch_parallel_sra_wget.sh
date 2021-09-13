#!/bin/bash
#SBATCH --time=01-00
#SBATCH --ntasks-per-node=1
#SBATCH --output=sbatch.out
#SBATCH --mail-type=ALL
#SBATCH --partition io

module load sratoolkit/2.8.2
module load parallel/20190722

HELP='false'

MAJEL_TIME="04-00"
MAJEL_MEM="60gb"
MAJEL_NTASKS="1"

RSYNC_TIME="08:00:00"
RSYNC_MEM="4gb"

GENOME='hg38'
GENOME_PATH='/datasets/work/hb-meth-atlas/work/pipeline_data/Genomes/'
ALIGNER_THREADS='6'
MAJEL_DIR='/datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main'

#function checks if mandatory arguments have been set
check_argument() {
    if [ -z "$2" ];then
        echo "Error: --${1}= argument not set"
        exit 1
   fi
   echo "--${1}=$2"
}

#grab input argument values
while [ $# -gt 0 ]; do
  case "$1" in
    --project-dir=*)
      PROJECT_DIR="${1#*=}"
      ;;
    --sample-name=*)
      SAMPLE_NAME="${1#*=}"
      ;;
    --sra-array=*)
      if [[ -z $SRA_ARRAY ]]; then
          SRA_ARRAY="${1#*=}"
      else
          SRA_ARRAY+=" ${1#*=}"
      fi
      ;;
    --majel-time=*)
      MAJEL_TIME="${1#*=}"
      ;;
    --majel-mem=*)
      MAJEL_MEM="${1#*=}"
      ;;
    --majel-ntasks=*)
      MAJEL_NTASKS="${1#*=}"
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
      MAJEL_ARGS="${1#*=}"
      ;;
    --help*)
      help="${1#*p}"
      ;;
    *)
  esac
  shift
done

if [ -z $HELP ];then
    printf '\n'
    printf '%s\n' 'usage: initilize_majel_submission.sh [--help] [--project-dir=<path>] [--sample-name=<name>] [--mail-user=<email>] [sra-array=<list>]'
    printf '\n'
    printf '%s\n' 'Mandatory arguments:'
    printf '%s\n' '  --project-dir=         sets the path to the project directory containing the sample directory (e.g. --project-dir=/scratch1/usr001/PRJNA123456)'
    printf '%s\n' '  --sample-name=         sets name of the sample to run through Majel.py pipeline (e.g. --sample-name=Tissue_Subtissue_CancerType_SampleInfo_SAMN12345678)'
    printf '%s\n' '                         the sample name must correspond to a directory in the project directory and conform to the Majel.py naming conventions'
    printf '%s\n' '                         the sample directory must contain a data/ directory containing either SRA (.sra) or FASTQ (.fq or .fastq) files'
    printf '%s\n' '                         i.e. /path/to/PROJECT_NAME/SAMPLE_NAME/data/file.sra'
    printf '%s\n' '  --mail-user=           sets email for SLURM notifications'
    printf '%s\n' '  --sra-array=           sets the list of SRAs to be downloaded (e.g. --SRA-array="SRR1234567 SRR1234568" OR --SRA-array=SRR123456{7..8} )'
    printf '%s\n' '                         Note: as per the example a list of SRAs must be contained within quotation marks'
    printf '\n'
    printf '%s\n' 'Optional arguments:'
    printf '%s\n' '  --majel_time=          sets --time= allocated to sbatch_majel_submission_AJ.sh     (default: --majel-time=04-00)'
    printf '%s\n' '  --majel-ntasks=        sets --ntasks-per-node= for sbatch_majel_submission_AJ.sh   (default: --majel-ntasks=20)'
    printf '%s\n' '  --majel-mem=           sets --mem= allocated to sbatch_majel_submission_AJ.sh      (default: --majel-mem=60gb'
    printf '%s\n' '  --rsync-time=          sets time allocated to sbatch_io_SyncProcessedData_AJ.sh    (default: --rsync-time=08:00:00)'
    printf '%s\n' '  --rsync-mem=           sets memory allocated to sbatch_io_SyncProcessedData_AJ.sh  (default: --rsync-mem=4gb'
    printf '%s\n' '  --genome=              used to alter --genome argument for Majel.py (e.g. --majel-genome=hg19)'
    printf '%s\n' '                         default: hg38'
    printf '%s\n' '  --genome-path=         used to alter --genome_path argument for Majel.py (e.g. --majel-genome-path=/path/to/Genomes/)'
    printf '%s\n' '                         default: /datasets/work/hb-meth-atlas/work/pipeline_data/Genomes/'
    printf '%s\n' '  --aligner-threads=     used to alter --aligner_threads for Majel.py (e.g. --majel-threads=4)'
    printf '%s\n' '                         default: 6'
    printf '%s\n' '  --majel-dir=           used to alter path to Majel.py'
    printf '%s\n' '                         default: /datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main'
    printf '%s\n' '  --majel-args=          used to add additional arguments to Majel.py (e.g. --majel-args="--pbat --is_paired_end False"'
    printf '%s\n' '                         Note: this list of additional arguments must be contained within quotation marks'
    printf '\n'
    exit 1
fi

check_argument "project-dir" $PROJECT_DIR
check_argument "sample-name" $SAMPLE_NAME
check_argument "mail-user" $EMAIL
check_argument "sra-array" "$SRA_ARRAY"

LOG_FILE=$PROJECT_DIR/$SAMPLE_NAME/slurm_submission_stdout.log
if [ ! -f $LOG_FILE ]; then
    > $LOG_FILE
fi

if [[ ! -z $MAJEL_ARGS ]]; then
    MAJEL_ARGS=$(echo "$MAJEL_ARGS " | tr -d '"')
fi

SCRIPT_DIR="$MAJEL_DIR/Batch_script_submission"

DATA_DIR="$PROJECT_DIR/$SAMPLE_NAME/data"
if [[ ! -d $DATA_DIR ]];then
    mkdir -p $DATA_DIR
fi

cd $DATA_DIR
SRA_ARRAY=($SRA_ARRAY)
printf '%s\n' "${SRA_ARRAY[@]}" | parallel -j20 'eval "wget -O {}.sra $(srapath {})"' &> $PROJECT_DIR/$SAMPLE_NAME/slurm_SRA_download.log

COUNT=0
for SRA in ${SRA_ARRAY[@]}; do
    if [[ -s $SRA.sra ]]; then
        echo "$SRA.sra exists and is not empty" | tee -a $LOG_FILE
    else
        echo "$SRA.sra doesn't exist or is empty" | tee -a $LOG_FILE
        COUNT=$((COUNT+1))
    fi
done

if [[ "$COUNT" -gt 0 ]]; then
    echo "Error: $COUNT SRA file(s) specified in --sra-array= doesn't exist or is empty. Exiting pipeline." | tee -a $LOG_FILE
    exit 1
fi

cd $PROJECT_DIR/$SAMPLE_NAME
SUBMISSION="sbatch --time=$MAJEL_TIME --ntasks-per-node=$MAJEL_NTASKS --mem=$MAJEL_MEM --job-name=MAJEL:${SAMPLE_NAME} \
--mail-user=$EMAIL $SCRIPT_DIR/2_sbatch_majel_submission.sh --sample-name=$SAMPLE_NAME --project-dir=$PROJECT_DIR \
--rsync-time=$RSYNC_TIME --rsync-mem=$RSYNC_MEM --mail-user=$EMAIL --genome=$GENOME --genome-path=$GENOME_PATH \
--majel-dir=$MAJEL_DIR --majel-args=\"${MAJEL_ARGS}\" &> slurm_majel_stdout.log"
TIME=$(date '+%B %d %T %Z %Y')
printf '%s\n\n' "$TIME> $SUBMISSION" | tee -a $LOG_FILE
eval $SUBMISSION 

