#!/bin/bash
#SBATCH --time=01-00
#SBATCH --ntasks-per-node=20
#SBATCH --output=sbatch.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrew.johnston@csiro.au
#SBATCH --job-name=SRR9070198_SRA_dl
#SBATCH --partition io

module load sratoolkit/2.8.2
module load parallel/20190722

help='false'

MAJEL_TIME="04-00"
MAJEL_MEM="60gb"
MAJEL_NTASKS="20"

RSYNC_TIME="08:00:00"
RSYNC_MEM="4gb"

GENOME='hg38'
GENOME_PATH='/datasets/work/hb-meth-atlas/work/pipeline_data/Genomes/'
ALIGNER_THREADS='6'
SCRIPT_DIR='/datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main/Batch_script_submission'
MAJEL_DIR='/datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main'
MAJEL_ARGS=''

#function checks if mandatory arguments have been set
check_argument() {
    if [ -z "$2" ]
    then
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
      if [ -z $SRA_ARRAY ]
      then
          SRA_ARRAY="${1#*=} "
      else
          SRA_ARRAY+="${1#*=} "
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
      MAJEL_ARGS="majel-args=${1#*=} "
      ;;
    --help*)
      help="${1#*p}"
      ;;
    *)
  esac
  shift
done

if [ -z $help ]
then
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
    printf '%s\n' '  --script-dir=          used to alter path to the Batch_script_submission/ directory'
    printf '%s\n' '                         default: /datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main/Batch_script_submission'
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

SRA_ARRAY=($SRA_ARRAY)

check_argument "majel-time" $MAJEL_TIME
check_argument "majel-ntaskts" $MAJEL_NTASKS
check_argument "majel-mem" $MAJEL_MEM
check_argument "rsync-time" $RSYNC_TIME
check_argument "rsync-mem" $RSYNC_MEM
check_argument "aligner-threads" $ALIGNER_THREADS
check_argument "genome" $GENOME
check_argument "genome-path" $GENOME_PATH
check_argument "script-dir" $SCRIPT_DIR
check_argument "majel-dir" $MAJEL_DIR

DATA_DIR="$PROJECT_DIR/$SAMPLE_NAME/data"

if [ ! -d $DATA_DIR ]
then
    mkdir -p $DATA_DIR
fi

cd $DATA_DIR

printf '%s\n' "${SRA_ARRAY[@]}" | parallel -j20 'eval "wget -O {}'.sra' $(srapath {})"'

eval "$SCRIPT_DIR/initilize_majel_submission.sh --skip-prompt --project-dir=$PROJECT_DIR --sample-name=$SAMPLE_NAME --mail-user=$EMAIL \
--majel_time=$MAJEL_TIME --majel-ntasks=$MAJEL_NTASKS --majel-mem=$MAJEL_MEM --rsync-time=$RSYNC_TIME --rsync-mem=$RSYNC_MEM \
--script-dir=$SCRIPT_DIR --majel-dir=$MAJEL_DIR --genome $GENOME --genome_path $GENOME_PATH --aligner_threads $ALIGNER_THREADS ${MAJEL_ARGS}"
