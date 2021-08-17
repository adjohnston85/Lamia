#!/bin/bash

HELP='false'
SKIP_PROMPT='false'

####Default settings####

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
    if [ -z $2 ]
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
    --skip-prompt*)
      SKIP_PROMPT="${1#*t}"
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
    printf '%s\n' 'usage: initilize_majel_submission.sh [--help] [--project-dir=<path>] [--sample-name=<name>] [--majel-time=<time>] [--majel-ntasks=<ntasks>] [--majel-mem=<size[units]> [--rsync-time=<time>] [--rsync-mem=<size[units]>'
    printf '\n'
    printf '%s\n' 'Mandatory arguments:'
    printf '%s\n' '  --project-dir=         sets the path to the project directory containing the sample directory (e.g. --project-dir=/scratch1/usr001/PRJNA123456)'
    printf '%s\n' '  --sample-name=         sets name of the sample to run through Majel.py pipeline (e.g. --sample-name=Tissue_Subtissue_CancerType_SampleInfo_SAMN12345678)'
    printf '%s\n' '                         the sample name must correspond to a directory in the project directory and conform to the Majel.py naming conventions'
    printf '%s\n' '                         the sample directory must contain a data/ directory containing either SRA (.sra) or FASTQ (.fq or .fastq) files'
    printf '%s\n' '                         i.e. /path/to/PROJECT_NAME/SAMPLE_NAME/data/file.sra'
    printf '%s\n' '  --mail-user=           sets email for SLURM notifications'
    printf '\n'
    printf '%s\n' 'Optional arguments:'
    printf '%s\n' '  --majel-time=          sets --time= allocated to sbatch_majel_submission_AJ.sh     (default: --majel-time=04-00)'
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
    printf '%s\n' '  --skip-prompt          skips verification step - use when user input is not possible (i.e. when submitting script with sbatch)'
    printf '\n'
    exit 1
fi

check_argument "project-dir" $PROJECT_DIR
check_argument "sample-name" $SAMPLE_NAME
check_argument "mail-user" $EMAIL

if [ ! -d "$PROJECT_DIR/$SAMPLE_NAME/data" ] 
then
    echo "Error: Directory $PROJECT_DIR/$SAMPLE_NAME/data DOES NOT exist" 
    exit 1
fi

count=`find $PROJECT_DIR/$SAMPLE_NAME/data -type f \( -iname \*.fastq.gz -o -iname \*.sra -o -iname \*.fq.gz \) | wc -l`
if [ $count = 0 ]
then
    echo "Error: No SRA (.sra) or FASTQ (.fastq.gz OR fq.gz) files exist in $PROJECT_DIR/$SAMPLE_NAME/data"
    exit 1
fi

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

printf '\n'
printf '%s\n' 'These parameters will result in the following Majel.py job:'
printf '\n'

echo "python3 $MAJEL_DIR/Majel.py --data_dir $PROJECT_DIR/$SAMPLE_NAME/data/ --sample_name $SAMPLE_NAME \
--genome $GENOME --genome_path $GENOME_PATH --aligner_threads $ALIGNER_THREADS ${MAJEL_ARGS} \
-v 3 -L $PROJECT_DIR/$SAMPLE_NAME/${SAMPLE_NAME}_majel.log &> slurm_majel_stdout.log"
printf '\n' 

if [[ $SKIP_PROMPT == "false" ]]
then
    read -p "Would you like to continue?(Y/N)" -n 1 -r
    printf '\n'
fi

cd $PROJECT_DIR/$SAMPLE_NAME

if [[ $REPLY =~ ^[Yy]$ ]] || [ -z $SKIP_PROMPT ]
then
    SUBMISSION="sbatch --time=$MAJEL_TIME --ntasks-per-node=$MAJEL_NTASKS --mem=$MAJEL_MEM --job-name=${SAMPLE_NAME}_majel \
--mail-user=$EMAIL $SCRIPT_DIR/sbatch_majel_submission_AJ.sh --sample-name=$SAMPLE_NAME --project-dir=$PROJECT_DIR \
--rsync-time=$RSYNC_TIME --rsync-mem=$RSYNC_MEM --mail-user=$EMAIL --genome=$GENOME --genome-path=$GENOME_PATH \
--script-dir=$SCRIPT_DIR --majel-dir=$MAJEL_DIR --majel-args='${MAJEL_ARGS}'"

    echo "The following job was submitted: 
$SUBMISSION" | tee initilize_majel.log

eval $SUBMISSION

else
    echo 'job was not submitted to sbatch' | tee initilize_majel.log
fi
