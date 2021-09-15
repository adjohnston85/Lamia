#!/bin/bash

HELP='false'
SKIP_PROMPT='false'

#function checks if arguments have been set, exit if values are empty
check_argument() {
    if [[ -z "$2" ]]; then
        echo "Error: --${1}= argument not set"
        exit 1
   fi
   echo "--${1}=$2"
}


#function to set or reset default values
set_defaults() {
    
    #by default this script can be run from a SAMPLE directory within a PROJECT directory ( e.g. /path/to/PROJECT/SAMPLE )
    PROJECT_DIR=$(dirname $PWD)
    SAMPLE_NAME=$(basename $PWD)
    CONCURRENT_JOBS=5

    MAJEL_TIME='04-00'
    MAJEL_MEM='60gb'
    MAJEL_NTASKS='20'

    RSYNC_TIME='08:00:00'
    RSYNC_MEM='4gb'

    GENOME='hg38'
    GENOME_PATH='/datasets/work/hb-meth-atlas/work/pipeline_data/Genomes/'
    ALIGNER_THREADS='6'

    #fetches the directory from which this script is located and run
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
    #the Majel.py script is located one directory up from to bash scripts directory
    MAJEL_DIR="$(dirname "$SCRIPT_DIR")"

    unset SRA_ARRAY
    unset MAJEL_ARGS
}


#function checks if arguments have been set and prints their values
check_argument() {
    if [[ -z "$2" ]]; then
        echo "Error: --${1}= argument not set" 
        exit 1
   fi
   echo "--${1}=$2" | tee -a $LOG_FILE
}


#function to grab input argument values
get_arguments() {
    while [[ $# -gt 0 ]]; do
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
        --sample-file=*)
          SAMPLE_FILE="${1#*=}"
          ;;
        --concurrent-jobs=*)
          CONCURRENT_JOBS="${1#*=}"
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
}


job_submission() {

    LOG_FILE=$PROJECT_DIR/$SAMPLE_NAME/slurm_submission_stdout.log

    #check the values mandatory arguments 
    check_argument "project-dir" $PROJECT_DIR
    check_argument "sample-name" $SAMPLE_NAME

    #print info on the sample file used to run multiple Majel jobs at once
    if [[ ! -z $SAMPLE_FILE ]]; then
        echo "--sample-file=$SAMPLE_FILE"
        echo "--concurrent-jobs=$CONCURRENT_JOBS"
    fi
    
    #check if this pipeline is being run on cluster with a slurm submission system, if so specifying an email is mandatory
    if command -v slurm; then
        check_argument "mail-user" $EMAIL
    fi
        
    #check if mandatory /data directory exists in path, if not make it
    DATA_DIR="$PROJECT_DIR/$SAMPLE_NAME/data"
    if [[ ! -d $DATA_DIR ]]; then
        mkdir -p $DATA_DIR
    fi

    #covert path variable to full path (i.e. ./destionation -> /full/path/to/destionation)
    PROJECT_DIR="$(cd $PROJECT_DIR && pwd)"

    #if SRAs were not specified for download make sure sequence files exist in /data directory
    if [[ -z $SRA_ARRAY ]]; then
        COUNT=`find $PROJECT_DIR/$SAMPLE_NAME/data -type f \( -iname \*.fastq.gz -o -iname \*.sra -o -iname \*.fq.gz \) | wc -l`
        if [ $COUNT = 0 ]; then
            echo "Error: $PROJECT_DIR/$SAMPLE_NAME/data must contain SRA (.sra) or FASTQ (.fastq.gz OR fq.gz) files if --sra-array= is not used to specify SRA file(s) for download. Exiting pipeline."
            exit 1
        fi
    else
        if [[ $SRA_ARRAY != *"}"* ]]; then
            #removes potential " arifacts that can arise when --sra-array= is specified in the SAMPLE_FILE as a string series
            SRA_ARRAY="\"$(echo $SRA_ARRAY | tr -d '"')\""
        else
            #converts a specified range ( e.g. --sra-array=SRR1{1..3} ) into a string series ( e.g. --sra-array="SRR1 SRR2 SRR3" )
            SRA_ARRAY=$SRA_ARRAY
        fi

        check_argument "sra-array" "$SRA_ARRAY"
    fi

    cd $PROJECT_DIR/$SAMPLE_NAME

    #print the values of all arguments for user to examine
    check_argument "majel-time" $MAJEL_TIME
    check_argument "majel-ntaskts" $MAJEL_NTASKS
    check_argument "majel-mem" $MAJEL_MEM
    check_argument "rsync-time" $RSYNC_TIME
    check_argument "rsync-mem" $RSYNC_MEM
    check_argument "aligner-threads" $ALIGNER_THREADS
    check_argument "genome" $GENOME
    check_argument "genome-path" $GENOME_PATH
    check_argument "majel-dir" $MAJEL_DIR
    
    #add space after additional majel arguments, if they exist
    if [[ ! -z $MAJEL_ARGS ]]; then
        MAJEL_ARGS=$(echo "$MAJEL_ARGS " | tr -d '"')
    fi
    
    #prints MAJEL_ARGS values without checking for existance (this argument is not required for subsequent steps)
    echo "--majel-args=\"${MAJEL_ARGS}\"" | tee -a $LOG_FILE

    #print time and script inputs
    TIME=$(date '+%B %d %T %Z %Y')
    printf '%s\n' "$TIME> $BASH_SOURCE --project-dir=$PROJECT_DIR --sample-name=$SAMPLE_NAME --mail-user=$EMAIL --majel-time=$MAJEL_TIME --majel-ntasks=$MAJEL_NTASKS \
--majel-mem=$MAJEL_MEM --rsync-tim=$RSYNC_TIME --rsync-mem=$RSYNC_MEM --aligner-threads=$ALIGNER_THREADS --genome=$GENOME --genome-path=$GENOME_PATH \
--majel-dir=$MAJEL_DIR --majel-args=\"${MAJEL_ARGS}\" --sra-array=\"$SRA_ARRAY\"" > $LOG_FILE

    if [[ -z $SRA_ARRAY ]]; then
        #set job submission variable for if --sra-array= was not declared (skips SRA download script)
        SUBMISSION="sbatch --time=$MAJEL_TIME --ntasks-per-node=$MAJEL_NTASKS --mem=$MAJEL_MEM --job-name=MAJEL:${SAMPLE_NAME} \
--mail-user=$EMAIL $SCRIPT_DIR/2_sbatch_majel_submission.sh --sample-name=$SAMPLE_NAME --project-dir=$PROJECT_DIR \
--rsync-time=$RSYNC_TIME --rsync-mem=$RSYNC_MEM --mail-user=$EMAIL --genome=$GENOME --genome-path=$GENOME_PATH \
--majel-dir=$MAJEL_DIR --majel-args=\"${MAJEL_ARGS}\" &>> $LOG_FILE"
    else
        #convert string series to proper array
        TMP_ARRAY=($SRA_ARRAY)
        LEN_ARRAY=$(expr ${#TMP_ARRAY[@]} - 1)
        #calculate number of cores based on number of SRA files to be downloaded (5 file downloads per core, as advised by IMT)
        CORES=$(($LEN_ARRAY/5+1))
 
        ##set job submission variable for if --sra-array= was declared
        SUBMISSION="sbatch --job-name=SRA_DL:${SAMPLE_NAME} --mail-user=$EMAIL --ntasks-per-node=$CORES $SCRIPT_DIR/1_sbatch_parallel_sra_wget.sh \
--majel-time=$MAJEL_TIME --majel-ntasks=$MAJEL_NTASKS --majel-mem=$MAJEL_MEM --sra-array=$SRA_ARRAY --sample-name=$SAMPLE_NAME \
--project-dir=$PROJECT_DIR --rsync-time=$RSYNC_TIME --rsync-mem=$RSYNC_MEM --mail-user=$EMAIL --genome=$GENOME \
--genome-path=$GENOME_PATH --majel-dir=$MAJEL_DIR --majel-args=\"${MAJEL_ARGS}\" &>> $LOG_FILE"
    fi

    #print job parameters and sbatch submission for user to check
    printf '\n%s\n\n' "These parameters will result in the following slurm submission:" | tee -a $LOG_FILE
    printf '%s\n\n' "$SUBMISSION" | tee -a $LOG_FILE

    printf '%s\n\n' "These parameters will result in the following Majel.py job:"
    printf '%s\n\n' "python3 $MAJEL_DIR/Majel.py --data_dir $PROJECT_DIR/$SAMPLE_NAME/data/ --sample_name $SAMPLE_NAME \
--genome $GENOME --genome_path $GENOME_PATH --aligner_threads $ALIGNER_THREADS ${MAJEL_ARGS}\
-v 3 -L $PROJECT_DIR/$SAMPLE_NAME/${SAMPLE_NAME}_majel.log &> slurm_majel_stdout.log"

    #skip user confirmation step if --skip-promt argument was set
    if [[ $SKIP_PROMPT == "false" ]]; then
        read -p "Would you like to continue?(Y/N)" -n 1 -r
        printf '\n'
    fi

    #print confirmation of user input
    TIME=$(date '+%B %d %T %Z %Y')
    if [[ $REPLY =~ ^[Yy]$ ]] || [[ $SKIP_PROMPT != 'false' ]]; then
        printf '%s\n\n' "$TIME> job was submitted to slurm" | tee -a $LOG_FILE
        eval $SUBMISSION
        printf '~%.0s' {1..150} | tee -a $LOG_FILE
        printf '\n\n' | tee -a $LOG_FILE
    else
        printf '%s\n\n' "$TIME> job was not submitted to slurm" | tee -a $LOG_FILE
    fi
}

#intilize argument variables
set_defaults
get_arguments "$@"

if [ -z $HELP ]; then
    printf '\n'
    printf '%s\n' 'usage: initilize_majel_submission.sh [--help] [--project-dir=<path>] [--sample-name=<name>] [--majel-time=<time>] [--majel-ntasks=<ntasks>] [--majel-mem=<size[units]> [--rsync-time=<time>] [--rsync-mem=<size[units]>'
    printf '\n'
    printf '%s\n' 'Mandatory arguments:'
    printf '%s\n' '  --project-dir=         sets path to the project directory containing the sample directory (e.g. --project-dir=/scratch1/usr001/PRJNA123456)'
    printf '%s\n' '  --sample-name=         sets name of the sample to run through Majel.py pipeline (e.g. --sample-name=Tissue_Subtissue_CancerType_SampleInfo_SAMN12345678)'
    printf '%s\n' '                         the sample name must correspond to a directory in the project directory and conform to the Majel.py naming conventions'
    printf '%s\n' '                         the sample directory must contain a data/ directory containing either SRA (.sra) or FASTQ (.fq or .fastq) files'
    printf '%s\n' '                         i.e. /path/to/PROJECT_NAME/SAMPLE_NAME/data/file.sra'
    printf '%s\n' '  --mail-user=           sets email for SLURM notifications'
    printf '\n'
    printf '%s\n' 'Optional arguments:'
    printf '%s\n' '  --sample-file=         sets path to a tab or comma delimited file containing sample information to run through pipeline (e.g. --sample-name=/scratch1/usr001/samples.txt)'
    printf '%s\n' '  --concurrent-jobs      if -sample-file= is delcared, this is the maximum number of Majel.py jobs submitted to slurm from the <sample-file> at any one time'
    printf '%s\n' '                         Note: for each line of this sample file any arguments generated or contained will override those declared globally'
    printf '%s\n' '  --sra-array=           sets the list of SRAs to be downloaded (e.g. --SRA-array="SRR1234567 SRR1234568" OR --SRA-array=SRR123456{7..8} )'
    printf '%s\n' '                         Note: as per the example a list of SRAs must be contained within quotation marks'
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
    printf '%s\n' '  --majel-dir=           used to alter path to Majel.py'
    printf '%s\n' '                         default: /datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main'
    printf '%s\n' '  --majel-args=          used to add additional arguments to Majel.py (e.g. --majel-args="--pbat --is_paired_end False"'
    printf '%s\n' '  --skip-prompt          skips verification step - use when user input is not possible (i.e. when submitting script with sbatch)'
    printf '\n'
    exit 1
fi

#check if SAMPLE_FILE was declared and if so run through this file in a job submission loop, 
#pauses when maximum number of concurrent jobs is reached and waits for a job to finished before submitting another
if [[ ! -z $SAMPLE_FILE ]]; then
    if [[ -f $SAMPLE_FILE ]]; then
        SKIP_PROMPT='true'

        CSV_FILE=$(cat $SAMPLE_FILE | tr "\\t" ",")
        JOB_ARRAY=()
        while read line
        do
            set_defaults

            CSV=(); while read -rd,; do CSV+=("$REPLY"); done <<<"$line,"
            SAMPLE_NAME=$(IFS=_ ; echo "${CSV[*]:0:4}" | sed 's/ [[:lower:]]/\U&/g' | tr -d ' ')
            if [[ $SAMPLE_NAME != "Tissue_Sub-Tissue_DiseaseStatus_SampleInfo" ]]; then
                ARGS=(); while read -rd\;; do ARGS+=("$REPLY"); done <<<"${CSV[6]}; "
                get_arguments "$@"
                get_arguments "${ARGS[@]}"

                if [[ ! -z ${CSV[4]} ]]; then
                    SRA_ARRAY="\"$(echo ${CSV[4]} | tr ';' ' ')\""
                fi

                job_submission
                JOB_ARRAY+=("$PROJECT_DIR/$SAMPLE_NAME")
            fi 
          
            i=0

            if [[ ${#JOB_ARRAY[@]} -ge $CONCURRENT_JOBS ]]; then
                printf '%s\n\n' "The maximum number of concurrent Majel.py jobs (--concurrent-jobs=$CONCURRENT_JOBS) has been reached. Waiting for a job to finish before continung with slurm submissions." | tee -a $LOG_FILE
            fi

            while [[ ${#JOB_ARRAY[@]} -ge $CONCURRENT_JOBS ]]
            do
                if grep -q 'methylseekrAndTDF' ${JOB_ARRAY[i]}/slurm_majel_stdout.log \
                || grep -q 'methylseekrAndTDF' ${JOB_ARRAY[i]}/slurm_submission_stdout.log \
                || grep -q 'Error:' ${JOB_ARRAY[i]}/slurm_submission_stdout.log; then
                    unset JOB_ARRAY[i]
                    JOB_ARRAY=("${JOB_ARRAY[@]}")
                fi
            
                ((i++))
            
                if [[ $i -ge ${#JOB_ARRAY[@]} ]]; then
                    i=0
                fi
                sleep 10m
            done

        done <<<"$CSV_FILE"

        printf '%s\n\n' "All jobs from $SAMPLE_FILE have been submitted. Exiting submission pipeline." | tee -a $LOG_FILE

    else
        echo "Error: the file \"$SAMPLE_FILE\" does not exist. Exiting pipeline"
    fi
else
    #submit a single job if SAMPLE_FILE was not declared
    job_submission
fi
