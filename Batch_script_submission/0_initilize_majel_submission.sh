#!/bin/bash

HELP='false'
SKIP_PROMPT='false'

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
    SYNC_TO='/datasets/work/hb-meth-atlas/work/Data/level_2/public'

    GENOME='hg38'
    GENOME_PATH='/datasets/work/hb-meth-atlas/work/pipeline_data/Genomes/'
    ALIGNER_THREADS='6'

    #fetches the directory from which this script is located and run
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
    #the Majel.py script is located one directory up from to bash scripts directory
    MAJEL_DIR="$(dirname "$SCRIPT_DIR")"

    unset FILE_LIST
    unset DATA_DIR
    unset MAJEL_ARGS
}


#function checks if arguments have been set and prints their values
check_argument() {
   if [[ -z "$2" ]]; then
        printf '%s\n\n' "Error: --${1}= argument not set" | tee -a $LOG_FILE
	FAIL='true'
   fi
   printf '%s\n' "--${1}=$2" | tee -a $LOG_FILE $PARAMETERS_FILE
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
        --file-list=*)
          if [[ -z $FILE_LIST ]]; then
              FILE_LIST="${1#*=}"
          else
              FILE_LIST+=" ${1#*=}"
          fi
          ;;
	 --data-dir=*)
           DATA_DIR="${1#*=}"
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

set_fixes() {   
    if hash slurm 2> /dev/null; then
        SUB_PREFIX="$1"
        unset SUB_SUFFIX
    else
        SUB_PREFIX='nohup '
	SUB_SUFFIX=' >/dev/null 2>&1 </dev/null &'
    fi
}


job_submission() {
    
    SAMPLE_DIR="$PROJECT_DIR/$SAMPLE_NAME"
    mkdir -p "$SAMPLE_DIR/data"

    LOG_FILE="$SAMPLE_DIR/0_initilize_majel_submission.log"
    > $LOG_FILE

    PARAMETERS_FILE="$SAMPLE_DIR/run_parameters_$SAMPLE_NAME.txt"
    > $PARAMETERS_FILE

    #check the values of mandatory arguments 
    check_argument "project-dir" $PROJECT_DIR
    check_argument "sample-name" $SAMPLE_NAME

    #print info on the SAMPLE_FILE used to run multiple Majel jobs at once
    if [[ ! -z $SAMPLE_FILE ]]; then
        printf '%s\n' "--sample-file=$SAMPLE_FILE"
        printf '%s\n' "--concurrent-jobs=$CONCURRENT_JOBS"
    fi
    
    #check if this pipeline is being run on cluster with a slurm submission system, if so specifying an email is mandatory
    if hash slurm 2> /dev/null; then
        check_argument "mail-user" $EMAIL
    fi
        
    #covert path variable to full path (i.e. ./destionation -> /full/path/to/destionation)
    PROJECT_DIR="$(cd $PROJECT_DIR && pwd)"

    #if SRAs were not specified for download make sure sequence files exist in /data directory
    if [[ ! -z $DATA_DIR ]] && [[ ! -z $FILE_LIST ]]; then
        for FILE_PREFIX in $FILE_LIST
        do
	    FILES=$(find $DATA_DIR/ -type f -iname "$FILE_PREFIX*.fq.gz" -o -iname "$FILE_PREFIX*.fastq.gz" -o -iname "$FILE_PREFIX*.sra")

	    for FILE in $FILES
            do
                echo $(basename $FILE)
		ln -s $FILE $PROJECT_DIR/$SAMPLE_NAME/data/$(basename $FILE)
            done

            check_argument "file-list" "$FILE_LIST"
            check_argument "data-dir" "$DATA_DIR"
        done 
    elif [[ -z $FILE_LIST ]]; then
        COUNT=`find $PROJECT_DIR/$SAMPLE_NAME/data -type f \( -iname \*.fastq.gz -o -iname \*.sra -o -iname \*.fq.gz \) | wc -l`
        if [ $COUNT = 0 ]; then
            printf '%s\n' "Error: IF --file-list= is not used to specify SRA file(s) for download" | tee -a $LOG_FILE	
	    printf '%s\n' "Error: OR --file-list= and --data-dir= are not used to specify files for soft linking" | tee -a $LOG_FILE
            printf '%s\n\n' "Error: THEN $PROJECT_DIR/$SAMPLE_NAME/data must contain SRA (.sra) or FASTQ (.fastq.gz OR fq.gz) files" | tee -a $LOG_FILE
	    FAIL='true'
        fi
    else
        if [[ $FILE_LIST != *"}"* ]]; then
            #removes potential " arifacts that can arise when --file-list= is specified in the SAMPLE_FILE as a string series
            FILE_LIST="$(echo $FILE_LIST | tr -d '"')"
        else
            #converts a specified range ( e.g. --file-list=SRR1{1..3} ) into a string series ( e.g. --file-list="SRR1 SRR2 SRR3" )
            FILE_LIST=$FILE_LIST
        fi

        check_argument "file-list" "$FILE_LIST"
    fi

    cd $PROJECT_DIR/$SAMPLE_NAME

    #print the values of all arguments for user to examine
    check_argument "majel-time" $MAJEL_TIME
    check_argument "majel-ntaskts" $MAJEL_NTASKS
    check_argument "majel-mem" $MAJEL_MEM
    check_argument "rsync-time" $RSYNC_TIME
    check_argument "rsync-mem" $RSYNC_MEM
    check_argument "sync-to" $SYNC_TO
    check_argument "aligner-threads" $ALIGNER_THREADS
    check_argument "genome" $GENOME
    check_argument "genome-path" $GENOME_PATH
    check_argument "majel-dir" $MAJEL_DIR
    
    #add space after additional majel arguments, if they exist
    [[ -z $MAJEL_ARGS ]] || MAJEL_ARGS=$(echo "$MAJEL_ARGS " | tr -d '"')
    
    #prints MAJEL_ARGS values without checking for existance (this argument is not required for subsequent steps)
    printf '%s\n\n' "--majel-args=\"${MAJEL_ARGS}\"" | tee -a $LOG_FILE $PARAMETERS_FILE

    #print time and script inputs
    TIME=$(date '+%B %d %T %Z %Y')
    printf '%s'     "$TIME> $BASH_SOURCE --project-dir=$PROJECT_DIR --sample-name=$SAMPLE_NAME --mail-user=$EMAIL --majel-time=$MAJEL_TIMES " | tee -a $LOG_FILE
    printf '%s'     " --majel-ntasks=$MAJEL_NTASKS $MAJEL_MEM --rsync-tim=$RSYNC_TIME --rsync-mem=$RSYNC_MEM --aligner-threads=$ALIGNER_THREADS --genome=$GENOME " | tee -a $LOG_FILE
    printf '%s\n'   " --genome-path=$GENOME_PATH --majel-dir=$MAJEL_DIR --majel-args=\"${MAJEL_ARGS}\" --file-list=\"$FILE_LIST\"" | tee -a $LOG_FILE

    if [[ -z $FILE_LIST ]]; then

        set_fixes "sbatch --time=$MAJEL_TIME --ntasks-per-node=$MAJEL_NTASKS --mem=$MAJEL_MEM --job-name=MAJEL:${SAMPLE_NAME} --mail-user=$EMAIL "
 
        #set job submission variable for if --file-list= was not declared (skips SRA download script)
        SUBMISSION="${SUB_PREFIX}$SCRIPT_DIR/2_sbatch_majel_submission.sh --sample-name=$SAMPLE_NAME --project-dir=$PROJECT_DIR "
        SUBMISSION+="--rsync-time=$RSYNC_TIME --rsync-mem=$RSYNC_MEM --mail-user=$EMAIL --genome=$GENOME --genome-path=$GENOME_PATH "
        SUBMISSION+="--sync-to=$SYNC_TO --majel-dir=$MAJEL_DIR --majel-args=\"${MAJEL_ARGS}\"$SUB_SUFFIX"
    else
        #convert string series to proper array
        TMP_ARRAY=($FILE_LIST)
        LEN_ARRAY=$(expr ${#TMP_ARRAY[@]} - 1)
        #calculate number of cores based on number of SRA files to be downloaded (5 file downloads per core, as advised by IMT)
        CORES=$(($LEN_ARRAY/5+1))
 
        set_fixes "sbatch --job-name=SRA_DL:${SAMPLE_NAME} --mail-user=$EMAIL --ntasks-per-node=$CORES "

        ##set job submission variable for if --file-list= was declared
        SUBMISSION="${SUB_PREFIX}$SCRIPT_DIR/1_sbatch_parallel_sra_wget.sh --majel-time=$MAJEL_TIME --majel-ntasks=$MAJEL_NTASKS --majel-mem=$MAJEL_MEM "
        SUBMISSION+="--file-list=$FILE_LIST --sample-name=$SAMPLE_NAME --project-dir=$PROJECT_DIR --rsync-time=$RSYNC_TIME --rsync-mem=$RSYNC_MEM "
        SUBMISSION+="--sync-to=$SYNC_TO --mail-user=$EMAIL --genome=$GENOME --genome-path=$GENOME_PATH --majel-dir=$MAJEL_DIR --majel-args=\"${MAJEL_ARGS}\"$SUB_SUFFIX"
    fi

    #print job parameters and sbatch submission for user to check
    printf '%s\n\n' "These parameters will result in the following submission:" | tee -a $LOG_FILE
    printf '%s\n\n' "$SUBMISSION" | tee -a $LOG_FILE

    printf '%s\n\n' "These parameters will result in the following Majel.py job:"
    printf '%s'     "python3 $MAJEL_DIR/Majel.py --data_dir $PROJECT_DIR/$SAMPLE_NAME/data/ --sample_name $SAMPLE_NAME " | tee -a $PARAMETERS_FILE
    printf '%s'     "--genome $GENOME --genome_path $GENOME_PATH --aligner_threads $ALIGNER_THREADS ${MAJEL_ARGS}" | tee -a $PARAMETERS_FILE
    printf '%s\n\n' "-v 3 -L $PROJECT_DIR/$SAMPLE_NAME/${SAMPLE_NAME}_majel.log &> slurm_majel_stdout.log" | tee -a $PARAMETERS_FILE

    #skip user confirmation step if --skip-promt argument was set
    if [[ $SKIP_PROMPT == "false" ]]; then
        read -p "Would you like to continue?(Y/N)" -n 1 -r
        printf '\n'
    fi

    #confirmation of user input
    TIME=$(date '+%B %d %T %Z %Y')
    if [[ $REPLY =~ ^[Yy]$ ]] || [[ $SKIP_PROMPT != 'false' ]] && [[ -z $FAIL ]]; then
        printf '%s\n\n' "$TIME> job was submitted" | tee -a $LOG_FILE
        
        #submit job
        eval "$SUBMISSION" | tee -a $LOG_FILE
    else
	if [[ -z $FAIL ]]; then
            printf '%s\n\n' "$TIME> job was not submitted" | tee -a $LOG_FILE
	else
            printf '%s\n\n' "$TIME> job submission failed due to an error" | tee -a $LOG_FILE
            unset FAIL
	fi
    fi
        
    printf '\n' | tee -a $LOG_FILE
    printf '~%.0s' {1..150}
    printf '\n\n'

}


#intilize argument variables
set_defaults
get_arguments "$@"

if [ -z $HELP ]; then
    printf '\n'
    printf '%s'   'usage: 0_initilize_majel_submission.sh [--help] [--project-dir=<path>] [--sample-name=<name>] [--majel-time=<time>] [--majel-ntasks=<ntasks>] '
    printf '%s\n' '[--majel-mem=<size[units]> [--rsync-time=<time>] [--rsync-mem=<size[units]>]'
    printf '\n'
    printf '%s\n' 'Mandatory arguments:'
    printf '%s\n' '  --project-dir=         sets path to the project directory containing the sample directory (e.g. --project-dir=/scratch1/usr001/PRJNA123456)'
    printf '%s\n' '  --sample-name=         sets name of the sample to run through Majel.py pipeline (e.g. --sample-name=Tissue_Subtissue_CancerType_SampleInfo_SAMN12345678)'
    printf '%s\n' '                         the sample name must correspond to a directory in the project directory and conform to the Majel.py naming conventions'
    printf '%s\n' '                         the sample directory must contain a data/ directory containing either SRA (.sra) or FASTQ (.fq or .fastq) files'
    printf '%s\n' '                         i.e. /path/to/PROJECT_NAME/SAMPLE_NAME/data/file.sra'
    printf '%s\n' '  --mail-user=           sets email for SLURM notifications'
    printf '%s\n' '                         not required when running on a system that does not use SLURM (e.g. the SHIRO-RI workstation)'
    printf '\n'
    printf '%s\n' 'Optional arguments:'
    printf '%s\n' '  --sample-file=         sets path to a tab or comma delimited file containing sample information to run through pipeline (e.g. --sample-name=/scratch1/usr001/samples.csv)'
    printf '%s\n' '  --concurrent-jobs      if -sample-file= is delcared, this is the maximum number of Majel.py jobs submitted from the <sample-file> at any one time'
    printf '%s\n' '                         Note: for each line of this sample file any arguments generated or contained will override those declared globally'
    printf '%s\n' '  --file-list=           sets the list of sequence files for downloaded or softlinking (e.g. --file-list="SRR1234567 SRR1234568" OR --SRA-array=SRR123456{7..8} )'
    printf '%s\n' '                         Note: as per the example the list must be contained within quotation marks and sperated by spaces'
    printf '%s\n' '  --data-dir             sets the directory from where the files specified in --file-list= will be softlinked'
    printf '%s\n' '  --majel-time=          sets --time= allocated to sbatch_majel_submission_AJ.sh     (default: --majel-time=04-00)'
    printf '%s\n' '  --majel-ntasks=        sets --ntasks-per-node= for sbatch_majel_submission_AJ.sh   (default: --majel-ntasks=20)'
    printf '%s\n' '  --majel-mem=           sets --mem= allocated to sbatch_majel_submission_AJ.sh      (default: --majel-mem=60gb'
    printf '%s\n' '  --rsync-time=          sets time allocated to sbatch_io_SyncProcessedData_AJ.sh    (default: --rsync-time=08:00:00)'
    printf '%s\n' '  --rsync-mem=           sets memory allocated to sbatch_io_SyncProcessedData_AJ.sh  (default: --rsync-mem=4gb'
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
    printf '%s\n' '  --skip-prompt          skips user confirmation (y/n) step'
    printf '\n'
    exit 1
fi

LEDGER_TITLES="Data Type,Sample_Name,Organism,Tissue,Sub-Tissue,Disease Status,User Name,Date Processed,"
LEDGER_TITLES+="Citation,Repository,Input File(s),Notes,hg38,Directory location 1,Directory location 2,Notes,"

ledger_check() {
    
    LEDGER_FILE="$SCRIPT_DIR/$1_sample_ledger.csv"
    if [[ ! -f $LEDGER_FILE ]]; then
        echo -e "$LEDGER_TITLES" > "$LEDGER_FILE"
    fi

    echo -e "$2" >> "$LEDGER_FILE"

    SAMPLE_RECORD="$PROJECT_DIR/$SAMPLE_NAME/$1_record_$SAMPLE_NAME.csv"
    echo -e "$2" > "$SAMPLE_RECORD"
}

CHECK_PHRASES="methylseekrAndTDF Error:"
CHECK_FILES="slurm_majel_stdout.log 0_initilize_majel_submission.log 1_sbatch_parallel_sra_wget.log 2_sbatch_majel_submission.log"

completion_check() {

    CHECK='false'	
    for FILE in $CHECK_FILES; do
        for PHRASE in $CHECK_PHRASES; do
	    if grep -q -s "$PHRASE" $FILE; then
                CHECK='true'
		break 1
		break 2
            fi
        done
    done
}


submission_cycle() {
   
    i=0
    while [[ ${#SAMPLE_ARRAY[@]} -ge $1 ]]
    do
        completion_check
       
        if [[ $CHECK='true' ]]; then
            
            SAMPLE_NAME=$(basename ${SAMPLE_ARRAY[i]})
	    PROJECT_DIR=$(dirname ${SAMPLE_ARRAY[i]})
            PROJECT_NAME=$(basename ${PROJECT_DIR})
            
	    CSV=(); while read -rd,; do CSV+=("$(echo "${REPLY^}")"); done <<<"${CSV_ARRAY[i]},"
            DATE=$(date '+%Y-%m-%d')
	 
	    if [[ ! -z ${CSV[5]} ]]; then
                USER_NAME=${CSV[5]}
            elif [[ ! -z $EMAIL ]]; then
                USER_NAME="$(echo  "${EMAIL^}" | cut -d @ -f 1 | tr '.' ' ' | sed 's/ [[:lower:]]/\U&/g')"
            else
                USER_NAME=$USER
            fi
            LEDGER_INFO=",${SAMPLE_NAME},,${CSV[0]},${CSV[1]},${CSV[2]},$USER_NAME,${DATE},,${PROJECT_NAME},"
            LEDGER_INFO+="${CSV[4]},,${GENOME},${SYNC_TO}/${PROJECT_NAME}/${SAMPLE_NAME},,,"    
            
	    if grep -q -s "MethylSeekR and toTDF Completed" ${SAMPLE_ARRAY[i]}/slurm_majel_stdout.log; then
                ledger_check "completed" "$LEDGER_INFO"
            else
                ledger_check "failed" "$LEDGER_INFO"
            fi
    
            unset SAMPLE_ARRAY[i]
            SAMPLE_ARRAY=("${SAMPLE_ARRAY[@]}")
            unset CSV_ARRAY[i]
            CSV_ARRAY=("${CSV_ARRAY[@]}")
    
        fi
    
        ((i++))
    
        if [[ $i -ge ${#SAMPLE_ARRAY[@]} ]]; then
            i=0
        fi

        sleep 10m
    done
}


#check if SAMPLE_FILE was declared and if so run through this file in a job submission loop, 
#pauses when maximum number of concurrent jobs is reached and waits for a job to finished before submitting another
if [[ ! -z $SAMPLE_FILE ]]; then
    if [[ -f $SAMPLE_FILE ]]; then
        SKIP_PROMPT='true'

        CSV_FILE=$(cat $SAMPLE_FILE | tr "\\t" "," | sed 's/ [[:lower:]]/\U&/g')
        SAMPLE_ARRAY=()
	PROJECT_ARRAY=()
	CSV_ARRAY=()

        while read line
        do
            set_defaults

	    CSV=(); while read -rd,; do CSV+=("$(echo "${REPLY^}")"); done <<<"$line,"
            SAMPLE_NAME=$(IFS=_ ; echo "${CSV[*]:0:4}" | tr -d ' ')
            if [[ $SAMPLE_NAME != "Tissue_Sub-Tissue_DiseaseStatus_SampleInfo" ]]; then
                ARGS=(); while read -rd\;; do ARGS+=("$REPLY"); done <<<"${CSV[6]}; "
                get_arguments "$@"
                get_arguments "${ARGS[@]}"

                if [[ ! -z ${CSV[4]} ]]; then
                    FILE_LIST="\"$(echo ${CSV[4]} | tr ';' ' ')\""
                fi

                job_submission
                SAMPLE_ARRAY+=("$PROJECT_DIR/$SAMPLE_NAME")
		CSV_ARRAY+=("$line")
            fi 
          
            if [[ ${#SAMPLE_ARRAY[@]} -ge $CONCURRENT_JOBS ]]; then
                printf '%s'     "The maximum number of concurrent Majel.py jobs (--concurrent-jobs=$CONCURRENT_JOBS) has been reached. "
		printf '%s\n\n' "Waiting for a job to finish before continung with submissions."
            fi
            
            submission_cycle $CONCURRENT_JOBS

        done <<<"$CSV_FILE"

        submission_cycle 1
       
        printf '%s\n\n' "All jobs from $SAMPLE_FILE have completed. Exiting submission pipeline."

    else
        printf '%s\n\n' "Error: the file \"$SAMPLE_FILE\" does not exist. Exiting submission pipeline"
    fi
else
    #submit a single job if SAMPLE_FILE was not declared
    job_submission 
    SAMPLE_ARRAY+=("$PROJECT_DIR/$SAMPLE_NAME")
    SAMPLE_INFO=(""$(echo $SAMPLE_NAME | tr '_' ' ')"")
    line="$(echo "${SAMPLE_INFO[*]:0:3}" | tr ' ' ','),,,$(echo $FILE_LIST | tr ' ' ';' | tr -d '"'),,,"
    CSV_ARRAY+=($line)
    submission_cycle 1
fi
