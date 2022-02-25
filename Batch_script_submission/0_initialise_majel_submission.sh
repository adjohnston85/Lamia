#!/bin/bash
#SBATCH --time=05-00
#SBATCH --output=sbatch.out

HELP='false'
SKIP_PROMPT='false'
JOB_DIR=$PWD
TRANSFER_CHECK='false'

#function to set or reset default values
set_defaults() {
    
    CONCURRENT_JOBS=5

    RSYNC_TIME='08:00:00'
    RSYNC_MEM='512mb'
    SYNC_TO='/datasets/work/hb-meth-atlas/work/Data/level_2/public'
    SQLITE_DIR='/datasets/work/hb-meth-atlas/work/Data/level_2'

    GENOME='hg38'
    GENOME_PATH='/datasets/work/hb-meth-atlas/work/pipeline_data/Genomes'

    if hash slurm 2> /dev/null; then
        SCRIPT_DIR="/datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main/Batch_script_submission"
	NO_GENOME_TRANSFER='false'
    else
        #fetches the directory from which this script is located and run
        SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
	NO_GENOME_TRANSFER='true'
    fi

    #the Majel.py script is located one directory up from to bash scripts directory
    MAJEL_DIR="$(dirname "$SCRIPT_DIR")"

    DL_ATTEMPTS='5'

    RUN_FILES=()

    WHOLE_EXPERIMENT='false'
    EXPERIMENT_ACCESSION='false'
    SKIP_DL='false'
   
    unset TRIM_PROFILE
    unset MAJEL_MEM
    unset MAJEL_NTASKS
    unset MAJEL_TIME
    unset NO_RSYNC
    unset DL_ONLY
    unset RSYNC_ONLY
    unset SAMPLE_NAME
    unset PROJECT_NAME
    unset RUN_LIST
    unset RUN_DIR
    unset MAJEL_ARGS
    unset FAIL
    unset SQL_CHECK
}


#function checks if arguments have been set and prints their values
check_argument() {

    printf '%s\n' "$1$2" | tee -a $LOG_FILE $PARAMETERS_FILE
    if [[ -z "$2" ]] || [[ $2 == '\"\"' ]]; then
        printf '\n%s\n\n' "Error: ${1} argument not set" | tee -a $LOG_FILE
        FAIL='true'
    fi
}


#function to grab input argument values
get_arguments() {

    while [[ $# -gt 0 ]]; do
    case "$1" in
        --job-dir=*)
          JOB_DIR="${1#*=}"
          ;;
	--sqlite-dir=*)
          SQLITE_DIR="${1#*=}"
	  ;;
	--project-name=*)
          PROJECT_NAME="${1#*=}"
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
        --run-list=*)
          if [[ -z $RUN_LIST ]]; then
              RUN_LIST="${1#*=}"
          else
              RUN_LIST+=" ${1#*=}"
          fi
          ;;
        --whole-experiment)
          WHOLE_EXPERIMENT='true'
	  ;;
        --skip-dl)
          SKIP_DL='true'
          ;;
        --dl-only)
          DL_ONLY="--dl-only "
          ;;
	--rsync-only)
          RSYNC_ONLY="--rsync-only "
          ;;
        --no-rsync)
          NO_RSYNC="--no-rsync "
          ;; 
	--experiment-accession)
          EXPERIMENT_ACCESSION='true'
	  ;;
        --dl-attempts=*)
          DL_ATTEMPTS="${1#*=}"
          ;;
	 --run-dir=*)
          RUN_DIR="${1#*=}"
          ;;
        --majel-time=*)
          MAJEL_TIME="--majel-time=${1#*=} "
	  TIME="--time=${1#*=} "
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
        --no-genome-transfer)
          NO_GENOME_TRANSFER='true'
          ;;	  
        --majel-dir=*)
          MAJEL_DIR="${1#*=}"
          SCRIPT_DIR="$MAJEL_DIR/Batch_script_submission"
          ;;
        --majel-args=*)
          MAJEL_ARGS="${1#*=}"
          ;;
        --trim-profile=*)
          TRIM_PROFILE="--trim_profile ${1#*=} "
          ;;
        --skip-prompt)
          SKIP_PROMPT='true'
          ;;
        --help)
          unset HELP
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


#Check if reference genome files for alignement exist in the job directory, if not transfer them
transfer_ref_genome() {

    REF_LOG=$JOB_DIR/refgenome_rsync.out
    THREADS=8
    if hash slurm 2> /dev/null && [[ $NO_GENOME_TRANSFER == 'false' ]]; then
        > $REF_LOG
        sbatch --ntasks-per-node=$THREADS --time=00:20:00 --output=$REF_LOG --mail-type=ALL --mail-user=$EMAIL --wrap "module load rclone/1.55.1 ; rclone copy $GENOME_PATH/$GENOME/ $JOB_DIR/$GENOME/ --progress --multi-thread-streams=$THREADS ; echo 'transfer complete' >> $REF_LOG"
        echo 'transferring genome folder to --job-dir='
        while [[ $TRANSFER_CHECK == 'false' ]]; do
            if [[ ! -z $(grep "transfer complete" $REF_LOG) ]]; then
                TRANSFER_CHECK='true'
                sleep 10s
            fi
        done
    fi
}


job_submission() {
 
    [[ $NO_GENOME_TRANSFER == 'true' ]] || GENOME_PATH="$JOB_DIR"
    SAMPLE_NAME=$(echo $SAMPLE_NAME | tr -d ' ')
    PROJECT_NAME=$(echo $PROJECT_NAME | tr -d ' ')
    PROJECT_DIR="$JOB_DIR/$PROJECT_NAME"

    if [[ ! -z $RSYNC_ONLY ]]; then

        SUBMISSION="$MAJEL_DIR/majel_cleanup.sh $SAMPLE_NAME &>> slurm_majel_stdout.log"
        eval $SUBMISSION

        SUBMISSION="${RSYNC_PREFIX}$SCRIPT_DIR/3_sbatch_io_SyncProcessedData.sh --sync-to=$SYNC_TO --sync-from=$PROJECT_DIR/$SAMPLE_NAME &>> slurm_majel_stdout.log"
    else

        mkdir -p "$PROJECT_DIR/$SAMPLE_NAME/data"
   
        #create log files to overwrite any existing log files, as these are used to track job completion 	
	> $PROJECT_DIR/$SAMPLE_NAME/1_sbatch_parallel_sra_wget.log
	> $PROJECT_DIR/$SAMPLE_NAME/2_sbatch_majel_submission.log
	> $PROJECT_DIR/$SAMPLE_NAME/3_sbatch_io_SyncProcessedData.log
	
	LOG_FILE="$PROJECT_DIR/$SAMPLE_NAME/0_initialise_majel_submission.log"
        > $LOG_FILE
    
        PARAMETERS_FILE="$PROJECT_DIR/$SAMPLE_NAME/run_parameters_$SAMPLE_NAME.txt"
        > $PARAMETERS_FILE
    
        #add space after additional majel arguments, if they exist
        [[ -z $MAJEL_ARGS ]] || MAJEL_ARGS=$(echo "$MAJEL_ARGS " | tr -d '"')
    
        #print time and script inputs
        printf '%s'     "$BASH_SOURCE --project-dir=$PROJECT_DIR --sample-name=$SAMPLE_NAME --mail-user=$EMAIL $MAJEL_TIME " | tee -a $LOG_FILE
        printf '%s'     " --majel-ntasks=$MAJEL_NTASKS --majel-mem=$MAJEL_MEM --rsync-tim=$RSYNC_TIME --rsync-mem=$RSYNC_MEM --genome=$GENOME " | tee -a $LOG_FILE
        printf '%s'     " --genome-path=$GENOME_PATH --majel-dir=$MAJEL_DIR --majel-args=\"${MAJEL_ARGS}\" --run-list=\"$RUN_LIST\"" | tee -a $LOG_FILE
        printf '%s\n\n' " --run-dir=${RUN_DIR}" | tee -a $LOG_FILE 
    
        #check the values of mandatory arguments
        check_argument "--job-dir=" $JOB_DIR
        check_argument "--project-name=" $PROJECT_NAME 
        check_argument "--project-dir=" $PROJECT_DIR
        check_argument "--sample-name=" $SAMPLE_NAME
    
        #print info on the SAMPLE_FILE used to run multiple Majel jobs at once
        if [[ ! -z $SAMPLE_FILE ]]; then
            printf '%s\n' "--sample-file=$SAMPLE_FILE"
            printf '%s\n' "--concurrent-jobs=$CONCURRENT_JOBS"
        fi
        
        #check if this pipeline is being run on cluster with a slurm submission system, if so specifying an email is mandatory
        if hash slurm 2> /dev/null; then
            check_argument "--mail-user=" $EMAIL
            RSYNC_PREFIX="sbatch --time=$RSYNC_TIME --mem=$RSYNC_MEM --mail-user=$EMAIL --job-name=RSYNC:$SAMPLE_NAME "
        fi
            
        #covert path variable to full path (i.e. ./destionation -> /full/path/to/destionation)
        PROJECT_DIR="$(cd $PROJECT_DIR && pwd)"
    
        #if SRAs were not specified for download make sure sequence files exist in /data directory
        if [[ ! -z $RUN_DIR ]]; then
	    check_argument "--run-dir=" $RUN_DIR
	    if [[ -d $RUN_DIR ]]; then
	        if [[ ! -z $RUN_LIST ]]; then
		    check_argument "--run-list=" "$RUN_LIST"
                    for FILE_PREFIX in $RUN_LIST; do
		        RUN_FILES+=($(find $RUN_DIR/ -regextype posix-extended -regex ".*/($FILE_PREFIX).*_?[rR]?[12]?\.(fq|fastq)\.?(gz)?"))
		        RUN_FILES+=($(find $RUN_DIR/ -regextype posix-extended -regex ".*/($FILE_PREFIX)"))
	            done
                else
                    RUN_FILES+=($(find $RUN_DIR/ -regextype posix-extended -regex ".*/.*_?[rR]?[12]?\.(fq|fastq)\.?(gz)?"))
                fi
    
	        if [[ ! -z $RUN_FILES ]]; then
                    for FILE in ${RUN_FILES[@]}; do
	                ln -sf $FILE $PROJECT_DIR/$SAMPLE_NAME/data/$(basename $FILE)
                    done
                else
                    printf '%s\n' "Error: the run \"${FILE_PREFIX}\" specified in --run-list= does not exist in the --run-dir=$RUN_DIR directory" | tee -a $LOG_FILE
	            FAIL='true'
	        fi
            else
                printf '%s\n' "Error: the directory in --run-dir=$RUN_DIR does not exist" | tee -a $LOG_FILE
                FAIL='true'
                
                check_argument "--run-list=" "\"$RUN_LIST\""
                check_argument "--run-dir=" "$RUN_DIR"
            fi
    
        elif [[ -z $RUN_LIST ]]; then
	        RUN_FILES=($(find $PROJECT_DIR/$SAMPLE_NAME/data -regextype posix-extended -regex ".*/.*_?[rR]?[12]?\.(fq|fastq|sra)\.?(gz)?"))
	    if [[ ${#RUN_FILES[@]} == 0 ]]; then
                printf '%s\n'   "Error: IF --run-list= is not used to specify SRA file(s) for download" | tee -a $LOG_FILE	
	        printf '%s\n'   "Error: OR --run-list= and --run-dir= are not used to specify files for soft linking" | tee -a $LOG_FILE
                printf '%s\n\n' "Error: THEN $PROJECT_DIR/$SAMPLE_NAME/data must contain SRA (.sra) or FASTQ (.fastq.gz OR fq.gz) files" | tee -a $LOG_FILE
	        FAIL='true'
            fi
        else
            check_argument "--run-list=" "\"$RUN_LIST\""
	    printf '%s\n' "--run-dir=${RUN_DIR}" | tee -a $LOG_FILE $PARAMETERS_FILE
        fi
        
        #set defaults for majel memory and cores, or base one on the other if only one is stipulated 	
        if [[ -z $MAJEL_MEM ]] && [[ -z $MAJEL_NTASKS ]]; then
            MAJEL_MEM='64gb'
            MAJEL_NTASKS=32
        elif [[ -z $MAJEL_MEM ]]; then
            MAJEL_MEM="$(( $MAJEL_NTASKS * 2 ))gb"
        elif [[ -z $MAJEL_NTASKS ]]; then
            MAJEL_NTASKS=$(( ${MAJEL_MEM//[!0-9]/} / 2 ))
        fi            

        #print the values of all arguments for user to examine
        check_argument "--majel-ntaskts=" $MAJEL_NTASKS
        check_argument "--majel-mem=" $MAJEL_MEM
        check_argument "--rsync-time=" $RSYNC_TIME
        check_argument "--rsync-mem=" $RSYNC_MEM
        check_argument "--sync-to=" $SYNC_TO
        check_argument "--genome=" $GENOME
        check_argument "--genome-path=" $GENOME_PATH
        check_argument "--majel-dir=" $MAJEL_DIR
	check_argument "--trim-profile=" "$(echo $TRIM_PROFILE | cut -d' ' -f2)"
           
        if [[ $FAIL == 'true' ]]; then
            printf '%s\n\n' "job submission failed due to an error" | tee -a $LOG_FILE
            return
        fi
    
        cd $PROJECT_DIR/$SAMPLE_NAME
    
        if [[ -z $RUN_LIST ]] || [[ ! -z $RUN_DIR ]] || [[ $SKIP_DL == "true" ]]; then

            if [[ -z $MAJEL_TIME ]]; then
                #plotting size of sequence files versus majel processing time yielded the following linear equation for 32 cores and 64gb memory: y = 0.0004x + 0.8605
                FILE_SIZE=$(find $PROJECT_DIR/$SAMPLE_NAME/data -regextype posix-extended -regex ".*/*\.(fq|fastq|sra)\.?(gz)?" -print0 | du -L --files0-from=- -cm | grep total | cut -f1)
		#time for sbatch majel submission is calculated using this linear equation, adjusting for number of cores and adding 50% contingency
                TIME="--time=$(( (4 * $FILE_SIZE + 8605) * 15 / 100000 * 32 / $MAJEL_NTASKS)):00:00 "
            fi
            
	    check_argument "" $TIME 
            set_fixes "sbatch $TIME--ntasks-per-node=$MAJEL_NTASKS --mem=$MAJEL_MEM --job-name=MAJEL:${SAMPLE_NAME} --mail-user=$EMAIL "
     
            #set job submission variable for if --run-list= was not declared (skips SRA download script)
            SUBMISSION="${SUB_PREFIX}$SCRIPT_DIR/2_sbatch_majel_submission.sh --sample-name=$SAMPLE_NAME --project-dir=$PROJECT_DIR "
            SUBMISSION+="--rsync-time=$RSYNC_TIME --rsync-mem=$RSYNC_MEM --mail-user=$EMAIL --genome=$GENOME --genome-path=$GENOME_PATH "
            SUBMISSION+="--majel-threads=$MAJEL_NTASKS --sync-to=$SYNC_TO --majel-dir=$MAJEL_DIR $NO_RSYNC--majel-args=\"${TRIM_PROFILE}${MAJEL_ARGS}\"$SUB_SUFFIX"
        else
            #convert string series to proper array
            TMP_ARRAY=($RUN_LIST)
            LEN_ARRAY=$(expr ${#TMP_ARRAY[@]} - 1)
            #calculate number of cores based on number of SRA files to be downloaded (5 file downloads per core, as advised by IMT)
            CORES=$(($LEN_ARRAY/5+1))
     
            set_fixes "sbatch --job-name=SRA_DL:${SAMPLE_NAME} --mail-user=$EMAIL --ntasks-per-node=$CORES "
    
            ##set job submission variable for if --run-list= was declared
            SUBMISSION="${SUB_PREFIX}$SCRIPT_DIR/1_sbatch_parallel_sra_wget.sh $MAJEL_TIME--majel-ntasks=$MAJEL_NTASKS --majel-mem=$MAJEL_MEM "
            SUBMISSION+="--run-list=\"$RUN_LIST\" --dl-attempts=$DL_ATTEMPTS --sample-name=$SAMPLE_NAME --project-dir=$PROJECT_DIR --rsync-time=$RSYNC_TIME --rsync-mem=$RSYNC_MEM "
            SUBMISSION+="--sync-to=$SYNC_TO --mail-user=$EMAIL --genome=$GENOME --genome-path=$GENOME_PATH --majel-dir=$MAJEL_DIR $NO_RSYNC$DL_ONLY--majel-args=\"${TRIM_PROFILE}${MAJEL_ARGS}\"$SUB_SUFFIX"
        fi
        
        #prints MAJEL_ARGS values without checking for existance (this argument is not required for subsequent steps)
	printf '%s\n\n' "--majel-args=\"${MAJEL_ARGS}\"" | tee -a $LOG_FILE $PARAMETERS_FILE

        if [[ -z $RUN_FILES ]]; then
            printf '%s\n\n' "Majel will run on following SRA files that will be downloaded to $PROJECT_DIR/$SAMPLE_NAME/data:" | tee -a $LOG_FILE
            RUN_FILES=($RUN_LIST)
            RUN_FILES=( "${RUN_FILES[@]/%/.sra}" )
        else
            printf '%s\n\n' "Majel will run on the following files located in $PROJECT_DIR/$SAMPLE_NAME/data:" | tee -a $LOG_FILE
        fi

        for FILE in ${RUN_FILES[@]}; do
            printf '%s' "$(basename $FILE) " | tee -a $LOG_FILE
        done

        printf '\n\n' | tee -a $LOG_FILE

        #print job parameters and sbatch submission for user to check
        printf '%s\n\n' "These parameters will result in the following submission:" | tee -a $LOG_FILE
        printf '%s\n\n' "$SUBMISSION" | tee -a $LOG_FILE
    
        if [[ -z $DL_ONLY ]] && [[ -z $RSYNC_ONLY ]]; then
            printf '%s\n\n' "These parameters will result in the following Majel.py job:" | tee -a $LOG_FILE
            printf '%s'     "python3 $MAJEL_DIR/Majel.py --data_dir $PROJECT_DIR/$SAMPLE_NAME/data/ --sample_name $SAMPLE_NAME " | tee -a $LOG_FILE
            printf '%s'     "--genome $GENOME --genome_path $GENOME_PATH --threads $MAJEL_NTASKS $TRIM_PROFILE$MAJEL_ARGS" | tee -a $LOG_FILE
            printf '%s\n\n' "-v 3 -L $PROJECT_DIR/$SAMPLE_NAME/${SAMPLE_NAME}_majel.log &> slurm_majel_stdout.log" | tee -a $LOG_FILE
        fi
    
        #skip user confirmation step if --skip-promt argument was set
        if [[ $SKIP_PROMPT == "false" ]]; then
            read -p "Would you like to continue?(Y/N)" -n 1 -r
            printf '\n'
        fi
   fi
  
    #confirmation of user input
    if [[ $REPLY =~ ^[Yy]$ ]] || [[ $SKIP_PROMPT != 'false' ]] && [[ -z $FAIL ]]; then
        if [[ $NO_GENOME_TRANSFER == 'false' ]] && [[ $TRANSFER_CHECK == 'false' ]]; then
            transfer_ref_genome
        fi 
        
        printf '%s\n\n' "Job was submitted on $(date '+%B %d %Y at %T %Z')" | tee -a $LOG_FILE
         
        #submit job
        eval "$SUBMISSION" | tee -a $LOG_FILE
    else
        printf '%s\n\n' "Job was not submitted" | tee -a $LOG_FILE
        exit 1
    fi
    
    unset QUEUED_ARRAY[0]
    QUEUED_ARRAY=("${QUEUED_ARRAY[@]}")
    SAMPLE_ARRAY+=("$PROJECT_DIR/$SAMPLE_NAME")
    CSV_ARRAY+=("$line")

    printf '~%.0s' {1..150}
    printf '\n\n'
}


LEDGER_TITLES="Data Type,Sample_Name,Organism,Tissue,Sub-Tissue,Disease Status,User Name,Date Processed,"
LEDGER_TITLES+="Citation,Repository,Input File(s),Notes,genome,Directory location 1,Directory location 2,Notes,"

ledger_check() {

    LEDGER_INFO=",${SAMPLE_NAME},,${CSV[0]},${CSV[1]},${CSV[2]},$USER_NAME,${DATE},,${PROJECT_NAME},"
    LEDGER_INFO+="${CSV[3]},,${GENOME},${SYNC_TO}/${PROJECT_NAME}/${SAMPLE_NAME},,,"
    
    LEDGER_FILE="$PROJECT_DIR/${PROJECT_NAME}_$1_sample_ledger.csv"
    if [[ ! -f $LEDGER_FILE ]]; then
        echo -e "$LEDGER_TITLES" > "$LEDGER_FILE"
    fi

    echo -e "$LEDGER_INFO" >> "$LEDGER_FILE"

    SAMPLE_RECORD="$PROJECT_DIR/$SAMPLE_NAME/$1_record_$SAMPLE_NAME.csv"
    echo -e "$LEDGER_INFO" > "$SAMPLE_RECORD"
}


completion_check() {
    
    if [[ ! -z $DL_ONLY ]]; then
        CHECK_PHRASES=("dl-only completed" "Error:")
        CHECK_FILES="0_initialise_majel_submission.log 1_sbatch_parallel_sra_wget.log slurm_majel_stdout.log"
    else
        CHECK_PHRASES=("MethylSeekR and toTDF Completed" "Error:" "methylseekrAndTDF did not complete" "sending incremental file list")
        CHECK_FILES="0_initialise_majel_submission.log 1_sbatch_parallel_sra_wget.log 2_sbatch_majel_submission.log 3_sbatch_io_SyncProcessedData.log"
    fi

    CHECK='false'
    for FILE in $CHECK_FILES; do
        for PHRASE in "${CHECK_PHRASES[@]}"; do
            if grep -q -s "$PHRASE" $PROJECT_DIR/$SAMPLE_NAME/$FILE; then
                CHECK='true'
                break 2
            fi
        done
    done
}


report_jobs() {

    printf '%s\n' "Queued samples:"
    for SAMPLE in "${QUEUED_ARRAY[@]}"; do
        printf '%s\n' "$SAMPLE"
    done

    printf '\n%s\n' "Samples currently running:"
    for SAMPLE in "${SAMPLE_ARRAY[@]}"; do
        printf '%s\n' "$SAMPLE"
    done

    printf '\n%s\n' "Failed samples:"
    for SAMPLE in "${FAIL_ARRAY[@]}"; do
        printf '%s\n' "$SAMPLE"
    done
    
    printf '\n%s\n' "Samples completed successfully:"
    for SAMPLE in "${SUCCESS_ARRAY[@]}"; do
        printf '%s\n' "$SAMPLE"
    done
    printf '\n'
}


submission_cycle() {

    i=0
    while [[ ${#SAMPLE_ARRAY[@]} -ge $1 ]]; do
        SAMPLE_NAME=$(basename ${SAMPLE_ARRAY[i]})
        PROJECT_DIR=$(dirname ${SAMPLE_ARRAY[i]})
        PROJECT_NAME=$(basename $PROJECT_DIR)
        completion_check

        if [[ $CHECK == 'true' ]]; then

            CSV=(); while read -rd,; do CSV+=("$(echo "${REPLY^}")"); done <<<"${CSV_ARRAY[i]},"
            DATE=$(date '+%Y-%m-%d')

            if [[ ! -z ${CSV[4]} ]]; then
                USER_NAME=${CSV[4]}
            elif [[ ! -z $EMAIL ]]; then
                USER_NAME="$(echo  "${EMAIL^}" | cut -d @ -f 1 | tr '.' ' ' | sed 's/ [[:lower:]]/\U&/g')"
            else
                USER_NAME=$USER
            fi

            if grep -q -s "MethylSeekR and toTDF Completed" ${SAMPLE_ARRAY[i]}/slurm_majel_stdout.log || \
               grep -q -s "sending incremental file list" ${SAMPLE_ARRAY[i]}/3_sbatch_io_SyncProcessedData.log || \
	       grep -q -s "dl-only completed" ${SAMPLE_ARRAY[i]}/1_sbatch_parallel_sra_wget.log; then
                ledger_check "completed"
                SUCCESS_ARRAY+=("${SAMPLE_ARRAY[i]}")
            elif grep -q -s "dl-only completed" ${SAMPLE_ARRAY[i]}/1_sbatch_parallel_sra_wget.log; then
                SUCCESS_ARRAY+=("${SAMPLE_ARRAY[i]}")
            else
                ledger_check "failed"
                FAIL_ARRAY+=("${SAMPLE_ARRAY[i]}")
            fi

            unset SAMPLE_ARRAY[i]
            SAMPLE_ARRAY=("${SAMPLE_ARRAY[@]}")
            unset CSV_ARRAY[i]
            CSV_ARRAY=("${CSV_ARRAY[@]}")

            if [[ ! -z $2 ]]; then
                report_jobs
            fi
        fi

        ((i++))

        if [[ $i -ge ${#SAMPLE_ARRAY[@]} ]]; then
            i=0
        fi

        sleep 10s
    done
}


sql_dl() {
    echo $1

    if [[ $1 == " " ]]; then
        printf '%s\n' "Error: $2 could not be located in $SQLITE_DIR/SRAmetadb.sqlite."
        printf '%s\n' "The latest version of SRAmetadb.sqlite can be downloaded from https://s3.amazonaws.com/starbuck1/sradb/SRAmetadb.sqlite.gz"

        SQ_URL='https://s3.amazonaws.com/starbuck1/sradb/SRAmetadb.sqlite.gz'
        SQ_BACKUP='https://gbnci-abcc.ncifcrf.gov/backup/SRAmetadb.sqlite.gz'

	SQL_CHECK='fail'
    fi
}


get_run_list() {

    FIRST_RUN=($(echo $1 | tr ';' ' '))
    #if --experiment-accession was specified, use this to procure then run files
    if [[ $EXPERIMENT_ACCESSION == 'true' ]]; then
        RUN_LIST=$(sqlite3 $SQLITE_DIR/SRAmetadb.sqlite "select run_accession from run where experiment_accession='$FIRST_RUN'")
	printf '%s\n' "--experiment-accession option selected. Locating SRA files from $FIRST_RUN specified in --run-list in $SQLITE_DIR/SRAmetadb.sqlite..."
        sql_dl "$RUN_LIST " $EXPERIMENT_ACCESSION
	[[ $SQL_CHECK != 'fail' ]] || printf '%s\n\n' "If $EXPERIMENT_ACCESSION cannot be located in the latest SRAmetadb.sqlite version, this job will need to be run by specifying all SRA files for download in --run-list="
    #if --whole-experiment option was set then use the SRA file specified in the table to procure all run files from the same experiment
    elif [[ $WHOLE_EXPERIMENT == 'true' ]]; then
        EXPERIMENT_ACCESSION=$(sqlite3 $SQLITE_DIR/SRAmetadb.sqlite "select experiment_accession from run where run_accession='$FIRST_RUN'")
	printf '%s\n' "--whole-experiment option selected. $FIRST_RUN (the first run accession specified in --run-list=) will be used to locate other SRA files from the same experiment" 
	printf '%s\n' "Locating other SRA files in $SQLITE_DIR/SRAmetadb.sqlite with the same experiment accession $EXPERIMENT_ACCESSION..."
	[[ -z $EXPERIMENT_ACCESSION ]] || RUN_LIST=$(sqlite3 $SQLITE_DIR/SRAmetadb.sqlite "select run_accession from run where experiment_accession='$EXPERIMENT_ACCESSION'")
	sql_dl "$EXPERIMENT_ACCESSION " "$FIRST_RUN"
	[[ $SQL_CHECK != 'fail' ]] || printf '%s\n\n' "If $FIRST_RUN cannot be located in the latest SRAmetadb.sqlite version then --whole-experiment cannot be used and this job will need to be run by specifying all SRA files for download in --run-list="
    else
        RUN_LIST=$(echo "$1" | tr ';' ' ')
    fi

    #At this point the run list may be in a string or range (e.g. SRR123{1..4}) format. This coverts range to string.
    RUN_LIST=$(eval 'echo ${RUN_LIST[@]} | tr -d "\""')
    RUN_ARRAY=($RUN_LIST)
}


#intilize argument variables
set_defaults
get_arguments "$@"

if [ -z $HELP ]; then
    set_defaults
    printf '\n%s\n'   'usage: 0_initialise_majel_submission.sh [--help] [--trim-profile=<profile>] [--job-dir=<path>] [--sample-name=<name>] [--mail-user=<email>]'
    printf '%s\n'     '                                        [--project-name=<name>] [--run-list=<list of run accessions (SRAs) OR fastq file prefixes>]'
    printf '%s\n'     '                                        [--run-dir=<path>] [--whole-experiment] [--experiment-accession] [--sample-file=<file>]'
    printf '%s\n'     '                                        [--concurrent-jobs] [--dl-attempts=<integer>]  [--dl-only] [--skip-dl] [--majel-time=<time>]'
    printf '%s\n'     '                                        [--majel-ntasks=<integer>] [--majel-mem=<integer+mb/gb>] [--rsync-time=<time>]'
    printf '%s\n'     '                                        [--rsync-mem=integer+Mb/Gb>] [--sync-to=<path>] [--genome=<genome>] [--genome-path=<path>]'
    printf '%s\n'     '                                        [--majel-dir=<path>] [--majel-args=<Majel.py arguments>] [--sqlite-dir=<path>]'
    printf '%s\n\n'   '                                        [--skip-prompt] [--no-genome-transfer] [--rsync-only] [--no-rsync]'
    printf '%s\n'     'arguments:'
    printf '%s\n'     '  --help                 show this help message and exit'
    printf '%s\n'     "  --trim-profile=        Sets the profile for number of base pairs trimmed from 5'and 3' ends of sequence reads (post adapter trimming)"
    printf '%s\n'     "                         Options are: swift, em-seq, & no-trim, a single integer to cut from all ends, or a comma seperated list of 4 integers (R1 5', R2 5', R1 3', R2 3')"
    printf '%s\n'     "  --job-dir=             sets path to the job directory where the project and sample directories will reside, e.g. --job-dir=/scratch1/usr001 (default: your current working directory)"
    printf '%s\n'     '  --sample-name=         sets name of the sample to run through Majel.py pipeline (e.g. --sample-name=Tissue_Subtissue_CancerType_SampleInfo_SAMN12345678)'
    printf '%s\n'     '  --mail-user=           sets email for SLURM notifications'
    printf '%s\n'     '  --project-name=        sets name of the project. This will be used to create a project directory (e.g. --project-name=PRJN12334)'
    printf '%s\n'     '  --run-list=            sets the list of run accessions (SRAs) or prefixes of FASTQ files for downloaded or soft linking (e.g. --run-list="SRR1234567 SRR1234568" OR --run-list=SRR123456{7..8} )'
    printf '%s\n'     '                         Note: as per the example the list must be contained within quotation marks and sperated by spaces'
    printf '%s\n'     '  --run-dir=             sets the directory where the file name prefixes specified in --run-list= will be used to locate and soft link FASTQ files'   
    printf '%s\n'     '  --whole-experiment     specify a single run accession (SRA) in run-list=. All other runs from the same experiment will be located, downloaded, and processed'  
    printf '%s\n'     '  --experiment-accession specify an experiment accession (rather than a run accession) in --run-list=. This accession is used to locate, download, and process all runs (SRAs) included in this experiment'
    printf '%s\n\n'   '  --sample-file=         sets path to a tab or comma delimited file containing sample information to run through pipeline (e.g. --sample-name=/scratch1/usr001/samples.csv)'
    printf '%s\n'     '                         <sample-file> layout:'
    printf '%s\n\n'   '                         Tissue          Subtissue      SampleType     <run-list>      Username       Argument(1)       Argument(2)        Argument(3)  ...   Argument(n)'
    printf '%s\n'     '                         SRA <sample-file> example:'
    printf '%s\n\n'   '                         colon,colon,Normal_adjacent_tissue_P6,SRR949213 SRR949214 SRR949215,Andrew Johnston'
    printf '%s\n'     '                         FASTQ <sample-file> example:'
    printf '%s\n\n'   '                         breast,breast,normal tissue,breast_run1;breast_run2,Andrew Johnston,--project-name=BRE001,--run-dir=/path/to/files,--sample-name=Breast_NOS_NormalTissue_BRE001'
    printf '%s\n'     '                         Note: arguments are generally not required in this file for run accessions (SRAs), as <sample-name> and <project-name> will be determined from SRAmetadb.sqlite, if possible'
    printf '%s\n'     '                               arguments that apply to all samples included in the <sample-file> such as <job-dir> and <mail-user> can be declared globally'    
    printf '%s\n'     '                               any arguments declared in this file will override those declared globally'
    printf '%s\n'     '                               accessions or file prefixes specified in Column 4 must be sperated by spaces or semicolons'
    printf '%s\n\n'   '                               columns 6 and above are used to place arguments that apply specifically to the sample job within the row'
    printf '%s\n'     "  --concurrent-jobs=     if -sample-file= is delcared, this is the maximum number of Majel.py jobs submitted from the <sample-file> at any one time (default: --concurrent-jobs=$CONCURRENT_JOBS)"
    printf '%s\n'     "  --dl-attempts=         sets the number of failed attempts to download an SRA file before the pipeline exits on an error (default: -dl-attempts=$DL_ATTEMPTS)"
    printf '%s\n'     '                         e.g. if --dl-attempts=1 the pipeline will not reattempt failed SRA downloads'
    printf '%s\n'     '  --dl-only              downloads SRA files but skips subsequent steps including running through Majel.py'
    printf '%s\n'     '  --rsync-only           jumps to final rsync step enacted by 3_sbatch_io_SyncProcessedData.sh'
    printf '%s\n'     '  --no-rsync             files will not be rsynced after 2_sbatch_majel_submission.sh completes'
    printf '%s\n'     '  --skip-dl              skips the SRA download step and goes directly to Majel submission. For use when SRA files have already been downloaded'
    printf '%s\n'     '  --majel-time=          sets --time= allocated to 2_sbatch_majel_submission.sh (default: --majel-time= is set based on size of seqence files in data/ directory)'
    printf '%s\n'     "  --majel-ntasks=        sets --ntasks-per-node= for 2_sbatch_majel_submission.sh and number of cores used by Majel.py (default: --majel-ntasks=$MAJEL_NTASKS)"
    printf '%s\n'     "  --majel-mem=           sets --mem= allocated to 2_sbatch_majel_submission.sh (default: --majel-mem=$MAJEL_MEM"
    printf '%s\n'     '                         Note: if only one of--majel-mem= or --majel-ntasks= is set, the other will automatically adjust by a factor of 2 (i.e. memory [gb] is set to twice the ntasks [cores]'
    printf '%s\n'     "  --rsync-time=          sets time allocated to 3_sbatch_io_SyncProcessedData.sh (default: --rsync-time=$RSYNC_TIME)"
    printf '%s\n'     "  --rsync-mem=           sets memory allocated to 3_sbatch_io_SyncProcessedData.sh (default: --rsync-mem=$RSYNC_MEM"
    printf '%s\n'     "  --sync-to=             sets path to directory being synced to (Default: --sync-to=$SYNC_TO)"
    printf '%s\n'     "  --genome=              used to alter --genome argument for Majel.py (default: --majel-genome=$GENOME)"
    printf '%s\n'     "  --genome-path=         used to alter --genome_path argument for Majel.py (default: --majel-genome-path=$GENOME_PATH)"
    printf '%s\n'     "  --majel-dir=           used to alter path to Majel.py (default: $MAJEL_DIR)"
    printf '%s\n'     '  --majel-args=          used to add additional arguments to Majel.py (e.g. --majel-args="--pbat --is_paired_end False")'
    printf '%s\n'     '                         Note: this list of additional arguments must be contained within quotation marks'
    printf '%s\n'     '  --sqlite-dir=          directory containing SRAmetadb.sqlite'
    printf '%s\n'     '  --skip-prompt          skips user confirmation (y/n) step'
    printf '%s\n'     '  --no-genome-transfer   disables files from --genome-path= being transferred to the --job-dir=, which occurs by default to greatly speed up SLURM job initiation' 
    printf '\n'
    exit 1
fi

#check if SAMPLE_FILE was declared and if so run through this file in a job submission loop, 
#pauses when maximum number of concurrent jobs is reached and waits for a job to finished before submitting another
if [[ ! -z $SAMPLE_FILE ]]; then
    [[ -f $SAMPLE_FILE ]] || (printf '%s\n\n' "Error: the file \"$SAMPLE_FILE\" does not exist. Exiting submission pipeline" ; exit 1)
         
    SKIP_PROMPT='true'
    CSV_FILE=$(cat $SAMPLE_FILE | tr "\\t" ",")
    SAMPLE_ARRAY=()
    PROJECT_ARRAY=()
    CSV_ARRAY=()
    FAIL_ARRAY=()
    SUCCESS_ARRAY=()
    readarray -t QUEUED_ARRAY < $SAMPLE_FILE
    
    while read line; do
        CSV=(); while read -rd,; do CSV+=("$REPLY"); done <<<"$line,"
        set_defaults
        get_arguments "$@"
        get_arguments "${CSV[@]:5}"

	SAMPLE_NO=$(cat $SAMPLE_FILE | sed '/^\s*$/d' | wc -l)
	[[ $SAMPLE_NO -ge $CONCURRENT_JOBS ]] || CONCURRENT_JOBS=$SAMPLE_NO
 
        #if the --run-list was not specified then grab the run files from the --sample-file table 
        if [[ -z $RUN_LIST ]]; then
            RUN_LIST=${CSV[3]}
        fi
 
        get_run_list "$RUN_LIST"
 
        #If --sample-name was not specified then make the sample name from sample-file table. Use SRA to determine experiment accession and project
        if [[ -z $SAMPLE_NAME ]]; then
            SAMPLE_NAME=$(IFS=_ ; CELLS="${CSV[*]:0:3}"; echo "${CELLS^}" | sed 's/ [[:lower:]]/\U&/g' | tr -d ' ' | sed 's/_[[:lower:]]/\U&/g')
 
            if [[ -z $RUN_DIR ]]; then
                EXPERIMENT_ACCESSION=$(sqlite3 $SQLITE_DIR/SRAmetadb.sqlite "select experiment_accession from experiment where experiment_accession in (select experiment_accession from run where run_accession='${RUN_ARRAY[0]}')")
		if [[ $EXPERIMENT_ACCESSION ]]; then
                    SAMPLE_NAME+="_$EXPERIMENT_ACCESSION"
                else
		    SAMPLE_NAME+="_${RUN_ARRAY[0]}"
                fi
            else
                SAMPLE_NAME+="_${RUN_ARRAY[0]}"
            fi
 
            if [[ -z $PROJECT_NAME ]]; then
                SQL_QUERY="select study_accession from study where study_accession in (select study_accession from experiment where "
                SQL_QUERY+="experiment_accession in (select experiment_accession from run where run_accession='${RUN_ARRAY[0]}'))"
                printf '%s\n' "--project-name= option was not specfied. ${RUN_ARRAY[0]} (the first run accession specified in --run-list=) will be used to locate the study accession for use as the project name"
                printf '%s\n' "Locating study accession in $SQLITE_DIR/SRAmetadb.sqlite..."
                PROJECT_NAME=$(sqlite3 $SQLITE_DIR/SRAmetadb.sqlite "$SQL_QUERY")
	        
                sql_dl "$PROJECT_NAME " ${RUN_ARRAY[0]}
                [[ ! -z $PROJECT_NAME ]] || printf '%s\n' "If ${RUN_ARRAY[0]} cannot be located in the latest SRAmetadb.sqlite version, this job will need to be run by specifying --project-name="
 
   	        if [[ $SQL_CHECK == 'fail' ]]; then
                    FAIL_ARRAY+=("${SAMPLE_ARRAY[i]}")
                    printf '~%.0s' {1..150}
                    printf '\n\n'
                    continue
                fi
            fi
        fi

        job_submission
 
        if [[ ${#SAMPLE_ARRAY[@]} -ge $CONCURRENT_JOBS ]]; then
            printf '%s'     "The maximum number of concurrent Majel.py jobs (--concurrent-jobs=$CONCURRENT_JOBS) has been reached. "
            printf '%s\n\n' "Waiting for a sample to finish before continuing with submissions."
 
            report_jobs
        fi
 
        submission_cycle $CONCURRENT_JOBS
 
    done <<<"$CSV_FILE"
 
    submission_cycle 1 0
 
    printf '%s\n\n' "All jobs from $SAMPLE_FILE have completed. Exiting submission pipeline."

else
    #submit a single job if SAMPLE_FILE was not declared
    get_run_list "$RUN_LIST"
    job_submission 
    SAMPLE_INFO=(""$(echo $SAMPLE_NAME | tr '_' ' ')"")
    line="$(echo "${SAMPLE_INFO[*]:0:3}" | tr ' ' ','),,,$(echo $RUN_LIST | tr ' ' ';' | tr -d '"'),,,"
    CSV_ARRAY+=($line)
fi
