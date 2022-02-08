#!/bin/bash
#SBATCH --partition=io
#SBATCH --time=08:00:00
#SBATCH --output=MappedData_rsync.out
#SBATCH --mail-type=ALL

help='false'

SYNC_FROM=$PWD
SYNC_TO='/datasets/work/hb-meth-atlas/work/Data/level_2/public'

while [ $# -gt 0 ]; do
  case "$1" in
    --sync-from=*)
      SYNC_FROM="${1#*=}"
      ;;
    --sync-to=*)
      SYNC_TO="${1#*=}"
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
    printf '%s\n' 'Optional arguments:'
    printf '%s\n' '  --sync-from=  sets path to directory being synced from (e.g. --sync-from=/scratch1/usr001/PROJECT_NAME/SAMPLE_NAME)'
    printf '%s\n' '                Note: Do not add a / to the end of this path'
    printf '%s\n' '                Default set to current working directory'
    printf '%s\n' '  --sync-to=    sets path to directory being synced to (Default: --sync-to=/datasets/work/hb-meth-atlas/work/Data/level_2/public)'
    printf '%s\n' '                Note: the ultimate and penultimate directories specified in the --sync-from= argument will be added to this path (e.g. PROJECT_NAME/SAMPLE_NAME)'
    printf '%s\n' '                if these directories do not exist they will be created'
    
    exit 1
fi

SYNC_FROM="$(cd $SYNC_FROM && pwd)"
SAMPLE_NAME=$(basename "$SYNC_FROM")
PROJECT_FOLDER=$(basename "$(dirname $SYNC_FROM)")

SUBMISSION="rsync -avzh --whole-file --remove-source-files --no-g --chmod=Dg=rwxs --exclude '3_sbatch_io_SyncProcessedData.log' $SYNC_FROM $SYNC_TO/$PROJECT_FOLDER/ &>> 3_sbatch_io_SyncProcessedData.log"
TIME=$(date '+%B %d %T %Z %Y')
printf '%s\n\n' "$TIME> $SUBMISSION" > 3_sbatch_io_SyncProcessedData.log
eval $SUBMISSION
cp $SYNC_FROM/3_sbatch_io_SyncProcessedData.log $SYNC_TO/$PROJECT_FOLDER/$SAMPLE_NAME/3_sbatch_io_SyncProcessedData.log
