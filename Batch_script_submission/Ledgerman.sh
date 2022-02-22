LEDGER_TITLES="Data Type,Sample_Name,Organism,Tissue,Sub-Tissue,Disease Status,User Name,Date Processed,"
LEDGER_TITLES+="Citation,Repository,Input File(s),Notes,genome,Directory location 1,Directory location 2,Notes,"
SQLITE_DIR='/datasets/work/hb-meth-atlas/work/Data/level_2'

SAMPLE_NAME="$(basename $1)"
TISSUE="$(echo $SAMPLE_NAME | cut -d'_' -f1)"
SUBTISSUE="$(echo $SAMPLE_NAME | cut -d'_' -f2)"
DISEASE="$(echo $SAMPLE_NAME | cut -d'_' -f3)"
EXPERIMENT_ACCESSION="${SAMPLE_NAME##*_}"
RUN_LIST=$(sqlite3 $SQLITE_DIR/SRAmetadb.sqlite "select run_accession from run where experiment_accession='$EXPERIMENT_ACCESSION'")

if [[ ! $RUN_LIST ]]; then
    RUN_LIST=$(sqlite3 $SQLITE_DIR/SRAmetadb.sqlite "select run_accession from run where experiment_accession in (select experiment_accession from experiment where sample_accession in (select sample_accession from sample where sample_alias='SEXPERIMENT_ACCESSION'))")
fi

if [[ ! $RUN_LIST ]]; then
    RUN_LIST=$EXPERIMENT_ACCESSION
fi

RUN_ARRAY=($RUN_LIST)

printf -v joined '%s; ' "${RUN_ARRAY[@]}"

SQL_QUERY="select study_accession from study where study_accession in (select study_accession from experiment where "
SQL_QUERY+="experiment_accession in (select experiment_accession from run where run_accession='${RUN_ARRAY[0]}'))"
PROJECT_NAME=$(sqlite3 $SQLITE_DIR/SRAmetadb.sqlite "$SQL_QUERY")
DATE="$(date +%F -r $(find $1 -name *_sd.bam))"

LEDGER_INFO=",${SAMPLE_NAME},,$TISSUE,$SUBTISSUE,$DISEASE,,${DATE},,https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=${PROJECT_NAME},"
LEDGER_INFO+="$joined,,,,$1,,,"

LEDGER_FILE="$(dirname $1)/${PROJECT_NAME}_completed_sample_ledger.csv"

printf "%s\n" "$LEDGER_INFO" >> "$LEDGER_FILE"

