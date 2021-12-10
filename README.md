# Majel Whole Genome Bisulfite Sequencing (WGBS) processing pipeline

This data processing pipeline (nicknamed Majel) will process both paired- and single-end SRA files (from NCBI SRA database)\
or paired-end FASTQ files. Fastq files can be obtained from some public data repositories, directly from your sequencing \
provider or using bcl2fastq (on illumina output). SRA files can be obtained from [NCBI SRA database](https://www.ncbi.nlm.nih.gov/) using your prefered \
method (e.g. wget).\
\
Majel will automate all subsequent standard processing steps, including conversion from sra to fastq (fastq-dump), read \
trmming & FASTQC (trim-galore), mapping (bismark), duplicate marking (picard), methylation bias assesment & CpG methylation \
calling (methyldackel), borwser tracks, UMR/LMR/PMD calls, as well as produce some summary statistics.

![Majel Process](./images/MajelFlowchart.png)

## Output 
The pipeline will output the following.

* Processed bam file & index (aligned, sorted and duplicate marked)
* Some summary statistics (average genomic coverage, read counts, duplication rate)
* CpG methylation calls in bedGraph format (with --merge-context)
* QC files
   * Methylation bias plots
   * FASTQC reports
* MethylSeekR UMR, LMR and PMD calls plus QC plots
* TDF file for viewing in IGV

\
A large number of intermediate files are also produced, these can be removed (and the final directory structure formed) using the \
included cleanup script\
```> cleanup.sh sampleID```

These intermediate files can be very large when working with high-coverage datasets. This can be especially true during the \
bismark mapping steps or when running from SRA files. When running from SRA, cleanup.sh will not remove fastq files produced\
by fastq-dump. The user must remove these manually. This is to prevent the accidental deletion of raw fastq files when running on\
user generated data. This can also be prevented by good data practices (e.g. soft linking raw data).

## Using Majel
Majel makes use of the Python pipelining module ruffus (see [ruffus docs](http://www.ruffus.org.uk/)). Expect a long walltime on\
high coverage datasets (>5 days).\
\
Genome data files (e.g. indexes and CpG islands) will need to be created before running. These should be stored at a logical\
path in a folder named after the genome build. CpG island files will be USCS track files (bed like) names [GENOME]_CpGislands.txt.\

Majel help is available using the '--help' flag.

Basic usage:\
cd to data folder\
python3 Path_to_majel_wgbspipline/Majey.py --data_dir Path_to_data/ --genome hg19/hg_38/mm10 --file_type sra/FASTQ --sampleID Name_of_tissue --genomePath Path_to_genome_folder/ -v 3 -L Path_to_data/Log_file \

```
usage: Majel.py [-h] [--verbose [VERBOSE]] [--version] [-L FILE] [-T JOBNAME]
                [-j N] [--use_threads] [-n] [--touch_files_only]
                [--recreate_database] [--checksum_file_name FILE]
                [--flowchart FILE] [--key_legend_in_graph]
                [--draw_graph_horizontally] [--flowchart_format FORMAT]
                [--forced_tasks JOBNAME] [--genome GENOME]
                [--data_dir DATA_DIR] [--sample_name SAMPLE_NAME]
                [--aligner_threads ALIGNER_THREADS] [--pbat]
                [--is_paired_end IS_PAIRED_END] [--genome_path GENOME_PATH]

Majel.py - Automated WGBS processing pipeline for sra and fastq file types

optional arguments:
  -h, --help            show this help message and exit
  --genome GENOME       Genome for alignment. Must be a directory in your
                        --genome_path. Defaults to hg38 (hg38)
  --data_dir DATA_DIR   Directory for input fastq/sra files
  --sample_name SAMPLE_NAME
                        Sample name. Used for output file names and should
                        match directory name. Will determine colouring of TDF
                        track file and must conform to following convention
                        "Tissue_SubTissue_HealthStatus_Identifier"
  --aligner_threads ALIGNER_THREADS
                        Speed up alignment by increasing number of threads.
                        Values depends on the aligner. Defaults to 5
  --pbat                Specify when aligning pbat library
  --is_paired_end IS_PAIRED_END
                        Is the libarary paired end (defaults to True)
  --genome_path GENOME_PATH
                        Path to genome folder, must contain a --genome
                        directory

Common options:
  --verbose [VERBOSE], -v [VERBOSE]
                        Print more verbose messages for each additional
                        verbose level.
  --version             show program's version number and exit
  -L FILE, --log_file FILE
                        Name and path of log file

pipeline arguments:
  -T JOBNAME, --target_tasks JOBNAME
                        Target task(s) of pipeline.
  -j N, --jobs N        Allow N jobs (commands) to run simultaneously.
  --use_threads         Use multiple threads rather than processes. Needs
                        --jobs N with N > 1
  -n, --just_print      Don't actually run any commands; just print the
                        pipeline.
  --touch_files_only    Don't actually run the pipeline; just 'touch' the
                        output for each task to make them appear up to date.
  --recreate_database   Don't actually run the pipeline; just recreate the
                        checksum database.
  --checksum_file_name FILE
                        Path of the checksum file.
  --flowchart FILE      Don't run any commands; just print pipeline as a
                        flowchart.
  --key_legend_in_graph
                        Print out legend and key for dependency graph.
  --draw_graph_horizontally
                        Draw horizontal dependency graph.
  --flowchart_format FORMAT
                        format of dependency graph file. Can be 'pdf', 'svg',
                        'svgz' (Structured Vector Graphics), 'pdf', 'png'
                        'jpg' (bitmap graphics) etc
  --forced_tasks JOBNAME
                        Task(s) which will be included even if they are up to
                        date.
```
 
## Required software
Majel was written using the following packages
* FastQC v0.11.5
* bismark-0.18.1 https://www.bioinformatics.babraham.ac.uk/projects/bismark/
* Methyldackel-0.3.0 (using HTSlib version 1.2.1) https://github.com/dpryan79/MethylDackel
* trim_galore-0.4.3 https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
* fastq-dump : 2.8.2 (from sra toolkit) https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
* samtools-1.4.1 (using htslib 1.4.1) http://www.htslib.org/doc/samtools.html
* picard MarkDuplicates version 2.9.4-1-gcda9516-SNAPSHOT https://broadinstitute.github.io/picard/command-line-overview.html
* igvtools (from IGV Version 2.3.95)
* bedtools 2.26.0

\
Majel will call the following python modules
* ruffus
* os
* subprocess
* time
* shlex
* re
* pandas
* sys

\
Majel requires the following R packages
* littler
* optparse
* MethylSeekR
* data.table
* GenomicRanges
* rtracklayer
* BSgenome
* BSgenome.Hsapiens.UCSC.hg19
* BSgenome.Hsapiens.UCSC.hg38
* BSgenome.Mmusculus.UCSC.mm10
* parallel

## Majel submission scripts

The Batch_script_submission directory contains four BASH wrapper scripts used to simplify the process of submitting sequence files to Majel.py and SLURM.
These BASH submission scripts allow users to input settings without altering the BASH submission script code.
These wrapper scripts also produce detailed logs of parameter inputs and pipeline outputs and errors. 
If SLURM is present these scripts are submitted as jobs using the 'sbatch' command, otherwise they are run as background jobs using 'nohup' and '&'.

The scripts are as follows:

### 0_initilize_majel_submission.sh
This script is the master controller and is used to take all user input and submit jobs to run through the entire processing pipeline. 
This is where users will set all of their Majel run settings and is used to run through the entire pipeline, including:
- sequnce file download or soft linking
- running sequnce files through Majel.py
- cleanup of temporary files
- rsycning

When running this script on a single sample this script can and should be executed without submitting to SLURM.
When run locally without submitting with the 'sbatch' command, the input parameters and details of the ensuing job submissions are displayed to a user.
At this point the user is prompted to confirm whether or not to continue with the submission.
This gives an opportunity to cancel the submission process before it is engaged. 

0_initilize_majel_submission.sh also allows users to submit file with a list of samples to submit to SLURM and/or run through the Majel pipeline.
Jobs are submitted based on the maximum number of concurrent jobs set by the user and the user prompt/confirmation function is disabled.
These jobs are then monitored and subsequent samples in the list are submitted when a job fails or completes.
If running on a system using SLURM, this script should be submitted with 'sbatch'.

Basic usage for previously generated FASTQ files:\
Path_to_majel_wgbspipline/Batch_script_submission/0_initilize_majel_submission.sh --job-dir=/path/to/pipeline_output
 --run-dir=/path/to/fastq_files --run-list="fastq_prefix_1 fastq_prefix_2" --sample-name=tissue_subtissue_condition_sampleID --project-name=project_ID --mail-user=email_address
 
 '''
usage: 0_initilize_majel_submission.sh [--help] [--job-dir=<path>] [--sample-name=<name>] [--mail-user=<email>]
                                       [--project-name=<name>] [--run-list=<list of run accessions (SRAs) OR fastq file prefixes>]
                                       [--run-dir=<path>] [--whole-experiment] [--experiment-accession] [--sample-file=<file>]
                                       [--concurrent-jobs] [--dl-attempts=<integer>]  [--dl-only] [--skip-dl] [--majel-time=<time>]
                                       [--majel-ntasks=<integer>] [--majel-mem=<integer+Mb/Gb>] [--rsync-time=<time>]
                                       [--rsync-mem=integer+Mb/Gb>] [--sync-to=<path>] [--genome=<genome>] [--genome-path=<path>]
                                       [--majel-dir=<path>] [--majel-args=<Majel.py arguments>] [--sqlite-dir=<path>]
                                       [--skip-prompt] [--no-genome-transfer]

arguments:
  --job-dir=             sets path to the job directory containing the project and sample directories (e.g. --job-dir=/scratch1/usr001)
  --sample-name=         sets name of the sample to run through Majel.py pipeline (e.g. --sample-name=Tissue_Subtissue_CancerType_SampleInfo_SAMN12345678)
  --mail-user=           sets email for SLURM notifications
  --project-name=        sets name of the project. This will be used to create a project directory (e.g. --project-name=PRJN12334)
  --run-list=            sets the list of run accessions (SRAs) or prefixes of FASTQ files for downloaded or soft linking (e.g. --run-list="SRR1234567 SRR1234568" OR --run-list=SRR123456{7..8} )
                         Note: as per the example the list must be contained within quotation marks and sperated by spaces
  --run-dir=             sets the directory where the file name prefixes specified in --run-list= will be used to locate and soft link FASTQ files
  --whole-experiment     specify a single run accession (SRA) in run-list=. All other runs from the same experiment will be located, downloaded, and processed
  --experiment-accession specify an experiment accession (rather than a run accession) in --run-list=. This accession is used to locate, download, and process all runs (SRAs) included in this experiment
  --sample-file=         sets path to a tab or comma delimited file containing sample information to run through pipeline (e.g. --sample-name=/scratch1/usr001/samples.csv)

                         <sample-file> layout:
                         Tissue          Subtissue      SampleType     <run-list>      Username       Argument(1)       Argument(2)        Argument(3)  ...   Argument(n)

                         SRA <sample-file> example:
                         colon,colon,Normal_adjacent_tissue_P6,SRR949213 SRR949214 SRR949215,Andrew Johnston

                         FASTQ <sample-file> example:
                         breast,breast,normal tissue,breast_run1;breast_run2,Andrew Johnston,--project-name=BRE001,--run-dir=/path/to/files,--sample-name=Breast_NOS_NormalTissue_BRE001

                         Note: arguments are generally not required in this file for run accessions (SRAs), as <sample-name> and <project-name> will be determined from SRAmetadb.sqlite, if possible.
                               arguments that apply to all samples included in the <sample-file> such as <job-dir> and <mail-user> can be declared globally
                               any arguments declared in this file will override those declared globally
                               accessions or file prefixes specified in Column 4 must be sperated by spaces or semicolons
                               columns 6 and above are used to place arguments that apply specifically to the sample job within the row

  --concurrent-jobs      if -sample-file= is delcared, this is the maximum number of Majel.py jobs submitted from the <sample-file> at any one time
  --dl-attempts=         sets the number of failed attempts to download an SRA file before the pipeline exits on an error (default: -dl-attempts=5)
                         e.g. if --dl-attempts=1 the pipeline will not reattempt failed SRA downloads
  --dl-only              downloads SRA files but skips subsequent steps including running through Majel.py
  --skip-dl              skips the SRA download step and goes directly to Majel submission. For use when SRA files have already been downloaded
  --majel-time=          sets --time= allocated to 2_sbatch_majel_submission_AJ.sh (default: --majel-time=04-00)
  --majel-ntasks=        sets --ntasks-per-node= for 2_sbatch_majel_submission_AJ.sh and number of cores used by Majel.py (default: --majel-ntasks=32)
  --majel-mem=           sets --mem= allocated to 2_sbatch_majel_submission_AJ.sh (default: --majel-mem=64gb
  --rsync-time=          sets time allocated to 3_sbatch_io_SyncProcessedData_AJ.sh (default: --rsync-time=08:00:00)
  --rsync-mem=           sets memory allocated to 3_sbatch_io_SyncProcessedData_AJ.sh (default: --rsync-mem=512mb
  --sync-to=             sets path to directory being synced to (Default: --sync-to=/datasets/work/hb-meth-atlas/work/Data/level_2/public)
  --genome=              used to alter --genome argument for Majel.py (default: --majel-genome=hg38)
  --genome-path=         used to alter --genome_path argument for Majel.py (default: --majel-genome-path=/datasets/work/hb-meth-atlas/work/pipeline_data/Genomes/)
  --majel-dir=           used to alter path to Majel.py (default: /datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main)
  --majel-args=          used to add additional arguments to Majel.py (e.g. --majel-args="--pbat --is_paired_end False")
                         Note: this list of additional arguments must be contained within quotation marks
  --sqlite-dir=          directory containing SRAmetadb.sqlite
  --skip-prompt          skips user confirmation (y/n) step
  --no-genome-transfer   disables files from --genome-path= being transferred to the --job-dir=, which occurs by default to greatly speed up SLURM job initiation

### 1_sbatch_parallel_sra_wget.sh
This script is called by 0_initilize_majel_submission.sh if the user inputs run accessions (SRAs) for download.
SRAs specified by the user are downloaded using aria2c (much faster than wget) and will automatically reattempt failed downloads without losing progress.
For redundancy, wget is also attempted if aria2c fails.

### 2_sbatch_majel_submission.sh
This script is called either by 0_initilize_majel_submission.sh if no SRAs are specified for download or by 1_sbatch_parallel_sra_wget.sh upon completion of SRA downloads.
2_sbatch_majel_submission.sh is primarily a wrapper script for Majel.py

### 3_sbatch_io_SyncProcessedData.sh
This script will be called by 2_sbatch_majel_submission.sh upon the successful completion of a Majel.py job and is used to sync data to the DNA methylation atlas or other location for long-term storage.

Each of these BASH scripts can be used independently of previous steps.
For example, users can run the 3_sbatch_io_SyncProcessedData.sh script in a folder to sync files without having first run scripts 0-2.
