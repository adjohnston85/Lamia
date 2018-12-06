#Majel Whole Genome Bisulfite Sequencing (WGBS) processing pipeline\

This data processing pipeline (nicknamed Majel) will process both paired- and single-end SRA files (from NCBI SRA database)\
or paired-end FASTQ files. Fastq files can be obtained from some public data repositories, directly from your sequencing \
provider or using bcl2fastq (on illumina output). SRA files can be obtained from [NCBI SRA database](https://www.ncbi.nlm.nih.gov/) using your prefered \
method (e.g. wget).\
\
Majel will automate all subsequent standard processing steps, including conversion from sra to fastq (fastq-dump), read \
trmming & FASTQC (trim-galore), mapping (bismark), duplicate marking (picard), methylation bias assesment & CpG methylation \
calling (methyldackel), as well as produce some summary statistics.\
\
![Majel Process](/images/MajelFlowchart.png)
\
##Output\
The pipeline will output the following.\

* Processed bam file & index (aligned, sorted and duplicate marked)
* Some summary statistics (average genomic coverage, read counts, duplication rate)
* CpG methylation calls in bedGraph format (with --merge-context)
* QC files
   * Methylation bias plots
   * FASTQC reports
\
A large number of intermediate files are also produced, these can be removed (and the final directory structure formed) using the \
included cleanup script\
	cleanup.sh sampleID\
These intermediate files can be very large when working with high-coverage datasets. This can be especially true during the \
bismark mapping steps or when running from SRA files. When running from SRA, cleanup.sh will not remove fastq files produced\
by fastq-dump. The user must remove these manually. This is to prevent the accidental deletion of raw fastq files when running on\
user generated data. This can also be prevented by good data practices (e.g. soft linking raw data).\
\
##Using Majel\
Majel makes use of the Python pipelining module ruffus (see [ruffus docs](http://www.ruffus.org.uk/)). Expect a long walltime on\
high coverage datasets (>5 days).\
\
Majel help is available using the '--help' flag.\
\
usage: Majel.py [-h] [--verbose [VERBOSE]] [--version] [-L FILE] [-T JOBNAME]\
                [-j N] [--use_threads] [-n] [--touch_files_only]\
                [--recreate_database] [--checksum_file_name FILE]\
                [--flowchart FILE] [--key_legend_in_graph]\
                [--draw_graph_horizontally] [--flowchart_format FORMAT]\
                [--forced_tasks JOBNAME] [--genome GENOME]\
                [--data_dir DATA_DIR] [--sampleID SAMPLEID]\
                [--aligner_threads ALIGNER_THREADS] [--pbat PBAT]\
                [--file_type FILE_TYPE] [--isPairedEnd ISPAIREDEND]\
                [--bowenPath BOWENPATH]\
\
Majel.py - Automated WGBS processing pipeline\
\
optional arguments:\
  -h, --help            show this help message and exit\
  --genome GENOME       Genome reads are aligned too. Check genome files for\
                        your prefered aligner are in /media/bowen_work/pipelin\
                        e_data/majel_wgbspipline/main/Genomes. Defaults to\
                        hg38 (hg38)\
  --data_dir DATA_DIR   Directory for fastq files\
  --sampleID SAMPLEID   Sample name. Used for directory name and output file\
                        names\
  --aligner_threads ALIGNER_THREADS\
                        Speed up alignment by increasing number of threads.\
                        Values depends on the aligner. Defaults to 5\
  --pbat PBAT           True if library method uses post-bis adapter tagging\
  --file_type FILE_TYPE\
                        Starting file type (sra or fastq)\
  --isPairedEnd ISPAIREDEND\ 
                        Is the libarary paired end (defaults to True)\
  --bowenPath BOWENPATH\
                        path to bowen volume (defaults to HPC path\
                        '/OSM/CBR/HB_FSP_TBI/work/'\
\
Common options:\
  --verbose [VERBOSE], -v [VERBOSE]\
                        Print more verbose messages for each additional\
                        verbose level.\
  --version             show program's version number and exit\
  -L FILE, --log_file FILE\
                        Name and path of log file\
\
pipeline arguments:\
  -T JOBNAME, --target_tasks JOBNAME\
                        Target task(s) of pipeline.\
  -j N, --jobs N        Allow N jobs (commands) to run simultaneously.\
  --use_threads         Use multiple threads rather than processes. Needs\
                        --jobs N with N > 1\
  -n, --just_print      Don't actually run any commands; just print the\
                        pipeline.\
  --touch_files_only    Don't actually run the pipeline; just 'touch' the\
                        output for each task to make them appear up to date.\
  --recreate_database   Don't actually run the pipeline; just recreate the\
                        checksum database.\
  --checksum_file_name FILE\
                        Path of the checksum file.\
  --flowchart FILE      Don't run any commands; just print pipeline as a\
                        flowchart.\
  --key_legend_in_graph\
                        Print out legend and key for dependency graph.\
  --draw_graph_horizontally\
                        Draw horizontal dependency graph.\
  --flowchart_format FORMAT\
                        format of dependency graph file. Can be 'pdf', 'svg',\
                        'svgz' (Structured Vector Graphics), 'pdf', 'png'\
                        'jpg' (bitmap graphics) etc\
  --forced_tasks JOBNAME\
                        Task(s) which will be included even if they are up to\
                        date.\
 
##Required software\
Majel was written using the following packages\
* FastQC v0.11.5
* bismark-0.18.1 https://www.bioinformatics.babraham.ac.uk/projects/bismark/
* Methyldackel-0.3.0 (using HTSlib version 1.2.1) https://github.com/dpryan79/MethylDackel
* trim_galore-0.4.3 https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
* fastq-dump : 2.8.2 (from sra toolkit) https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
* samtools-1.4.1 (using htslib 1.4.1) http://www.htslib.org/doc/samtools.html
* picard MarkDuplicates version 2.9.4-1-gcda9516-SNAPSHOT https://broadinstitute.github.io/picard/command-line-overview.html
\
Majel will call the following python modules\
* ruffus
* os
* subprocess
* time
* shlex
* re
* pandas
