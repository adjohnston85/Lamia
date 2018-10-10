Majel WGBS Pipeline

This pipeline (nicknamed Majel) will process paired and single end SRA files or paired end FASTQ files.

Output
The pipeline will output;
	Processed bam file & index (aligned, sorted and duplicate marked)
	Some summary statistics (average genomic coverage, read counts, duplication rate)
	CpG methylation calls in bedGraph format (with --merge-context)
	Methylation bias plots
	FASTQC reports

A large number of intermediate files are also produced. These take considerably more space than the final output, particularly towards the end of mapping 
or when using SRA as input. Due to the way ruffus works and the occasional need for error recovery, most intermediate files cannot be deleted on the fly. 
Instead, these files can be removed and the final files put into a tidy directory structure using the cleanup script. When using SRA as input, the user is 
responsible for cleaning up the intermediate fastq files from fastq-dump. 

Wall time can be long (>5 days) for high coverage files, especially if the input is SRA. A standard 20-30x WGBS run from fastq will take 2-3 days.

The pipleine has been developed using the follwoing software;
	FastQC v0.11.5
	bismark-0.18.1 see https://www.bioinformatics.babraham.ac.uk/projects/bismark/
	Methyldackel-0.3.0 (using HTSlib version 1.2.1)	see https://github.com/dpryan79/MethylDackel
	trim_galore-0.4.3 see https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
	fastq-dump : 2.8.2 (from sra toolkit) see https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
	samtools-1.4.1 (using htslib 1.4.1) see http://www.htslib.org/doc/samtools.html
	picard MarkDuplicates version 2.9.4-1-gcda9516-SNAPSHOT see https://broadinstitute.github.io/picard/command-line-overview.html
