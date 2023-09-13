# Majel Snakemake Pipeline README

The Majel Snakemake pipeline automates the analysis of DNA methylation sequencing data. This README provides an in-depth overview of the pipeline, its rules, and instructions for usage.  
<br>

## Table of Contents

- [Introduction](#introduction)
- [Pipeline Components](#pipeline-components)
- [Usage](#usage)
- [Configuration Options](#configuration-options)
- [Pipeline Rules](#pipeline-rules)
- [Customization](#customization)
- [Output](#output)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)  
<br>

## Introduction

The Majel pipeline streamlines the analysis of DNA methylation sequencing data by automating SRA downloading, fastq trimming, sequence alignment, deduplication, methylation calling, variant detection and various QC analyses. The pipeline utilizes the Snakemake workflow management system to ensure efficient and reproducible execution, allow customizable HPC resource usage and foster modularity.  
<br>

## Pipeline Components

The pipeline comprises the following main components:

- **Rules:** Each analysis step is represented by a Snakemake rule, defining inputs, outputs, parameters, and shell commands. Rules are seperated into related tasks with filenames structured as `XX_task_description.smk`.

- **Custom Scripts:** In addition to Snakemake rules, custom scripts are used to complete specific tasks. These scripts, along with splitting rules into smk files, enhance modularity and readability, making it easier to manage this complex analysis workflow.

- **Configuration:** The pipeline's behavior is governed by the `snakemake/config/config.yaml` file. This configuration file specifies reference genome paths, sample details, and other configuration options, detailed below.

- **Workflow Control:** The `Snakefile` orchestrates the workflow by setting which rules to include and managing their dependencies.  
<br>

## Usage

To run the Majel pipeline:
1. Setup the slurm profile.

   ```bash
   mkdir -p ~/.config/snakemake/slurm/

   ln -s /datasets/work/hb-meth-atlas/work/pipeline_data/majel_wgbspipline/main/snakemake/config/slurm/config.yaml ~/.config/snakemake/slurm/config.yaml
   ```
2. Activate the shared Conda install. 

   ```bash
   source /datasets/work/hb-meth-atlas/work/pipeline_data/majel_conda/bin/activate

   mamba init

   mamba activate majel
   ```

   Alternatively, install the majel environment (snakemake/workflow/envs/majel.yaml) on your own Mambaforge install.  

3. Copy the snakemake directory into your own project folder (or clone Majel from the BitBucket repository).  

4. Navigate to the 'snakemake/workflow' directory and execute the following command:

   ```bash
   snakemake --use-conda [--profile <string>] [--cores <integer>] [--dryrun] [--touch] [<other snakemake options>]
                         [--config [account=<string>] [cleanup=<boolean>] [conda_prefix=<path>] [data_dir=<path>] [email=<string>]
                                   [file_prefixes=<semi-colon seperated list>] [genome=<string>] [genome_dir=<path>]
                                   [library_type=<string>] [maxins=<integer>] [non_directional=<boolean>] [output_path=<path>]
                                   [project_dir=<string OR relative path>] [roi_bed=<file path>] [rsync=<path>] [run_file=<file path>]
                                   [trim_lengths=<string, integer or semi-colon seperated list of integers] [umi_len=<integer>]
                                   [umi_loc=<string>] [umi_prefix=<string>] [whole_experiment=<boolean>]
                         ]
   ```  
<br>

## Configuration Options

Placed after the snakemake `--config` option. Default options are set in the `snakemake/config/config.yaml` file. Options are specified without dashes and use equal signs to bridge with their values (e.g. `--config option1=value1 option2=value2`):

- **account**: Specifies the project identifier for HPC job allocation. Project codes can be listed on the HPC using the get_project_codes command. (e.g., `account=OD-01234`)

- **cleanup**: This option is a boolean (true/false) flag that controls whether cleanup operations are performed after pipeline successful completion. Cleanup operations involve removing temporary files, organizing the final output and file compression. (e.g., `cleanup=True`)

- **conda_prefix**: Defines the path to the directory where Snakemake installs Conda environments. (e.g., `conda_prefix="/path/to/conda/envs`)

- **data_dir**: Specifies the directory containing input data (i.e., fastq files). (e.g., `data_dir=/path/to/input/data`)

- **email**: Provides an email address for notifications during pipeline execution. (e.g., `email=username@domain.com`)

- **file_prefixes**: Specifies a list of file name prefixes, full file names or SRAs. Used for identifying fastq files within the input data directory (specified using the `data_dir` option) or SRA files for download using sra-tools. Use semicolon (;) as the list delimiter.  (e.g., `file_prefixes=sample_A_R;sample_B_R` OR `file_prefixes=relative/path/sampleA/sample_A_R1.fq.gz;relative/path/sampleA/sample_A_R2.fq.gz` OR `file_prefixes=SRA1234567;SRA1234568`)

- **genome**: Defines the reference genome (and directory name) used for alignment and analysis. (e.g., `genome=hg38`)

- **genome_dir**: Specifies the path for reference genome directories. (e.g., `genome_dir=/path/to/reference/genome/directories`)

- **library_type**: Indicates the type of library preparation for sequencing data. Values include bs-seq, swift or em-seq. (e.g., `library_type=em-seq`)

- **maxins**: Sets the maximum insert size for paired-end reads. This determines the distance Bismark will search for alignment pairs. (e.g., `maxins=1000`)

- **non_directional**: This option is a boolean (true/false) flag that indicates whether the sequencing data is non-directional. (e.g., `non_directional=True`)

- **output_path**: Defines the directory for pipeline output. Each sample with be given a directory within this directory. (e.g., `output_path=/path/to/output/directory`)

- **project_dir**: Defines an additional directory for pipeline output reliative to the `output_path`. Used for sorting sample output into projects. (e.g. `project_dir=SRX1234567`)

- **roi_bed**: Specifies a BED file defining regions of interest (ROIs) in the genome. Used to calculate coverage and conversion statistics and defines the capture regions for Picard hybrid-selection (HS) metrics. (e.g., `roi_bed=/path/to/roi_bed.bed`)

- **rsync**: Defines the path where the final pipeline output will be transferred. (e.g. `rsync=/path/to/final/destination`)

- **run_file**: Specifies a tab-separated values (TSV) file with 3 columns for sample_name, data_dir, and file_prefixes values. Addtional columns can include any other configuration option, but must use the option=value configuration. (e.g., `run_file=/path/to/tsv/samples.tsv`)

- **trim_lengths**: Specifies the profile for number of base pairs trimmed from 5'and 3' ends of sequence reads (post adapter trimming). Defaults to value of library_type option. Options are: swift, em-seq, & no-trim, a single integer to cut from all ends, or a semi-colon seperated list of 4 integers (R1 5', R2 5', R1 3', R2 3'). (e.g., `trim_lengths=em-seq` OR `trim_lengths=10` OR `trim_lengths=10;10;15;10`)

- **umi_len**: Sets the length of Unique Molecular Identifiers (UMIs) use by Gencore for deduplication. (e.g., `umi_len=8`)

- **umi_loc**: Defines the location of UMIs within the fastq data for relocation by Fastp. Options include: index1, index2, read1, read2, per_index, per_read. (e.g., `umi_loc=per_read`)

- **umi_prefix**: Defines the prefix used by Fastp to paste in front of the UMI and by Gencore to identify the UMI. (e.g. `umi_prefix=UMI`)

- **whole_experiment**: This option is a boolean (true/false) flag that controls whether all SRAs with the same experiment accession (i.e., different runs of the same sample) will be identified in an SQL database, downloaded and processed. If specified, only one SRA needs to be provided in the of the SRA specified in the file_prefixes option. Note: causes a delay in pipeline initation resulting from the SQL database search. (e.g. `whole_experiment=True`)

Use these configuration options to customize and configure your Majel pipeline for your specific sequencing data and analysis requirements.  

Snakemake will automatically manage the creation of Conda environments and execution of rules.  
<br>

## Pipeline Rules

### Snakefile
**all**: establishes all output files from all rules  
<br>
### 00_transfer_ref_genome.smk
**transfer_ref_genome**: transfers reference genome files to working directory when using an HPC job scheduling system  
<br>
### 01_softlink_fastq.smk
**softlink_fastq**: softlinks fastq files derived from `data_dir` and `file_prefixes` `--config` options  
<br> 
### 01_sra_download.smk
**sra_download**: downloads run accessions specified in `file_prefixes` `--config` options

**sra_to_fastq**: converts sra files to fastqs  
<br>
### 02_trim_fastq.smk
**move_umi**: moves umis from read to name uisng Fastp

**trim_fastq**: trims fastqs using trim_galore

**merge_fastq**: merges r1 fastqs together and r2 fastqs together  
<br>
### 03_align_fastq.smk
**bismark_align**: aligns fastq sequences using Bismark to produce a BAM file

**sort_bam**: sorts BAM file produced by Bismark  
<br>
### 03_bismark2report.smk
**bismark_deduplicate**: used only for bismark2report: deduplicates BAM file

**bismark_methylation**: used only for bismark2report: calls cytosine methylation

**bismark2report**: produces Bismark html report  
<br>
### 04_deduplicate_bam.smk
**deduplicate_bam**: deduplicates reads in sorted BAM using Gencore for UMI support

**merge_deduplicate_bams**: merges deduplicated CT and GA BAMs  
<br>
### 05_call_methylation.smk
**call_methylation**: calls cytosine methylation using MethylDackel for downstream analysis  
<br>
### 06_call_variants.smk
**mask_converted_bases**: masks base positions in BAM potentially affected by cytosine conversion

**call_variants**: calls variants from cytosine converted data  
<br>
### 07_calculate_statistics.smk
**calculate_coverage**: calculates sequencing coverage statistics

**calculate_conversion**: calculates cytosine conversion statistics  
<br>
### 08_methylseekr_and_TDF.smk
**methylseekr_and_TDF**: calls UMRs and LMRs with and without PMDs, produces TDF file for IGV  
<br>
### 09_majel_cleanup.smk
**cleanup**: removes temporary files, restructures ouput directories, and zips text files

**rsync**: moves final output to a specified directory  
<br>

## Customization

The pipeline is designed to be customizable:

- **Adding Steps:** Extend the pipeline by adding new rules for additional analysis steps, following the pattern of existing rules.
- **Configuration:** Modify `snakemake/config/config.yaml` or use `--config` to tailor parameters, paths, and settings to your specific analysis.
- **Workflow Modification:** Adjust the `Snakefile` to control rule execution or introduce conditional logic.  
<br>

## Output

Upon successful execution, the pipeline generates various outputs:

- Alignment files
- Methylation calls
- Variant calls
- Coverage and conversion statistics
- QC files
- Log files  

Default output path is the `snakemake/workflow` directory but can be customized as needed using the  `snakemake/config/config.yaml` file or `--config` option.  
<br>

## Troubleshooting

- **Dependencies:** Ensure all required software and tools are installed. Snakemake will manage Conda environments specified in the `snakemake/worflow/envs` directory.
- **Configuration:** Double-check paths, filenames, and parameters in the `snakemake/config/config.yaml` and `run_file` file.
- **Error Handling:** Review logs generated in the `logs` directory for informative error messages in case of failures.  
<br>

## Contributing

Contributions are welcome! Feel free to suggest improvements, new features, or report issues by opening an issue or submitting a pull request.  
<br>

## License

This pipeline is released under the MIT License.
<br>

