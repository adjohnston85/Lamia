# Import required libraries and modules
import itertools
import math
import pandas as pd
import re
import shlex
import shutil
import socket
import sqlite3
import subprocess
import sys
import time
import xmltodict
from statistics import mean
from pathlib import Path
import json
import fnmatch

# Load the configuration file specified in the YAML file located relative to the workflow directory
configfile: os.path.join(workflow.basedir, "../config/config.yaml")

# Set the conda environment prefix if it exists in the configuration
if "conda_prefix" in config:
    workflow.conda_prefix = config["conda_prefix"]

# Check for the required fields in the configuration file
if "sample_name" not in config and "run_file" not in config: 
    sys.exit("ERROR: sample_name or run_file must be specified")

# Specify the number of chunks for each chromosome
config["nchunks"] = {
    "chr1":5, "chr2":5, "chr3":4, "chr4":4, "chr5":3, "chr6":3, "chr7":3, "chr8":3, "chr9":2, "chr10":2,
    "chr11":2, "chr12":2, "chr13":2, "chr14":2, "chr15":2, "chr16":1, "chr17":1, "chr18":1, "chr19":1, "chr20":1, "chr21":1,
    "chr22":1, "chrX":3, "chrY":1
    }

# Create the full path for the genome based on the genome directory and name
config["genome_path"] = os.path.join(config["genome_dir"],  config["genome"])

# Process file_prefixes from the config. If not available, set to ["*"]
config["file_prefixes"] = ["*"] if "file_prefixes" not in config else list(filter(None, re.split(' |;', config["file_prefixes"])))

# List of UMIs, assuming it's manageable to place here directly; otherwise, consider an external file
config["umi_list"] = (
    "GAGACGAT TTCCAAGG CGCATGAT ACGGAACA CGGCTAAT GCTATCCT "
    "TGGACTCT ATCCAGAG CTTAGGAC GTGCCATA TCGCTGTT TTCGTTGG "
    "AAGCACTG GTCGAAGA ACCACGAT GATTACCG GCACAACT GCGTCATT "
    "GAAGGAAG ACTGAGGT TGAAGACG GTTACGCA AGCGTGTT GATCGAGT "
    "TTGCGAAG CTGTTGAC GATGTGTG ACGTTCAG TTGCAGAC CAATGTGG "
    "ACGACTTG ACTAGGAG"
)

# If output_path is not set, use the current working directory
if "output_path" not in config: config["output_path"] = os.getcwd()

# Function to convert variable types
def convert_var(var):
    # If the variable is already of type bool, int, list, or dict, return it as is
    if isinstance(var, bool) or isinstance(var, int) or isinstance(var, list) or isinstance(var, dict):
        return var
    # If the variable is a string, try to convert it to its logical or numerical form
    elif isinstance(var, str):
        if var.lower() == "true":
            return True
        elif var.lower() == "false":
            return False
        elif var.lower() == "none":
            return None
        else:
            try:
                return int(var)
            except ValueError:
                return var

# Function to populate sample details
def sample_config():
    D_sample_details[sample] = dict()  # Initialize dictionary for a sample
    # Link samples with their associated files and default parameters
    for key in config:
        D_sample_details[sample][key] = convert_var(config[key])

# Initialize the main dictionary to hold sample details
D_sample_details = dict()

# If 'run_file' is not in the configuration, rely on 'sample_name' in config
if "run_file" not in config:
    sample = config["sample_name"]
    sample_config()
else:
    # Check if the 'run_file' path is absolute. If not, create an absolute path
    if not os.path.isabs(config["run_file"]):
        config["run_file"] = os.path.join(config["output_path"], config["run_file"])
    # Open the 'run_file' and read it line-by-line
    with open(config["run_file"], 'r') as run_file:
        for i, line in enumerate(run_file):
            if not line.startswith("sample_name"):
                # Retrieve sample names and associated details from 'run_file'
                sample_info = re.split(r'[\t,]', line.strip())
                sample = sample_info[0]
                sample_config()
                
                # Overwrite default and command-line options with sample-specific parameters from 'run_file'
                if sample_info[1]:
                    D_sample_details[sample]["data_dir"] = sample_info[1]
                if sample_info[2]:
                    D_sample_details[sample]["file_prefixes"] = list(filter(None, re.split(' |;', sample_info[2])))
                if len(sample_info) > 3:
                    for column in sample_info[3:]:
                        variable = column.split("=")[0]
                        value = convert_var(column.split("=")[1])
                        D_sample_details[sample][variable] = value

# Function to print text to the console and write it to specified files
def print_and_write(text, output_file, mode):
    print(text)  # Print the text to the console
    # Open the output file and a log file (F_run_info) in append mode
    if output_file != F_run_info:
        with open(output_file, mode) as F_output, open(F_run_info, "a") as F_info:
            if output_file != F_run_info:  # Check if output_file is not the same as the log file
                F_output.write(text + "\n")  # Write the text to the output file
            F_info.write(text + "\n")  # Write the text to the log file (F_run_info)

# Get the number of samples from the D_sample_details dictionary
num_samples = len(D_sample_details)

# Initialize the log file for the Snakemake run
F_run_info = "snakemake_run_info.log"
with open(F_run_info, "w"):  # Open and immediately close the file to ensure it's empty
    pass

message = "Running Majel snakemake pipeline on {} sample{} with the following parameters:\n"

# Display and write the message with the number of samples
print_and_write(message.format(num_samples, "s" if num_samples > 1 else ""), F_run_info, "a")

# Initialize dictionaries for genomes, key replacements, and SRA sizes
D_genomes = {}
D_replace_keys = dict()
D_sra_sizes = dict()

# Iterate through each sample in the sample details dictionary
for sample in D_sample_details:
    # Define the SRA info file path for the current sample
    F_sra_info_path = "{}/.snakemake/{}_sra_info.txt".format(D_sample_details[sample]["output_path"], sample)
    
    # Initialize variable for the new sample name 
    # Used for updating the keys of D_sample_details if the sample name is modified by this loop
    new_sample_name = sample
    
    # Initialize dictionaries and lists for run files and trimmed fastq files
    D_sample_details[sample]["fq_files"] = dict()
    D_sample_details[sample]["sra_files"] = dict()
    D_trimmed_fqs = {"_fastp_1_val_1":[], "_fastp_2_val_2":[]}
    L_run_accessions = []
    L_sra_sizes = None
    
    # Iterate through the file prefixes of the sample
    for file_prefix in D_sample_details[sample]["file_prefixes"]:
        # Determine the sample path
        if "data_dir" in D_sample_details[sample] or os.path.isabs(file_prefix):
            if "data_dir" in D_sample_details[sample]:
                sample_path = os.path.join(D_sample_details[sample]["data_dir"], os.path.dirname(file_prefix))
            else:
                sample_path = os.path.dirname(file_prefix)

            # Extract base of the file prefix
            prefix_base = os.path.basename(file_prefix)

            # Compile regex for file prefix
            prefix_regex = re.compile(rf"{prefix_base}.*")  # Adjust this regex as needed

            # Loop through files in the sample path
            for base_file in os.listdir(sample_path):
                # Use regex to find files matching a specific pattern
                match = re.match(r".*_[rR]*(1|2).*\.(f(ast)*q)\.?(gz)?$", base_file)

                # If match found, populate dictionaries
                if match and prefix_regex.match(base_file):
                    D_sample_details[sample]["fq_files"][base_file] = os.path.join(sample_path, base_file)
                    fq_pair = match.group(1)
                    stem = base_file.split('.')[0]
                    D_trimmed_fqs[f"_fastp_{fq_pair}_val_{fq_pair}"].append(f"{stem}_fastp_{fq_pair}_val_{fq_pair}.fq.gz")
        else:
            # If SRA info file exists retrive info
            if os.path.exists(F_sra_info_path):
                with open(F_sra_info_path, 'r') as F_sra_info:
                    for line in F_sra_info:
                        columns = line.strip().split(',')
                        experiment_accession = columns[0]
                        study_accession = columns[1]
                        # Populate L_run_accessions from the file
                        L_run_accessions = columns[2].split(';')
                        L_sra_sizes = [int(mb) for mb in columns[3].split(';') if mb]
            # Otherwise, if 'whole_experiment' option is specified for a sample, get SRA info from SQL database
            elif D_sample_details[sample]["whole_experiment"]:
                print_and_write('Searching SQL database for all SRA files with the same expreriment accession as ' +
                    file_prefix +
                    ' (i.e., same sample from different sequencing runs or lanes)',
                    F_run_info, "a")
                # Connect to SQLite database
                con = sqlite3.connect("/datasets/work/hb-meth-atlas/work/Data/level_2/SRAmetadb.sqlite")
                cur = con.cursor()
                
                # Retrieve experiment accession for the given run accession
                cur.execute("SELECT experiment_accession FROM run WHERE run_accession=?", (file_prefix,))
                experiment_accessions = cur.fetchall()
                experiment_accession = experiment_accessions[0][0]
                
                # Retrieve all run accessions related to the experiment
                cur.execute("SELECT run_accession FROM run WHERE experiment_accession=?", (experiment_accession,))
                SQL_run_accessions = cur.fetchall()
                L_run_accessions = [run[0] for run in SQL_run_accessions]
                
                # Retrieve study accession related to the experiment
                cur.execute("SELECT study_accession FROM experiment WHERE experiment_accession=?", (experiment_accession,))
                study_accessions = cur.fetchall()
                study_accession = study_accessions[0][0]
                
                # Write SRA info to file
                Path(F_sra_info_path).parent.mkdir(parents=True, exist_ok=True) 
                with open(F_sra_info_path, 'w') as F_sra_info:
                    F_sra_info.write("{},{},{},".format(
                        experiment_accession, study_accession, ';'.join(L_run_accessions)
                    ))
            else:
                L_run_accessions.append(file_prefix)

            # Update 'project_dir' in sample details if SQL run accessions are available
            if D_sample_details[sample]["whole_experiment"] and "project_dir" not in D_sample_details[sample]:
                D_sample_details[sample]["project_dir"] = study_accession
                
                # Generate a new sample name based on the experiment accession
                split_sample = sample.split("_")
                if split_sample[-1] == file_prefix:
                    new_sample_name = "_".join(split_sample[:-1] + [experiment_accession])
                else:
                    new_sample_name = "{}_{}".format(sample, experiment_accession)
                
                # Store new sample name for replacement
                D_replace_keys[sample] = new_sample_name
                break

    # Sort the lists of files
    D_trimmed_fqs["_fastp_1_val_1"].sort()
    D_trimmed_fqs["_fastp_2_val_2"].sort()

    # Debugging output to ensure proper pairing
    for r1, r2 in zip(D_trimmed_fqs["_fastp_1_val_1"], D_trimmed_fqs["_fastp_2_val_2"]):
        if not r1.replace('_R1', '').replace('.fq.gz', '') == r2.replace('_R2', '').replace('.fq.gz', ''):
            print(f"Warning: Mismatched pair found - R1: {r1}, R2: {r2}")
        else:
            print(f"Paired: {r1} and {r2}")

    # Check if the sample name follows the naming convention
    if len(new_sample_name.split("_")) < 4:
       sys.exit('ERROR: sample_name "{}" is not within naming convention of tissue_subtissue_healthStatus_identifier')

    # Retrieve the project directory for the sample, if available
    project_dir = D_sample_details[sample].get("project_dir", None)
    
    # Update the genomes dictionary
    D_genomes[D_sample_details[sample]["genome"]] = D_sample_details[sample]["genome_path"]

    # Configure genome directory and path based on whether the workflow is run locally
    if not workflow.run_local:
        D_sample_details[sample]["genome_dir"] = D_sample_details[sample]["output_path"]
        D_sample_details[sample]["genome_path"] = os.path.join(D_sample_details[sample]["output_path"],
            D_sample_details[sample]["genome"])
        slurm = "slurm"
    else:
        slurm = ""

    # Configure rsync path
    rsync_path = D_sample_details[sample].get("rsync")
    if rsync_path:
        if not os.path.isabs(rsync_path):
            rsync_path = os.path.join(D_sample_details[sample]["output_path"], rsync_path)
        if project_dir:
            rsync_path = os.path.join(rsync_path, project_dir)
        D_sample_details[sample]["rsync"] = rsync_path

    # Update output path if project directory is available
    if project_dir:
        D_sample_details[sample]["output_path"] = os.path.join(D_sample_details[sample]["output_path"], project_dir)

    # Generate sample path and create necessary directories
    sample_path = os.path.join(D_sample_details[sample]["output_path"], new_sample_name)
    os.makedirs(os.path.join(sample_path, "logs", slurm), exist_ok=True)

    if not os.path.exists(F_sra_info_path):
        with open(F_sra_info_path, 'w') as F_sra_info:
            F_sra_info.write("{0},{0},{1},".format(
            "unknown", ";".join(L_run_accessions)
            ))

    # Update SRA sizes if available
    if L_sra_sizes:
        for i, run in enumerate(L_run_accessions):
            D_sra_sizes[run] = L_sra_sizes[i]

    # Iterate over each run accession to set file paths and get SRA sizes
    for run_accession in L_run_accessions:
        # Update file dictionaries with SRA paths
        D_sample_details[sample]["sra_files"][run_accession] = os.path.join(run_accession, run_accession + ".sra")
        D_sample_details[sample]["fq_files"][run_accession + "_1.fastq.gz"] = os.path.join(sample_path, "01_sequence_files",
            run_accession, run_accession + "_1.fastq.gz") 
        D_sample_details[sample]["fq_files"][run_accession + "_2.fastq.gz"] = os.path.join(sample_path, "01_sequence_files",
            run_accession, run_accession + "_2.fastq.gz")

        # Update D_trimmed_fqs dictionary with fq.gz file paths
        for label in D_trimmed_fqs:
            D_trimmed_fqs[label].append("{}{}{}.fq.gz".format(run_accession, "_1" if "1" in label else "_2", label))

        # If run_accession not in D_sra_sizes, attempt to get size using 'sra-stat'
        if run_accession not in D_sra_sizes:
            max_retries = 5
            retry_count = 0
            with open(F_sra_info_path, 'a') as F_sra_info:
                while retry_count < max_retries:
                    try:
                        # Try running sra-stat to get the XML stats
                        sra_stat = subprocess.check_output(['sra-stat', '--xml', '--quick', run_accession])
                        break  # Exit the loop if the command succeeds
                    except subprocess.CalledProcessError as e:
                        # Increment retry counter upon failure
                        retry_count += 1

                # Parse the XML output to fetch the SRA size
                D_sra_xml = xmltodict.parse(sra_stat)
                size_mb = int(D_sra_xml['Run']['Size']['@value']) // 1000**2
                D_sra_sizes[run_accession] = size_mb
                F_sra_info.write("{};".format(size_mb))

    # Update D_sample_details with run files and trimmed fastq files
    D_sample_details[sample]["trimmed_fqs"] = D_trimmed_fqs

    # Determine trim profile for the sample based on library_type or trim_profile
    if "trim_profile" not in D_sample_details[sample]:
        D_sample_details[sample]["trim_profile"] = D_sample_details[sample]["library_type"]
        
    trim_profile = str(D_sample_details[sample]["trim_profile"]).split("-")

    # Set trimming profiles based on the library type or user-input trim profile
    if trim_profile[0] == "em_seq":
        trim_profile = ["10", "10", "10", "10"]
    elif trim_profile[0] == "methyl_prism":
        trim_profile = ["4", "8", "4", "8"]
    elif trim_profile[0] in ["swift", "bs_seq"]:
        trim_profile = ["10", "15", "10", "10"]
    elif trim_profile[0] == "no_trim":
        trim_profile[0] = "0"
    else:
        # Validate if trim_profile consists of integers
        try:
            [int(x) for x in trim_profile]
        except:
            raise Exception(str(trim_profile) + " is an invalid --trim_profile. Must be a single integer, "
                            "dash (-) separated list of 4 integers, or one of the "
                            "following: em_seq, methyl_prism, swift, no_trim")
    
    # List of Trim_galore trimming options
    L_trim_options = ["--clip_R1 ", "--clip_R2 ", "--three_prime_clip_R1 ", "--three_prime_clip_R2 "]

    # Initialize the 'trim_lengths' attribute in D_sample_details for the current sample
    D_sample_details[sample]["trim_lengths"] = ""
    # Loop through to set the trim lengths for all four possible options (two reads, each with two trim sides)
    for x in range(4):
        trim_profile.append(trim_profile[0])
        if trim_profile[x] != "0":
            D_sample_details[sample]["trim_lengths"] += L_trim_options[x] + trim_profile[x] + " "

    # Check if 'roi_bed' (Region of Interest in BED format) is in D_sample_details
    if "roi_bed" in D_sample_details.get(sample, {}):
        D_sample_details[sample]["nchunks"] = {"roi": 1}

    # If 'maxins' (maximum insert size) is not set, set it based on the library type
    if "maxins" not in D_sample_details[sample]:
        if D_sample_details[sample]["library_type"] == "em-seq":
            D_sample_details[sample]["maxins"] = 1000
        else:
            D_sample_details[sample]["maxins"] = 500

# Loop over the keys to be replaced in D_sample_details
for old_key, new_key in D_replace_keys.items():
    # Replace the old key with the new key, while keeping the associated values unchanged
    D_sample_details[new_key] = D_sample_details.pop(old_key)

def generate_sample_outputs(sample_path, sample, details, cleanup_check, rsync_check):
    outputs = []

    # Section 1: Handling sequence files (FASTQ and SRA)
    sequence_files = [
        "{}/01_sequence_files/{}".format(sample_path, fq)
        for fq in details.get("fq_files", [])
    ] + [
        "{}/01_sequence_files/{}".format(sample_path, details["sra_files"][sra])
        for sra in details.get("sra_files", [])
    ] + [
        "{0}/01_sequence_files/{1}/{1}_{2}.fastq.gz".format(sample_path, run_accession, read)
        for read in ["1", "2"] for run_accession in details.get("sra_files", [])
    ]
    outputs.extend(sequence_files)

    # Section 2: Trimming FASTQ files
    trimmed_files = [
        "{}/02_trim_fastq/{}".format(sample_path, trimmed_fq)
        for label in details.get("trimmed_fqs", {})
        for trimmed_fq in details["trimmed_fqs"].get(label, [])
    ] + [
        "{}/02_trim_fastq/{}_r1.fq.gz".format(sample_path, sample),
        "{}/02_trim_fastq/{}_r2.fq.gz".format(sample_path, sample)
    ]
    outputs.extend(trimmed_files)

    # Section 3: Alignment with Bismark
    alignment_files = [
        "{}/03_align_fastq/{}_r1_bismark_bt2_pe.bam".format(sample_path, sample),
#        "{}/03_align_fastq/{}_r1_bismark_bt2_pe.nonCG_removed_seqs.bam".format(sample_path, sample),
#        "{}/03_align_fastq/{}_r1_bismark_bt2_pe.non-conversion_filtering.txt".format(sample_path, sample),
#:wq
#        "{}/03_align_fastq/{}_r1_bismark_bt2_pe.nonCG_filtered.bam".format(sample_path, sample),
        "{}/03_align_fastq/{}_r1_bismark_bt2_PE_report.txt".format(sample_path, sample),
        "{}/03_align_fastq/{}_r1_bismark_bt2_pe.nucleotide_stats.txt".format(sample_path, sample),
        "{}/03_align_fastq/{}_s.bam".format(sample_path, sample),
        "{}/03_align_fastq/{}_s.bam.bai".format(sample_path, sample),
        "{}/03_align_fastq/{}_insert_size_histogram.pdf".format(sample_path, sample),
        "{}/03_align_fastq/{}_insert_size_metrics.txt".format(sample_path, sample),
        "{}/03_align_fastq/{}_hs_metrics.txt".format(sample_path, sample)
    ]
    outputs.extend(alignment_files)

    # Section 4: BAM deduplication and further processing
    dedup_files = [
#       "{}/04_deduplicate_bam/{}_CT_genome.bam".format(sample_path, sample),
#       "{}/04_deduplicate_bam/{}_GA_genome.bam".format(sample_path, sample),
#       "{}/04_deduplicate_bam/{}_CT_genome_reverted.bam".format(sample_path, sample),
#       "{}/04_deduplicate_bam/{}_GA_genome_reverted.bam".format(sample_path, sample),
#        "{}/04_deduplicate_bam/{}_CT_genome_masked.bam".format(sample_path, sample),
#        "{}/04_deduplicate_bam/{}_GA_genome_masked.bam".format(sample_path, sample),
#       "{}/04_deduplicate_bam/{}_reverted_s.bam".format(sample_path, sample),
#       "{}/04_deduplicate_bam/{}_reverted_sd.bam".format(sample_path, sample),
#        "{}/04_deduplicate_bam/{}_masked_s.bam".format(sample_path, sample),
#        "{}/04_deduplicate_bam/{}_masked_sd.bam".format(sample_path, sample),
#       "{}/04_deduplicate_bam/{}_sd_flagstat.txt".format(sample_path, sample),
#       "{}/04_deduplicate_bam/{}_gencore.json".format(sample_path, sample),
#       "{}/04_deduplicate_bam/{}_gencore.html".format(sample_path, sample),
    ]
    outputs.extend(dedup_files)

    # Section 5: Methylation calling and visualization
    methylation_files = [
        "{}/05_call_methylation/{}_sd_CpG.bedGraph".format(sample_path, sample),
        "{}/05_call_methylation/{}_CHH_mbias.txt".format(sample_path, sample),
        "{}/05_call_methylation/{}_CHG_mbias.txt".format(sample_path, sample),
        "{}/05_call_methylation/{}_CpG_mbias.txt".format(sample_path, sample)
    ] + [
        "{}/05_call_methylation/{}{}_ROI_Conversion_{}.bedGraph".format(sample_path, sample, suffix, mG)
        for mG in ["CHH", "CHG", "CpG"]
        for suffix in ["_sd", "_s"]
    ] + [
        "{}/05_call_methylation/{}_{}_{}.svg".format(sample_path, sample, context, strand)
        for context in ["CHH", "CHG", "CpG"]
        for strand in ["OT", "OB"]
    ]
    outputs.extend(methylation_files)

    # Section 6: Variant calling
    if details.get("call_variants"):
        variant_files = [
            "{}/06_call_variants/{}_sd_calmd{}".format(sample_path, sample, suffix)
            for suffix in [".bam", ".bam.bai", "_masked.bam", "_masked.bam.bai"]
        ] + [
            "{}/06_call_variants/{}_sd.vcf".format(sample_path, sample)
        ]
        outputs.extend(variant_files)

    # Section 7: Calculate statistics
    stats_files = [
        "{}/07_calculate_statistics/{}_s_genomeCoverageBed.txt".format(sample_path, sample),
        "{}/07_calculate_statistics/{}_sd_genomeCoverageBed.txt".format(sample_path, sample),
        "{}/07_calculate_statistics/{}_s_conversion_and_coverage.txt".format(sample_path, sample),
        "{}/07_calculate_statistics/{}_sd_conversion_and_coverage.txt".format(sample_path, sample),
        "{}/07_calculate_statistics/{}_r1_bismark_bt2_pe.deduplicated_splitting_report.txt".format(sample_path, sample),
#        "{}/07_calculate_statistics/{}_r1_bismark_bt2_pe.deduplication_report.txt".format(sample_path, sample),
        "{}/07_calculate_statistics/{}_r1_bismark_bt2_PE_report.html".format(sample_path, sample)
    ]
    outputs.extend(stats_files)

    # Section 8: Further methylation analysis
    methyl_analysis_files = [
        "{}/08_methylseekr_and_TDF/{}_{}".format(sample_path, sample, suffix)
        for suffix in ["PMD.bed", "UMRLMR.bed", "wPMD_UMRLMR.bed", "sd_CpG.tdf"]
    ]

    outputs.extend(methyl_analysis_files)

    # Section 9: Cleanup and finalization
    # Outputs from the 'cleanup' rule in 09_majel_cleanup.smk
    #if cleanup_check:
    #    outputs.append(cleanup_check)  # Conditional output based on cleanup_check flag

    # Outputs from the 'rsync' rule
    #if rsync_check:
    #    outputs.append("{}/rsync_complete.txt".format(sample_path)) # Conditional output based on rsync_check flag

    return outputs

# Function to get all files associated with each sample and genome
# This function is used by rule 'all' in Snakemake workflow
def get_all_files(D_genomes, D_sample_details):
    # List to store sample files
    L_sample_files = []

    # If the workflow is not being run locally
    if not workflow.run_local:
        # Create a list of genome-related files for each genome
        for genome in D_genomes:
            L_genome_files = [genome+'.fa', genome+'.fa.fai', genome+'.genome', 'CHH_ROI.bed',
                'Bisulfite_Genome/'+genome+'.CT_conversion.fa',
                'Bisulfite_Genome/'+genome+'.GA_conversion.fa',
                'Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa',
                'Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa'
            ]

        # Add genome-related files to the sample files list
        for file_name in L_genome_files:
            L_sample_files.append(os.path.join(config["output_path"], genome, file_name))

        # Create log directories
        os.makedirs(os.path.join(config["output_path"], genome, "logs", "slurm"), exist_ok=True)

    # Loop through each sample in sample details
    for sample, details in D_sample_details.items():
        # Define sample-specific path
        sample_path = os.path.join(details["output_path"], sample)

        # Check if rsync is specified and create rsync_check accordingly
        rsync_path = details.get("rsync")
        if rsync_path:
            rsync_check = "{0}/{1}/stats/methyldackel/{1}_ROI_Conversion_CpG.bedGraph.gz".format(
                details["rsync"], sample
            )
            # Enable cleanup for the sample
            details["cleanup"] = True
        else:
            rsync_check = None

        # Check if cleanup is specified and create cleanup_check accordingly
        if details["cleanup"]:
            cleanup_check = "{}/{}_sd_CpG.bedGraph.gz".format(
                sample_path, sample
            )
            details["cleanup"] = cleanup_check
        else:
            cleanup_check = None

        # Define the path for the output file that contains configuration options
        output_file = "{}/{}_majel_config_options.txt".format(sample_path, sample)

        # Write sample name and details to the output file
        print_and_write("---------------------\nSample: {}\n---------------------".format(sample), output_file, "w")
        for option, value in sorted(details.items()):
            print_and_write("{}: {}".format(
                option, details[option]), output_file, "a"
            )

        # Handle the case where both --dryrun and --touch command-line options are specified
        if "--dryrun" in sys.argv and "--touch" in sys.argv:
            print_and_write("--dryrun and --touch options specified, creating dummy files", F_run_info, "a")

        # Write an empty line to the output file
        print_and_write("", output_file, "a")
		
        # Initialize an empty list for sample outputs
        details["sample_outputs"] = generate_sample_outputs(sample_path, sample, details, cleanup_check, rsync_check)

        # Check if the sample processing has already been completed
        if rsync_check:
            details["cleanup_inputs"] = []

            # Check if rsync has been completed
            if os.path.exists(rsync_check):
                # If rsync is done, clear sample outputs and log the information
                details["sample_outputs"] = []
                print_and_write(sample + ": cleanup and rsync previously completed, nothing else to do\n", F_run_info, "a")
        
        elif cleanup_check:
            # Check if cleanup has been completed
            if os.path.exists(cleanup_check):
                # If cleanup is done and rsync not specified, clear sample outputs and inputs, and log the information
                details["sample_outputs"] = []
                details["cleanup_inputs"] = []
                print_and_write(sample + ": cleanup previously completed and rsync not specified, nothing else to do\n", F_run_info, "a")
        else:
            # Clear cleanup inputs if neither cleanup nor rsync are specified
            details["cleanup_inputs"] = []
        
        # Extend the list of sample files with the outputs for this sample
        L_sample_files.extend(details["sample_outputs"])
        
        # Dry run option: Create empty files as placeholders for the expected outputs
        if "--dryrun" in sys.argv and "--touch" in sys.argv:
            for file_path in details["sample_outputs"]:
                # Create directories if they don't exist
                directory = os.path.dirname(file_path)
                os.makedirs(directory, exist_ok=True)
                
                # Create or touch the file
                with open(file_path, 'a') as dummy:
                    # Special case: add a note in the _PMD.bed file
                    if "_PMD.bed" in file_path:
                        dummy.write("MethylSeekR not run, coverage < 10x")
        
    return L_sample_files

# Function to allocate memory for each thread used by a rule
def get_mem_mb(wildcards, threads, mem_per_core):
    if threads == 1:
        return threads * mem_per_core * 2
    else:
        return threads * mem_per_core

# Function to determine the number of CPU cores to allocate to each sample
def get_cpus(min_cpus, max_cpus):
    # Calculate number of threads by dividing total available cores by the number of samples
    threads = int(workflow.cores / len(D_sample_details))
    
    # Limit the threads to the maximum allowed CPUs
    if threads > max_cpus:
        threads = max_cpus
    # Make sure the threads are at least the minimum required CPUs
    elif threads < min_cpus:
        threads = min_cpus
    
    return threads

# Function to get the file size in megabytes for a given list of file paths
def get_file_size_mb(wcs, L_file_paths):
    # Initialize size to 0 MB
    size_mb = 0
    
    # Loop through each file in the list
    for run_file in L_file_paths:
        # If the file is an SRA file, retrieve its size from the D_sra_sizes dictionary
        if run_file in D_sra_sizes:
            size_mb += D_sra_sizes[run_file]
        else:
            # If the file exists on the system, retrieve its size and convert to megabytes
            if os.path.exists(run_file):
                size_mb += int(os.path.getsize(run_file) / 1024**2)
                
    # Return size in MB, with a minimum size of 1 MB if the file doesn't exist or is empty
    return size_mb if size_mb > 0 else 1

# Function to estimate the time needed for job submission based on the size of the starting fastq files
def get_time_min(wcs, infiles, rule, threads):
    # Dictionary mapping rule names to baseline CPU minutes required for each rule
    D_cpu_mins = {
        "sra_download": 10,
        "sra_to_fastq": 500,
        "trim_fastq": 400,
        "bismark_align": 1400,
        "sort_bam": 45,
        "bismark_deduplicate": 50,
        "bismark_methylation": 30,
        "mask_converted_bases":600,
        "call_variants": 200,
        "call_methylation": 400,
        "bismark_methylation": 600,
        "deduplicate_bam": 5000,
        "split_deaminated_genomes": 500,
        "revert_deamination": 15000,
        "merge_sort_reverted": 50,
        "cleanup": 10,
        "rsync": 50,
    }

    # Default time if rule not in D_cpu_mins
    if rule not in D_cpu_mins:
        D_cpu_mins[rule] = 150

    # Initialize thread modifier for time calculation
    threads_modifier = threads
    diminisher = 0.0119047619047619
    # Modify threads_modifier based on the number of threads
    for x in range(0, threads):
        threads_modifier = threads_modifier - diminisher * x

    # Estimate time for rules prior to fastq merging
    if rule in ["move_umis", "sra_download", "sra_to_fastq", "trim_fastq"]:
        # If the input is an SRA file, get its size
        if "." not in infiles.r1:
            file_size_mb = get_file_size_mb(
                wcs, [infiles.r1]
            )
        else:
            # Get the file sizes of both fastq files
            L_file_paths = [infiles.r1, infiles.r2]
            file_size_mb = get_file_size_mb(wcs, L_file_paths)
    else:
        # For other rules, calculate the size based on all the starting run files
        L_starting_run_files = [D_sample_details[wcs.sample]["fq_files"][rf] for rf in D_sample_details[wcs.sample]["fq_files"]]
        file_size_mb = get_file_size_mb(wcs, L_starting_run_files)

    # Calculate estimated time in minutes per GB of data
    time_min = int(file_size_mb / 1000 * D_cpu_mins[rule] / int(threads_modifier))
#    print("Rule: " + rule + " file: " + str(infiles) + " size_mb: " + str(file_size_mb) + 
#        " threads_mod: " + str(threads_modifier) + " cpu_mins: " + str(D_cpu_mins[rule]) +
#        " time_mins: " + str(time_min))

    # Set minimum and maximum time limits for the job
    if rule == "bismark_align" and time_min < 30:
        return 30
    elif time_min < 5:
        return 5
    elif time_min > 10080:  # Maximum time limit is 7 days in minutes
        return 10080
    else:
        return time_min

# Function to determine how many threads each subprocess should get when running in parallel.
def get_parallel(wildcards, divider, threads):
    return threads // divider if threads > divider else 1

# Function used in the 'trim_fastq' rule to fetch UMI (Unique Molecular Identifier) parameters
# and construct a command-line argument string for 'fastp' tool.
def get_umi_info(sample):
    if "umi_len" in D_sample_details[sample]:  # Check if 'umi_len' exists in the sample details dictionary
        return "--umi --umi_prefix={} --umi_loc={} --umi_len={} ".format(
            D_sample_details[sample]["umi_prefix"],
            D_sample_details[sample]["umi_loc"],
            D_sample_details[sample]["umi_len"]
        )
    else:
        return ""  # Return an empty string if 'umi_len' does not exist

# Function to get the UMI prefix for a given sample.
# This function can be used by any rule that requires the UMI prefix.
def get_umi_prefix(sample):
    if "umi_prefix" in D_sample_details[sample]:  # Check if 'umi_prefix' exists in the sample details
        return "--umi_prefix={} ".format(
            D_sample_details[sample]["umi_prefix"],
        )
    else:
        return ""  # Return an empty string if 'umi_prefix' does not exist

# Function used by the 'call_variants' rule to construct a list of VCF (Variant Call Format) files.
# It iterates through chromosomes and their chunks to generate the full path for each file.
def calculate_vcf_files(wcs, file_dir, file_ext):
    file_list = []  # Initialize an empty list to store the file paths
    for chrom in D_sample_details[wcs.sample]["nchunks"]:
        for chunk in range(1, D_sample_details[wcs.sample]["nchunks"][chrom] + 1):
            # Construct the full path for each VCF file based on the chromosome and chunk
            bed = "{}/{}/06_call_variants/{}/{}.{}.region.{}{}".format(
                D_sample_details[wcs.sample]["output_path"],
                wcs.sample,
                file_dir,
                D_sample_details[wcs.sample]["genome"],
                chrom,
                str(chunk),
                file_ext,
            )
            file_list.append(bed)  # Append the path to the list
    return file_list  # Return the list of VCF file paths

# Transform MethylDackel bedgraph output into format readable by bismark2report
def transform_methylation(infile, outfile, D_conversion):
    total_Cs = sum(sum(value) for value in D_conversion.values())
    CpG_meth = D_conversion["CpG"][0]
    CHG_meth = D_conversion["CHG"][0]
    CHH_meth = D_conversion["CHH"][0]

    CpG_CtoT = D_conversion["CpG"][1]
    CHG_CtoT = D_conversion["CHG"][1]
    CHH_CtoT = D_conversion["CHH"][1]

    CpG_perc = format(CpG_meth/(CpG_meth+CpG_CtoT) * 100, '.1f')
    CHG_perc = format(CHG_meth/(CHG_meth+CHG_CtoT) * 100, '.1f')
    CHH_perc = format(CHH_meth/(CHH_meth+CHH_CtoT) * 100, '.1f')

    with open(outfile, 'w') as f:
        version = "MethylDackel"
        f.write(infile + "\n\nParameters used to extract methylation information:\n" + version +
            "\n\nFinal Cytosine Methylation Report\n=================================" +
            "\nTotal number of C's analysed:\t" + str(total_Cs) +
            "\n\nTotal methylated C's in CpG context:\t" + str(CpG_meth) +
            "\nTotal methylated C's in CHG context:\t" + str(CHG_meth) +
            "\nTotal methylated C's in CHH context:\t" + str(CHH_meth) +
            "\n\nTotal C to T conversions in CpG context:\t" + str(CpG_CtoT) +
            "\nTotal C to T conversions in CHG context:\t" + str(CHG_CtoT) +
            "\nTotal C to T conversions in CHH context:\t" + str(CHH_CtoT) +
            "\n\nC methylated in CpG context:\t" + CpG_perc +
            "%\nC methylated in CHG context:\t" + CHG_perc +
            "%\nC methylated in CHH context:\t" + CHH_perc + "%"
        )

# Create a text file to match the output of bismark deduplication for bismark2report
def process_dedup_json(infile, outfile, bam):
    # Read the JSON from the input file
    with open(infile, 'r') as f:
        data = json.load(f)

    # Extract necessary data from JSON
    total_reads = data["before_processing"]["total_reads"]
    duplication_histogram = data["before_processing"]["duplication_level_histogram"]

    # Do calculations
    alignments_analyzed = int(total_reads / 2)
    duplicated_alignments_removed = sum([duplication_histogram[i] * i for i in range(len(duplication_histogram))])
    percentage_duplicated_removed = (duplicated_alignments_removed / alignments_analyzed) * 100
    different_positions = sum(duplication_histogram[1:])
    leftover_sequences = alignments_analyzed - duplicated_alignments_removed
    percentage_leftover = (leftover_sequences / alignments_analyzed) * 100

    # Create the results as a list of strings
    results = [
        "\nTotal number of alignments analysed in {}:\t{}".format(bam, alignments_analyzed),
        "Total number duplicated alignments removed:\t{} ({:.2f}%)".format(duplicated_alignments_removed, percentage_duplicated_removed),
        "Duplicated alignments were found at:\t{} different position(s)\n".format(different_positions),
        "Total count of deduplicated leftover sequences:\t{} ({:.2f}% of total)".format(leftover_sequences, percentage_leftover)
    ]

    # Write the results to the output file
    with open(outfile, 'w') as f:
        for line in results:
            f.write(line + '\n')

# Transform MethylDackel mbias output into format readable by bismark2report 
def transform_mbias(infile, outfile, context, w_option):
    data = {}

    with open(infile, 'r') as f:
        # Skip the header
        next(f)
        
        for line in f:
            fields = line.strip().split()
            strand, read, position = fields[0], int(fields[1]), int(fields[2])
            nMethylated, nUnmethylated = int(fields[3]), int(fields[4])

            if read not in data:
                data[read] = {}

            if position not in data[read]:
                data[read][position] = {"methylated": 0, "unmethylated": 0}

            data[read][position]["methylated"] += nMethylated
            data[read][position]["unmethylated"] += nUnmethylated

    with open(outfile, w_option) as out:
        for read, positions in data.items():
            out.write("{} context (R{})\n".format(context, read))
            out.write("=" * (15 + len(str(read))) + "\n")
            out.write("position\tcount methylated\tcount unmethylated\t% methylation\tcoverage\n")
            
            for pos, counts in sorted(positions.items()):
                total = counts["methylated"] + counts["unmethylated"]
                methylation = (counts["methylated"] / total) * 100

                out.write("{}\t{}\t{}\t{:.2f}\t{}\n".format(pos, counts['methylated'], counts['unmethylated'], methylation, total))

            out.write("\n")

# Used by rule calculate_conversion to calculate the conversion rates
# from a bedGraph file.
def conversion_estimator(sample, context, D_conversion, suffix):
    # Constructs the full path to the bedGraph file based on the sample and context.
    bedgraph_name = "{0}/{1}/05_call_methylation/{1}{2}_ROI_Conversion_{3}.bedGraph".format(
        D_sample_details[sample]["output_path"], sample, suffix, context
    )
    # Open the bedGraph file for reading.
    with open(bedgraph_name, 'r') as F_bedgraph:
        # Initialize a dictionary to hold the conversion rates, mean_C_coverage, and total_Cs for each chromo (or region of interest).
        D_ROI_conversion = {}

        # Loop over each line in the bedGraph file.
        for i, line in enumerate(F_bedgraph):
            # Skip the first line, which is usually the header.
            if i == 0:
                continue

            # Split the line by tab characters into a list of columns.
            L_columns = line.strip().split("\t")

            # Determine the chromo (chromosome or region of interest) based on whether an 'roi_bed' is specified.
            if "roi_bed" in D_sample_details[sample]:
                chromo = "ROI"
            else:
                chromo = L_columns[0]

            # Extract methylation and unmethylation counts and calculate conversion and coverage.
            meth = float(L_columns[4])
            unmeth = float(L_columns[5])
            total_Cs = meth + unmeth
            conversion = meth / total_Cs
            D_conversion[context][0] += int(meth)
            D_conversion[context][1] += int(unmeth)

            # Update the dictionary with the new conversion, mean_C_coverage, and total_Cs data.
            if chromo in D_ROI_conversion:
                D_ROI_conversion[chromo]["conversion"][0] += unmeth
                D_ROI_conversion[chromo]["conversion"][1] += total_Cs
                D_ROI_conversion[chromo]["mean_C_coverage"].append(total_Cs)
                D_ROI_conversion[chromo]["total_Cs"] += total_Cs
            else:
                D_ROI_conversion[chromo] = {"conversion": [unmeth, total_Cs], "mean_C_coverage": [total_Cs], "total_Cs": total_Cs}

        # Initialize lists to hold the return values.
        L_values = []
        L_chromo = []

        # Loop over each chromo to calculate the mean conversion, mean_C_coverage, and total_Cs and populate the return lists.
        for chromo in D_ROI_conversion:
            unmeth = D_ROI_conversion[chromo]["conversion"][0]
            total_c = D_ROI_conversion[chromo]["conversion"][1]
            conversion = str(round(unmeth / total_c * 100, 2)) + "%"
            mean_C_coverage = str(round(mean(D_ROI_conversion[chromo]["mean_C_coverage"]), 2))
            total_Cs = D_ROI_conversion[chromo]["total_Cs"]
            L_values.extend([conversion, mean_C_coverage, str(total_Cs)])
            L_chromo.extend([chromo + '_con', chromo + '_mean_C_cov', chromo + '_total_Cs'])

    # Return the lists of conversion, mean_C_coverage, total_Cs values and their identifiers, along with the updated D_conversion dictionary.
    return [L_values, L_chromo], D_conversion


# Used by onsuccess and onerror in Snakefile to copy relevant log information
# to each sample's log files.
def copy_log_to_sample(log):
    # Initialize log extension and variables
    log_ext = '.log'
    L_lines = []
    current_sample = None
    rule_check = False

    # Open and read the main log file
    with open(log, 'r') as log_file:
        for line in log_file:
            # Stop reading if the local rule is encountered
            if "localrule all" in line:
                break
            # Update log extension to '.error' if an error or missing output is found
            elif any([line.startswith(x) for x in ["MissingOutputException", "Error"]]):
                log_ext = '.error'
            
            # Check for lines that start with 'rule' or 'Error' and update state variables
            if line.startswith("rule") or line.startswith("Error"):
                rule_check = True
                L_lines = [line]
                
            # Write to sample log if a rule block has been encountered
            elif line.startswith("[") and rule_check and current_sample:
                L_lines.append("\n\n")
                with open("{0}/{1}/{1}{2}".format(D_sample_details[current_sample]["output_path"], current_sample, log_ext), 'a') as sample_log:
                    sample_log.writelines(L_lines)
                L_lines = [line]
                rule_check = False
                
            # Gather lines related to the same rule or error
            elif rule_check:
                L_lines.append(line)
                for sample, sample_details in D_sample_details.items():
                    if re.search(r"output:.*{}.*".format(re.escape(sample)), line):
                        current_sample = sample
                        rsync_dir = sample_details.get("rsync")
                        if rsync_dir:
                            cur_dir = rsync_dir if os.path.exists(rsync_dir) else sample_details.get("output_path")
                        else:
                            cur_dir = sample_details.get("output_path")
                            
            # Default: Reset lines
            else:
                L_lines = [line]
                
    # Write to the sample-specific log file if an error has occurred
    if log_ext == '.error' and current_sample:
        with open("{0}/{1}/{1}{2}".format(cur_dir, current_sample, log_ext), 'a') as sample_log:
            sample_log.writelines(L_lines)

# Used by onsuccess in the Snakefile to remove temporary files and folders
# after the pipeline is finished running.
def remove_temp_files(sample, sample_details):
    # Check if cleanup or rsync options are enabled for this sample
    if sample_details.get("cleanup", False) or "rsync" in sample_details:
        # Iterate through the list of temporary input files to be cleaned
        for file_path in sample_details.get("cleanup_inputs", []):
            # Remove the file if it exists
            if os.path.exists(file_path):
                os.remove(file_path)
            # Replace output path with rsync path, if rsync option is enabled
            if sample_details["rsync"]:
                rsync_file_path = file_path.replace(sample_details["output_path"], sample_details["rsync"])

        # Define original and final output paths
        original_path = os.path.join(sample_details["output_path"], sample)
        final_path = sample_details["rsync"] if sample_details["rsync"] else original_path

        # Walk through the directory tree starting at final_path
        for root, dirs, files in os.walk(final_path):
            # Remove directories that start with "0"
            for dir_name in dirs:
                if dir_name.startswith("0"):
                    dir_path = os.path.join(root, dir_name)
                    shutil.rmtree(dir_path)
    
    # Remove original output path if rsync option is enabled and path exists
    if "rsync" in sample_details and os.path.exists(original_path):
        shutil.rmtree(original_path)
