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

configfile: os.path.join(workflow.basedir, "../config/config.yaml")
if "conda_prefix" in config:
    workflow.conda_prefix = config["conda_prefix"]

if "sample_name" not in config and "run_file" not in config: 
    sys.exit("ERROR: sample_name or run_file (TSV file with a list of sample_name and associated options) must be specified")

config["nchunks"] = {
    "chr1":5, "chr2":5, "chr3":4, "chr4":4, "chr5":3, "chr6":3, "chr7":3, "chr8":3, "chr9":2, "chr10":2,
    "chr11":2, "chr12":2, "chr13":2, "chr14":2, "chr15":2, "chr16":1, "chr17":1, "chr18":1, "chr19":1, "chr20":1, "chr21":1,
    "chr22":1, "chrX":3, "chrY":1
    }

config["genome_path"] = os.path.join(config["genome_dir"],  config["genome"])
config["file_prefixes"] = ["*"] if "file_prefixes" not in config \
                          else list(filter(None, re.split(',| |;', config["file_prefixes"])))
if "output_path" not in config: config["output_path"] = os.getcwd()

def convert_var(var):
    if isinstance(var, bool) or isinstance(var, int) or isinstance(var, list) or isinstance(var, dict):
        return var
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

def sample_config():
    D_sample_details[sample] = dict()

    # link samples with their associated files and default parameters    
    for key in config:
        D_sample_details[sample][key] = convert_var(config[key])

# initilize sample details dictionaries from sample_name or run_file
D_sample_details = dict()
if "run_file" not in config:
    sample = config["sample_name"]
    sample_config()
else:
    if not os.path.isabs(config["run_file"]):
        config["run_file"] = os.path.join(config["output_path"], config["run_file"])
    with open(config["run_file"], 'r') as run_file:
        for i, line in enumerate(run_file):
            if not line.startswith("sample_name"):
                # retreive sample names, their associated fastq or SRAs, and parameters from run_file
                sample_info = line.strip().split('\t')
                sample = sample_info[0]
                sample_config()

                # overwrites defaults and command-line --config options with sample-speficic parameters assigned in run_file
                if sample_info[1]:
                    D_sample_details[sample]["data_dir"] = sample_info[1]
                if sample_info[2]:
                    D_sample_details[sample]["file_prefixes"] = list(filter(None, re.split(',| |;', sample_info[2])))
                if len(sample_info) > 3:
                    for column in sample_info[3:]:
                        variable = column.split("=")[0]
                        value = convert_var(column.split("=")[1])
                        D_sample_details[sample][variable] = value

def print_and_write(text, output_file, mode):
    print(text)
    if output_file != F_run_info:
        with open(output_file, mode) as F_output, open(F_run_info, "a") as F_info:
            if output_file != F_run_info: 
                F_output.write(text + "\n")
            F_info.write(text + "\n")
 
num_samples = len(D_sample_details)
F_run_info = "snakemake_run_info.log"
with open(F_run_info, "w"):
    pass
message = "Running Majel snakemake pipeline on {} sample{} with the following parameters:\n"
print_and_write(message.format(num_samples, "s" if num_samples > 1 else ""), F_run_info, "a")

D_genomes = {}
D_replace_keys = dict()
D_sra_sizes = dict()
for sample in D_sample_details:
    F_sra_info_path = "{}/{}_sra_info.txt".format(D_sample_details[sample]["output_path"], sample)
    new_sample_name = sample
    D_run_files = dict()
    D_trimmed_fqs = {"_fastp_1":[], "_fastp_2":[], "_fastp_1_val_1":[], "_fastp_2_val_2":[]}
    L_run_accessions = []
    SQL_run_accessions = None
    L_sra_sizes = None
    if not "umi_len" in D_sample_details[sample]:
        D_sample_details[sample].pop("umi_prefix")
        D_sample_details[sample]["umi_loc"] = "no UMIs (umi_len not specified)"

    for file_prefix in D_sample_details[sample]["file_prefixes"]:
        if "data_dir" in D_sample_details[sample] or os.path.isabs(file_prefix):
            if "data_dir" in D_sample_details[sample]:
                sample_path = os.path.join(D_sample_details[sample]["data_dir"], os.path.dirname(file_prefix))
            else:
                sample_path = os.path.dirname(file_prefix)

            prefix_base = os.path.basename(file_prefix)
            for base_file in os.listdir(sample_path):
                match = re.match(r".*_[rR]*(1|2).*\.(f(ast)*q)\.?(gz)?$", base_file)
                if match and base_file.startswith(prefix_base):
                    D_run_files[base_file] = os.path.join(sample_path, base_file)
                    fq_pair = match.group(1)
                    stem = base_file.split('.')[0] 
                    D_trimmed_fqs["_fastp_{}".format(fq_pair)].append("{}_fastp_{}.fq.gz".format(stem, fq_pair))
                    D_trimmed_fqs["_fastp_{0}_val_{0}".format(fq_pair)].append("{0}_fastp_{1}_val_{1}.fq.gz".format(stem, fq_pair))
        else:
            SQL_run_accessions = None
            if "whole_experiment" in D_sample_details.get(sample, {}):
                if os.path.exists(F_sra_info_path):
                    with open(F_sra_info_path, 'r') as F_sra_info:
                        for line in F_sra_info:
                            columns = line.strip().split('\t')
                            experiment_accession = columns[0]
                            study_accession = columns[1]
                            SQL_run_accessions = columns[2].split(',')
                            L_sra_sizes = [int(mb) for mb in columns[3].split(',') if mb]
                else:
                    con = sqlite3.connect("/datasets/work/hb-meth-atlas/work/Data/level_2/SRAmetadb.sqlite")
                    cur = con.cursor()
                    cur.execute("SELECT experiment_accession FROM run WHERE run_accession=?", (file_prefix,))
                    experiment_accessions = cur.fetchall()
                    experiment_accession = experiment_accessions[0][0]
                    cur.execute("SELECT run_accession FROM run WHERE experiment_accession=?", (experiment_accession,))
                    SQL_run_accessions = cur.fetchall()
                    SQL_run_accessions = [run[0] for run in SQL_run_accessions]
                    cur.execute("SELECT study_accession FROM experiment WHERE experiment_accession=?", (experiment_accession,))
                    study_accessions = cur.fetchall()
                    study_accession = study_accessions[0][0]
                    with open(F_sra_info_path, 'w') as F_sra_info:
                        F_sra_info.write("{}\t{}\t{}\t".format(
                            experiment_accession, study_accession, ','.join(SQL_run_accessions)
                        ))
            
            if SQL_run_accessions:
                D_sample_details[sample]["project_dir"] = study_accession
                
                split_sample = sample.split("_")
                if split_sample[-1] == file_prefix:
                    new_sample_name = "_".join(split_sample[:-1] + [experiment_accession])
                else:
                    new_sample_name = "{}_{}".format(sample, experiment_accession)

                D_replace_keys[sample] = new_sample_name
                L_run_accessions = SQL_run_accessions
                break
            else:
                L_run_accessions.append(file_prefix)

    if len(new_sample_name.split("_")) < 4:
       sys.exit('ERROR: sample_name "{}" is not within naming convention of tissue_subtissue_healthStatus_identifier') 
    
    project_dir = D_sample_details[sample].get("project_dir", None)
    D_genomes[D_sample_details[sample]["genome"]] = D_sample_details[sample]["genome_path"]

    if not workflow.run_local:
        D_sample_details[sample]["genome_dir"] = D_sample_details[sample]["output_path"]
        D_sample_details[sample]["genome_path"] = os.path.join(D_sample_details[sample]["output_path"], D_sample_details[sample]["genome"])
        slurm = "slurm"
    else:
        slurm = ""

    if D_sample_details[sample]["rsync"]:
        if not os.path.isabs(D_sample_details[sample]["rsync"]):
            D_sample_details[sample]["rsync"] = os.path.join(D_sample_details[sample]["output_path"], D_sample_details[sample]["rsync"])
        if project_dir:
            D_sample_details[sample]["rsync"] = os.path.join(D_sample_details[sample]["rsync"], project_dir)
    
    if project_dir:
        D_sample_details[sample]["output_path"] = os.path.join(D_sample_details[sample]["output_path"], project_dir)
    sample_path = os.path.join(D_sample_details[sample]["output_path"], new_sample_name)
    os.makedirs(os.path.join(sample_path, "logs", slurm), exist_ok=True)
    os.makedirs(os.path.join(sample_path, "01_sequence_files/"), exist_ok=True)
    os.makedirs(os.path.join(sample_path, "06_call_variants/beds"), exist_ok=True)
    os.makedirs(os.path.join(sample_path, "06_call_variants/vcfs"), exist_ok=True)

    if os.path.exists(F_sra_info_path) and not SQL_run_accessions:
        with open(F_sra_info_path, 'r') as F_sra_info:
            for line in F_sra_info:
                columns = line.strip().split('\t')
                L_sra_sizes = [int(mb) for mb in columns[0].split(',') if mb]

    if L_sra_sizes:
        for i, run in enumerate(L_run_accessions):
            D_sra_sizes[run] = L_sra_sizes[i]

    for run_accession in L_run_accessions:
        D_run_files[run_accession] = os.path.join(run_accession,run_accession + ".sra")
        for label in D_trimmed_fqs:
            D_trimmed_fqs[label].append("{}{}{}.fq.gz".format(run_accession,"_1" if "1" in label else "_2",label))
        
        if run_accession not in D_sra_sizes:
            max_retries = 5
            retry_count = 0
            with open(F_sra_info_path, 'a') as F_sra_info:
                while retry_count < max_retries:
                    try:
                        sra_stat = subprocess.check_output(['sra-stat', '--xml', '--quick', run_accession])
                        break  # Exit the loop if the command succeeds
                    except subprocess.CalledProcessError as e:
                        # Handle the error or print an error message if needed
                        retry_count += 1
            
                D_sra_xml = xmltodict.parse(sra_stat)
                size_mb = int(D_sra_xml['Run']['Size']['@value']) // 1000**2
                D_sra_sizes[run_accession] = size_mb
                F_sra_info.write("{},".format(size_mb))

    D_sample_details[sample]["run_files"] = D_run_files
    D_sample_details[sample]["trimmed_fqs"] = D_trimmed_fqs

    if "trim_profile" not in D_sample_details[sample]:
        trim_profile = library_type
    else:
        trim_profile = str(D_sample_details[sample]["trim_profile"]).split(",")

    if trim_profile[0] == "em-seq":
        trim_profile = ["4","8","4","4"]
    elif trim_profile[0] in ["swift", "bs-seq"]:
        trim_profile = ["10","15","10","10"]
    elif trim_profile[0] == "no-trim":
        trim_profile[0] = "0"
    else:
        try:
            [int(x) for x in trim_profile]
        except:
            raise Exception("Invalid --trim_profile. Must be a single integer, "
                            "comma seperated list of 4 integers, or one of the "
                            "following: em-seq, swift, no-trim")

    L_trim_options = ["--clip_R1 ", "--clip_R2 ", "--three_prime_clip_R1 ", "--three_prime_clip_R2 "]

    D_sample_details[sample]["trim_lengths"] = ""
    for x in range(4):
        trim_profile.append(trim_profile[0])
        if trim_profile[x] != "0":
            D_sample_details[sample]["trim_lengths"] += L_trim_options[x] + trim_profile[x] + " "

    if "roi_bed" in D_sample_details.get(sample, {}):
        D_sample_details[sample]["nchunks"] = {"roi":1}
  
for old_key, new_key in D_replace_keys.items():
    D_sample_details[new_key] = D_sample_details.pop(old_key)

        
# Used by rule all
def get_all_files(D_genomes, D_sample_details):
    L_sample_files = []
   
    if not workflow.run_local:
        # rules/transfer_ref_genome.smk:
        # rule transfer_ref_genome outputs
        for genome in D_genomes:
            L_genome_files = [genome+'.fa', genome+'.fa.fai', genome+'.genome', 'CHH_ROI.bed',
                'Bisulfite_Genome/'+genome+'.CT_conversion.fa',
                'Bisulfite_Genome/'+genome+'.GA_conversion.fa',
                'Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa',
                'Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa'
            ]

        for file_name in L_genome_files:
            L_sample_files.append(os.path.join(config["output_path"], genome, file_name))

        os.makedirs(os.path.join(config["output_path"], genome, "logs", "slurm"), exist_ok=True)

    for sample in D_sample_details:
        sample_path = os.path.join(D_sample_details[sample]["output_path"], sample)
        
        if D_sample_details[sample]["rsync"]:
            rsync_check = "{0}/{1}/stats/methyldackel/{1}_ROI_Conversion_CpG.bedGraph.gz".format(
                D_sample_details[sample]["rsync"], sample
            )
            D_sample_details[sample]["cleanup"] = True
        else:
            rsync_check = None

        if D_sample_details[sample]["cleanup"]:
            cleanup_check = "{}/{}_sd_CpG.bedGraph.gz".format(
                sample_path, sample
            )
            D_sample_details[sample]["cleanup"] = cleanup_check
        else:
            cleanup_check = None 

        output_file = "{}/{}_majel_config_options.txt".format(sample_path, sample)
        print_and_write("---------------------\nSample: {}\n---------------------".format(sample), output_file, "w")
        for option, value in sorted(D_sample_details[sample].items()):
            print_and_write("{}: {}".format(
                option, D_sample_details[sample][option]), output_file, "a"
        )
        if "--dryrun" in sys.argv and "--touch" in sys.argv:
            print_and_write("--dryrun and --touch options specified, creating dummy files", F_run_info, "a")
        
        print_and_write("", output_file, "a")

        D_sample_details[sample]["sample_outputs"] = []
        D_sample_details[sample]["sample_outputs"].extend(
            itertools.chain.from_iterable([
                # softlink_fastq.smk:
                # rule softlink_fastq outputs
                ["{}/01_sequence_files/{}".format(sample_path,fq)
                 for fq in D_sample_details[sample]["run_files"]
                 if ".sra" not in D_sample_details[sample]["run_files"][fq]
                ],
                # sra_download.smk:
                # rule sra_download outputs
                ["{}/01_sequence_files/{}".format(sample_path, D_sample_details[sample]["run_files"][sra])
                 for sra in D_sample_details[sample]["run_files"]
                 if ".sra" in D_sample_details[sample]["run_files"][sra]
                ],
                # rule sra_to_fastq outputs
                ["{}/01_sequence_files/{}_{}.fastq.gz".format(sample_path, run_accession, read)
                 for read in ["1", "2"] for run_accession in D_sample_details[sample]["run_files"]
                 if ".sra" in D_sample_details[sample]["run_files"][run_accession]
                ],
                # trim_fastq.smk:
                # rule trim_fastq outputs
                ["{}/02_trim_fastq/{}".format(sample_path, trimmed_fq)
                    for label in D_trimmed_fqs
                    for trimmed_fq in D_sample_details[sample]["trimmed_fqs"][label]
                ],
                # rule merge_fastq outputs
                ["{}/02_trim_fastq/{}_r1.fq.gz".format(sample_path, sample),
                 "{}/02_trim_fastq/{}_r2.fq.gz".format(sample_path, sample),
                # bismark_align.smk:
                # rule bismark_align outputs
                "{}/03_align_fastq/{}_r1_bismark_bt2_pe.bam".format(sample_path, sample),
                "{}/03_align_fastq/{}_r1_bismark_bt2_PE_report.txt".format(sample_path, sample),
                "{}/03_align_fastq/{}_r1_bismark_bt2_pe.nucleotide_stats.txt".format(sample_path, sample),
                # rule sort_bam outputs
                "{}/03_align_fastq/{}_s.bam".format(sample_path, sample),
                "{}/03_align_fastq/{}_s.bam.bai".format(sample_path, sample),
                # bismark2report.smk:
                # rule bismark_deduplicate outputs
                "{}/03_align_fastq/{}_r1_bismark_bt2_pe.deduplicated.bam".format(sample_path, sample),
                "{}/03_align_fastq/{}_r1_bismark_bt2_pe.deduplication_report.txt".format(sample_path, sample),
                # rule bismark_methylation outputs
                "{}/03_align_fastq/{}_r1_bismark_bt2_pe.deduplicated.bismark.cov.gz".format(sample_path, sample),
                "{}/03_align_fastq/{}_r1_bismark_bt2_pe.deduplicated.bedGraph.gz".format(sample_path, sample),
                "{}/03_align_fastq/{}_r1_bismark_bt2_pe.deduplicated_splitting_report.txt".format(sample_path, sample),
                "{}/03_align_fastq/{}_r1_bismark_bt2_pe.deduplicated.M-bias.txt".format(sample_path, sample),
                # rule bismark2report outputs
                "{}/03_align_fastq/{}_r1_bismark_bt2_PE_report.html".format(sample_path, sample),
                # deduplicate_bam.smk:
                # rule deduplicate_bam outputs
                "{}/04_deduplicate_bam/{}_CT_genome.bam".format(sample_path, sample),
                "{}/04_deduplicate_bam/{}_GA_genome.bam".format(sample_path, sample),
                "{}/04_deduplicate_bam/{}_CT_genome_sd.bam".format(sample_path, sample),
                "{}/04_deduplicate_bam/{}_GA_genome_sd.bam".format(sample_path, sample),
                # rule merge_deduplicate_bams outputs
                "{}/04_deduplicate_bam/{}_sd.bam".format(sample_path, sample),
                "{}/04_deduplicate_bam/{}_sd.bam.bai".format(sample_path, sample),
                # call_methylation.smk:
                # rule call_methylation outputs
                "{}/05_call_methylation/{}_sd_CpG.bedGraph".format(sample_path, sample),
                ],
                ["{}/05_call_methylation/{}_ROI_Conversion_{}.bedGraph".format(sample_path,sample,mG)
                    for mG in ["CHH", "CHG", "CpG"]
                ],
                ["{}/05_call_methylation/{}_{}_{}.svg".format(sample_path, sample, context, strand)
                    for context in ["CHH", "CHG", "CpG"] for strand in ["OT", "OB"]
                ],

                ["{}/06_call_variants/{}_sd_calmd{}".format(sample_path, sample, suffix)
                    for suffix in [".bam", ".bam.bai", "_masked.bam", "_masked.bam.bai"]
                ],
                # rule call_variants outputs
                ["{}/06_call_variants/{}_sd.vcf".format(sample_path, sample)],

                # calculate_statistics.smk:
                # rule calculate_coverage outputs
                ["{}/07_calculate_statistics/{}_genomeCoverageBed.txt".format(sample_path, sample),
                # rule calculate_conversion outputs
                "{}/07_calculate_statistics/{}_sd_conversion_and_coverage.txt".format(sample_path, sample),
                ],
                # methylseekr_and_TDF.smk:
                # rule call_methylation outputs
                ["{}/08_methylseekr_and_TDF/{}_{}".format(sample_path, sample, suffix)
                    for suffix in ["PMD.bed", "UMRLMR.bed", "wPMD_UMRLMR.bed", "sd_CpG.tdf"]
                ],

                # rule cleanup
                [cleanup_check] if cleanup_check else [],
                # rule rsync
                ["{}/rsync_complete.txt".format(sample_path)] if rsync_check else []
                ]))

        # check if the sample has already been run and completed
        if rsync_check:
            D_sample_details[sample]["cleanup_inputs"] = []
            if os.path.exists(rsync_check):
                D_sample_details[sample]["sample_outputs"] = []
                print_and_write(sample + ": cleanup and rsync previously completed, nothing else to do\n", F_run_info, "a")
            elif os.path.exists(cleanup_check):
                 D_sample_details[sample]["sample_outputs"] = ["{}/rsync_complete.txt".format(sample_path)]
            else:
                D_sample_details[sample]["cleanup_inputs"] = D_sample_details[sample]["sample_outputs"][:-2]
        elif cleanup_check:
            if os.path.exists(cleanup_check):
                D_sample_details[sample]["sample_outputs"] = []
                D_sample_details[sample]["cleanup_inputs"] = []
                print_and_write(sample + ": cleanup previously completed and rsync not specified, nothing else to do\n", F_run_info, "a")
            else:
                D_sample_details[sample]["cleanup_inputs"] = D_sample_details[sample]["sample_outputs"][:-1]
        else:
            D_sample_details[sample]["cleanup_inputs"] = []
        
        L_sample_files.extend(D_sample_details[sample]["sample_outputs"])

        if "--dryrun" in sys.argv and "--touch" in sys.argv:
            for file_path in D_sample_details[sample]["sample_outputs"]:
                directory = os.path.dirname(file_path)
                os.makedirs(directory, exist_ok=True)
                with open(file_path, 'a') as dummy:
                    if "_PMD.bed" in file_path:
                        dummy.write("MethylSeekR not run, coverage < 10x")
    return L_sample_files

# allocate 2Gb of memory per thread used by a rule
def get_mem_mb(wildcards, threads):
    return threads * 2048 if threads > 1 else 4096

# split cores evenly across submitted samples
def get_cpus(min_cpus, max_cpus):
    threads = int(workflow.cores / len(D_sample_details))
    if threads > max_cpus: 
        threads = max_cpus
    elif threads < min_cpus:
        threads = min_cpus
    return threads

# retrieve file size in megabytes
def get_file_size_mb(wcs, L_file_paths):
    size_mb = 0
    for run_file in L_file_paths:
        if run_file in D_sra_sizes:
            size_mb += D_sra_sizes[run_file]
        else:
            fastq_path = D_sample_details[wcs.sample]["run_files"][run_file]
            size_mb += int(os.path.getsize(fastq_path) / 1000**2)
    return size_mb if size_mb > 0 else 1

# determine time for job submission based on starting fastq file sizes
def get_time_min(wcs, infiles, rule, threads):
    D_cpu_mins = {
        "sra_download":10,
        "sra_to_fastq":10,
        "move_umis":60,
        "bismark_deduplicate":360,
        "call_variants":480,
        "call_methylation":240,
        "bismark_methylation":480,
        "deduplicate_bam":240,
        "coverage_to_regions":60,
        "cleanup":1,
        "rsync":5,
        }

    if rule not in D_cpu_mins:
        D_cpu_mins[rule] = 30

    if rule in ["move_umis", "sra_download"]:
        # get file size of sra or fastq pair
        try:
            file_size_mb = get_file_size_mb(
                wcs, [os.path.basename(infiles.r1), os.path.basename(infiles.r2)]
            )
        except:
            file_size_mb = get_file_size_mb(
                wcs, [re.split(r"_|.sra", os.path.basename(infiles.r1))[0]]
            )
        # 2 cpu mins per GB
        time_min = int(file_size_mb/1000 * D_cpu_mins[rule] / threads)
    else:
        L_starting_run_files = [rf for rf in D_sample_details[wcs.sample]["run_files"]]
        file_size_sample = get_file_size_mb(wcs, L_starting_run_files)

        if rule == "bismark_align":
            # plotting size of sequence files versus sample processing time yielded the following linear equation:
            # y = 0.0004x + 0.8605 (for 32 cores and 64gb memory) where x is file size in kB
            # alignment time is calculated using this equation (adjusting for number of cores and adding 50% contingency)
            time_min = int((0.0004 * (1000 * file_size_sample) + 0.8605) * 32 / threads * 1.5)
        else:
            # cpu mins per GB
            time_min = int(file_size_sample/1000 * D_cpu_mins[rule] / threads)
    return time_min if time_min > 15 else 15

def get_parallel(wildcards, divider, threads):
    return threads // divider if threads > divider else 1

# Used by rule trim_fastq
# combine UMI parameters into an argument for fastp
def get_umi_info(sample):
    if "umi_len" in D_sample_details[sample]:
        return "--umi --umi_prefix={} --umi_loc={} --umi_len={} ".format(
            D_sample_details[sample]["umi_prefix"],
            D_sample_details[sample]["umi_loc"],
            D_sample_details[sample]["umi_len"]
        )
    else:
        return "" 

# Used by rule call_variants
def calculate_vcf_files(wcs,file_dir,file_ext):
    file_list = []
    for chrom in D_sample_details[wcs.sample]["nchunks"]:
        for chunk in range(1, D_sample_details[wcs.sample]["nchunks"][chrom] + 1):
            bed = "{}/{}/06_call_variants/{}/{}.{}.region.{}{}".format(
                D_sample_details[wcs.sample]["output_path"],
                wcs.sample,
                file_dir,
                D_sample_details[wcs.sample]["genome"],
                chrom,
                str(chunk),
                file_ext,
            )
            file_list.append(bed)
    return file_list

# Used by rule calculate_conversion
def conversion_estimator(sample, context):
    bedgraph_name = "{0}/{1}/05_call_methylation/{1}_ROI_Conversion_{2}.bedGraph".format(
        D_sample_details[sample]["output_path"], sample, context
    )
    with open(bedgraph_name, 'r') as F_bedgraph:
        D_ROI_conversion = {}
        for i, line in enumerate(F_bedgraph):
            if i == 0:
                continue

            L_columns = line.strip().split("\t")
            if "roi_bed" in D_sample_details[sample]:
                chromo = "ROI"
            else:
                chromo = L_columns[0]
            meth = float(L_columns[4])
            unmeth = float(L_columns[5])
            conversion = meth/(meth + unmeth)
            coverage = meth + unmeth
            if chromo in D_ROI_conversion:
                D_ROI_conversion[chromo]["conversion"].append(conversion)
                D_ROI_conversion[chromo]["coverage"].append(coverage)
            else:
                D_ROI_conversion[chromo] = {"conversion":[conversion], "coverage":[coverage]}

        L_values = []
        L_chromo = []
        for chromo in D_ROI_conversion:
            conversion = str(round(100 - mean(D_ROI_conversion[chromo]["conversion"])*100,2)) + "%"
            coverage = str(round(mean(D_ROI_conversion[chromo]["coverage"]),2))
            L_values.extend([conversion, coverage])
            L_chromo.extend([chromo + '_con', chromo + '_cov'])

    return L_values, L_chromo

# Used by onsuccess and onerror in Snakefile
# Copies the portions of snakemake.log relating to a sample to that sample's log files
def copy_log_to_sample(log):
    log_ext = '.log'
    L_lines = []
    current_sample = None
    rule_check = False

    with open(log, 'r') as log_file:
        for line in log_file:
            if "localrule all" in line:
                break 
            elif any([line.startswith(x) for x in ["MissingOutputException", "Error"]]):
                log_ext = '.error'
            
            if line.startswith("rule") or line.startswith("Error"):
                rule_check = True
                L_lines = [line]

            elif line.startswith("[") and rule_check and current_sample:
                L_lines.append("\n\n")
                with open("{0}/{1}/{1}{2}".format(D_sample_details[current_sample]["output_path"], current_sample, log_ext), 'a') as sample_log:
                    sample_log.writelines(L_lines)
                L_lines = [line]
                rule_check = False

            elif rule_check:
                L_lines.append(line)
                for sample, sample_details in D_sample_details.items():
                    if re.search(r"output:.*{}.*".format(re.escape(sample)), line):
                        current_sample = sample
                        rsync_dir = sample_details.get("rsync")
                        cur_dir = rsync_dir if rsync_dir and os.path.exists(rsync_dir) else sample_details.get("output_path")

            else:
                L_lines = [line]

    if log_ext == '.error' and current_sample:
        with open("{0}/{1}/{1}{2}".format(cur_dir, current_sample, log_ext), 'a') as sample_log:
            sample_log.writelines(L_lines)

# Used by onsuccess to perform final clean of temporary files
def remove_temp_files(sample, sample_details):
    if sample_details["cleanup"] or sample_details["rsync"]:
        for file_path in sample_details.get("cleanup_inputs", []):
            if os.path.exists(file_path):
                os.remove(file_path)
            if sample_details["rsync"]:
                rsync_file_path = file_path.replace(sample_details["output_path"], sample_details["rsync"])

        original_path = os.path.join(sample_details["output_path"], sample)
        final_path = sample_details["rsync"] if sample_details["rsync"] else original_path
        for root, dirs, files in os.walk(final_path):
            for dir_name in dirs:
                if dir_name.startswith("0"):
                    dir_path = os.path.join(root, dir_name)
                    shutil.rmtree(dir_path)
    if sample_details["rsync"] and os.path.exists(original_path):
        shutil.rmtree(original_path)
