#!usr/bin/env python

from ruffus import *
import ruffus.cmdline as cmdline
from ruffus.cmdline import MESSAGE
import os
from subprocess import Popen, PIPE
import time
import shlex
import re
import pandas as pd
import sys
from pathlib import Path
import math

#define all paths to sofware (Popen will not use $PATH)
def getFunctionPath(target):
    proc = Popen(shlex.split('which %s' % target), stdout = PIPE, stderr = PIPE)
    functionPath = re.sub("\\n", "",proc.communicate()[0].decode("utf-8"))
    return functionPath


# get important directory locations
pipeline_path = os.path.dirname(os.path.abspath(__file__))

# Parse command line arguments
parser = cmdline.get_argparse(description="""
                                        Majel.py - Automated WGBS processing pipeline for
                                        sra and fastq file types
                                        """)

# Add command line arguments
#removing multiple alingers for the time being
#parser.add_argument("--aligner", help="Select prefered aligner [bismark, bwameth]. Defaults to bismark",
#                          default="bismark")
parser.add_argument("--genome", help="Genome for alignment. Must be a directory in your --genome_path. Defaults to hg38 (hg38)",
                          default="hg38")
parser.add_argument("--data_dir", help="Directory for input fastq/sra files",
                    default = "./data/")
parser.add_argument("--sample_name", help='Sample name. Used for output file names and should match directory name. Will determine colouring of TDF track file and must conform to following convention "Tissue_SubTissue_HealthStatus_Identifier"')
parser.add_argument("--threads", help="Speed up alignment and other processes by increasing number of threads. Defaults to 20 (this number is divided by 5 for Bismark, as 4 aligner threads uses ~20 cores and ~40GB of RAM)",
                          default=4)
parser.add_argument("--trim_profile", help="Sets the profile for number of base pairs trimmed from 3' and 5' ends of sequence reads (post adapter trimming). Options are: bs-seq, em-seq, & no-trim (defaults to bs-seq)",
                          default="bs-seq")
parser.add_argument("--pbat", action='store_true',
                          help="Specify when aligning pbat library")
parser.add_argument("--non_directional", action='store_true',
                          help="Specify to instruct Bismark to use all four alignment outputs")
parser.add_argument("--is_paired_end",help = "Is the libarary paired end (defaults to True)",
                          default="True")
parser.add_argument("--genome_path",help = "Path to genome folder, must contain a --genome directory",
                    default="/datasets/work/hb-meth-atlas/work/pipeline_data/Genomes/")
options = parser.parse_args()

# standard python logger which can be synchronised across concurrent Ruffus tasks
logger, logger_mutex = cmdline.setup_logging(__name__, options.log_file, options.verbose)
options.logger = logger
#fixing aligner to bismark but keeping old code for when walt is added
options.aligner = "bismark"
# Utility function
class JobFailException(Exception):  # For when a cluster job fails
     def __init__(self, value):
          self.parameter = value

     def __str__(self):
          return repr(self.parameter)
     

# Prepend a timestamp to a debug message
def timestamp(msg):
     return time.strftime("%c") + ": " + msg


def make_target_dir(path):
     if os.path.exists(path):
          pass
     else:
          os.mkdir(path)
          

def execute_cmd(cmd, cwd=None):
     """
     Execute the external command and get its exitcode, stdout and stderr.
     """
     args = shlex.split(cmd)
     proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=cwd)
     out, err = proc.communicate()
     exitcode = proc.returncode
     return exitcode, out, err


def checkBam(bam_prefix):
     args = shlex.split('samtools quickcheck %s.bam' % bam_prefix)
     proc = Popen(args, stdout=PIPE, stderr=PIPE)
     out, err = proc.communicate()
     exitcode = proc.returncode
     if exitcode or err:
          logger.log(MESSAGE, timestamp("'%s.bam' appears to be a corrupt BAM!" % bam_prefix))
          sys.exit(1)
          
     else:
          logger.log(MESSAGE, timestamp("'%s.bam' passed quickcheck" % bam_prefix))


def checkOutput(files, task):
    if isinstance(files,str):
        files = [files]
    
    dirfiles = os.listdir(".")
    for file in files:
        if os.path.basename(file) not in dirfiles:
            logger.log(MESSAGE, file)
            logger.log(MESSAGE, dirfiles)
            logger.log(MESSAGE, timestamp("%s - '%s' appears to be missing!" % (task, file)))
            sys.exit(1)


def genome_select(dirpath, target_genome):
     """
     Check genome_path and select index.
     """
     genomedirs = os.listdir(dirpath)
     if os.path.isdir(dirpath):
          if target_genome in genomedirs:
              genome_file = [dirpath + "/" + genome for genome in genomedirs if genome == target_genome][0]
              logger.log(MESSAGE, timestamp("Genome Path - '%s'" % genome_file))
              return genome_file
          else:
                logger.log(MESSAGE, timestamp("'%s' genome not available" % target_genome))
                sys.exit(1)
     else:
          logger.log(MESSAGE, timestamp("genome_path '%s' in not a vaild directory" % dirpath))
          sys.exit(1)
     

def choose_alignCommand(target_genome, nThreads, fq_files, base_name, isPaired):
    alignerPath = getFunctionPath(options.aligner)
    if isPaired == "True":
        aligners = {'bismark':"%s --genome %s --parallel %s -1 %s -2 %s" % tuple([alignerPath] + [target_genome] + [nThreads] + fq_files),
                    'walt':'%s -t %s -i %s -1 %s -2 %s -o %s.sam' % tuple([alignerPath] + [nThreads] + [target_genome] + fq_files + [base_name])}
    else:
        aligners = {'bismark':"%s --genome %s --parallel %s --se %s" % tuple([alignerPath] + [target_genome] + [nThreads] + [fq_files]),
                    'walt':'%s -t %s -i %s -r %s -o %s.sam' % tuple([alignerPath] + [nThreads] + [target_genome] + [fq_files] + [base_name])}
    return aligners
     

def aligner_select(target_genome, nThreads, fq_files, base_name, pairedReads):
    aligners = choose_alignCommand(target_genome, nThreads, fq_files, base_name, pairedReads)
    align_cmd = aligners.get(options.aligner)
    if options.non_directional == True or options.pbat == True:
        align_cmd = align_cmd + ' ' + getDirection()
    logger.log(MESSAGE,  timestamp("Align Command - '%s'" % align_cmd))
    return align_cmd

def bamMerge(bam_files, file_prefix):
     merge_cmd = '%s merge -@ %s -O BAM %s.bam %s' % (getFunctionPath("samtools"), options.threads, file_prefix, ' '.join(bam_files))
     logger.log(MESSAGE, timestamp('merging %s files into 1' % len(bam_files)))
     os.system(merge_cmd)


def bam_name_fix(base_name):
     samtoolsPath = getFunctionPath("samtools")
     if options.is_paired_end == "True":
          aligners={'bismark':'mv %s_r1_bismark_bt2_pe.bam %s.bam' % (base_name,base_name),
                        'walt':'%s view -b -@ %s -o %s.bam %s.sam' % (samtoolsPath, options.threads, base_name,base_name)}
     else:
          aligners={'bismark':'mv %s_r1_bismark_bt2.bam %s.bam' % (base_name,base_name),
                        'walt':'%s view -b -@ %s -o %s.bam %s.sam' % (samtoolsPath, options.threads, base_name,base_name)}
     fix_cmd = aligners.get(options.aligner)
     logger.log(MESSAGE,  timestamp("Correct filename - '%s'" % fix_cmd))
     return fix_cmd


def getDirection():
    if options.non_directional == True:
        aligners={'bismark':'--non_directional',
                     'walt':''}
        logger.log(MESSAGE,  timestamp('Aligning in NON_DIRECTIONAL mode'))
    elif options.pbat == True:
        aligners={'bismark':'--pbat',
                     'walt':'-P'}
        logger.log(MESSAGE,  timestamp('Aligning in PBAT mode'))
    return aligners.get(options.aligner)

#Walt code to be fixed in the future
#def align_walt(input_files, base_name, target_genome):
#     logger.log(MESSAGE, timestamp('Unzipping %s and %s for WALT' % tuple(input_files[0])))
#     cmd1 = 'gunzip --keep %s' % input_files[0][0]
#     cmd2 = 'gunzip --keep %s' % input_files[0][1]
#     logger.log(MESSAGE, timestamp('Unzip r1: %s' % cmd1))
#     os.system(cmd1)
#     logger.log(MESSAGE, timestamp('Unzip r2: %s' % cmd2))
#     os.system(cmd2)
#     fq_files = []
#     for file in os.listdir(os.getcwd()):
#          if file.endswith('fq'):
#                fq_files.append(os.getcwd() + '/' + file)
#     fq_files.sort()
#     walt_cmd = aligner_select(target_genome, options.threads, [fq_files], base_name)
#     return walt_cmd

#if data directory contains SRA files, convert these to fastqs
def prepFastQ(input_files,logger, logger_mutex):
    dataFiles=[file for file in os.listdir(options.data_dir)]
    for sra in input_files:
        if Path(sra).stem + "_1.fastq.gz" not in dataFiles and Path(sra).stem + "_2.fastq.gz" not in dataFiles:
            cmd = ("%s %s --split-files --threads %s --gzip --outdir %s --sra-id %s --tmpdir ./" % (getFunctionPath("python"),pipeline_path + '/parallel-fastq-dump.py',
                   options.threads,options.data_dir,sra))
            logger.log(MESSAGE,  timestamp("Running Command: %s" % cmd))
            exitcode, out, err = execute_cmd(cmd)

            with logger_mutex:
                logger.log(MESSAGE,"Extracted fastq from sra")

sraFiles = []
for root, dirs, files in os.walk(options.data_dir):
    for file in files:
        #append the file name to the list
        if file.endswith('sra'):
            sraFiles.append(os.path.join(root,file))     

if sraFiles:
    prepFastQ(sraFiles, logger, logger_mutex)

#locate fastq files in data directory to run through pipeline
def detect_input_type(input_files):
    if all([re.match(".+\\.f(ast)?q(.gz)?", file) for file in input_files]) and input_files:
        logger.log(MESSAGE, timestamp("Running from fastq"))
    elif any([re.match(".+\\.f(ast)?q(.gz)?", file) for file in input_files]) and input_files:
        logger.log(MESSAGE, timestamp("Data directory contains a mixture of file types, including fastqs. Running from fastq."))
        infiles = [file for file in input_files if re.match(".+\\.f(ast)?q(.gz)?", file)]
    else:
        logger.log(MESSAGE, timestamp("No fastq files exist in data directory"))
        sys.exit(1)    
    
infiles=[os.path.join(options.data_dir, file) for file in os.listdir(options.data_dir) if os.path.isfile(os.path.join(options.data_dir, file))]
infiles.sort()
detect_input_type(infiles)

#pair fastq files and prepare for ruffus 
if options.is_paired_end == "True":
    infiles = [[r1,r2] for r1,r2 in zip(infiles[::2], infiles[1::2])]

#break if sample_name is not within convention
tmp_pd = pd.read_csv(pipeline_path + "/data/TissueToEmbryoMap.csv")
test = options.sample_name.split("_")
if len(test) >= 4 and test[0].lower() in tmp_pd.tissue.str.lower().tolist():
    logger.log(MESSAGE,timestamp("Sample name is good, proceeding to pipeline"))
else:
    logger.log(MESSAGE,timestamp("Check sample name before proceeding"))


@collate(infiles, formatter("(.*/)*(?P<file_details>.+_[rR]*)[12](?P<set_number>(_...)*)\.f(ast)*q(.gz)*"), 
                                     ['{file_details[0]}1{set_number[0]}_val_1.fq.gz',
                                      '{file_details[0]}2{set_number[0]}_val_2.fq.gz'], logger, logger_mutex)

def trim_fastq(input_files, output_paired_files, logger, logger_mutex):
    trim_profile = str(options.trim_profile).split(",")
    L_trim_lengths = []
    L_trim_options = ["--clip_R1 ", "--clip_R2 ", "--three_prime_clip_R1 ", "--three_prime_clip_R2 "]
    
    if trim_profile[0] == "em-seq":
        trim_profile[0] = "8"
    elif trim_profile[0] == "bs-seq":
        trim_profile[0] = "10"
    elif trim_profile[0] == "no-trim":
        trim_profile[0] = "0"
    else:
        try:
            int(trim_profile[0])
        except:
            raise Exception("Invalid --trim_profile. Must be a single integer, a comma seperated list of 4 integers, or one of the following: em-seq, bs-seq, no-trim")
    
    for x in range(4):
        trim_profile.append(trim_profile[0])
        L_trim_lengths.append(L_trim_options[x] + trim_profile[x])

    trimPath = getFunctionPath("trim_galore")
    if options.is_paired_end == "True":
        if all([len(files) == 2 for files in input_files]):
            threads = int(options.threads) // 4 if int(options.threads) > 4 else 1
          
            cmd=('%s --fastqc --fastqc_args "--noextract" --gzip --cores %s %s %s %s %s --paired %s %s' % tuple([trimPath] + [str(threads)] + L_trim_lengths + input_files[0]))
            exitcode, out, err = execute_cmd(cmd)
        else:
            raise Exception("Unpaired files in input.")
        
    else:
        raise Exception("Fastq processing only supports paired end sequencing")
    
    checkOutput(output_paired_files, "trim_fastq")
     
    with logger_mutex:
        logger.log(MESSAGE,  timestamp("trim_galore done"))

      
@merge([trim_fastq], [options.sample_name + "_r1.fq.gz", options.sample_name + "_r2.fq.gz"] if options.is_paired_end == "True" else options.sample_name + "_r1.fq.gz", logger, logger_mutex)

def merge_fastq(input_files, output_files, logger, logger_mutex):
    
    logger.log(MESSAGE, timestamp("Merging files before mapping"))
    if options.is_paired_end == "True":
        r1 = " " .join([files[0] for files in input_files])
        r2 = " " .join([files[1] for files in input_files])
        logger.log(MESSAGE, timestamp("Merging Read 1 files: %s" % r1))
        cmd_r1 = "cat %s > %s" % (r1, output_files[0])
        os.system(cmd_r1)
        logger.log(MESSAGE, timestamp("Merging Read 2 files: %s" % r2))
        cmd_r2 = "cat %s > %s" % (r2, output_files[1])
        os.system(cmd_r2)
    else:
        r1 = " " .join([files[0] for files in input_files])
        cmd_r1 = "cat %s > %s" % (r1, output_files)
        os.system(cmd_r1)
          
    checkOutput(output_files, "trim_fastq") 
    with logger_mutex:
        logger.log(MESSAGE,  timestamp("merge_fastq done"))
     
  
@transform(merge_fastq, regex(r"_r[12].fq.gz"), ".bam", logger, logger_mutex)

def align_fastq(input_files, output_file, logger, logger_mutex):
    genome_file=genome_select(options.genome_path, options.genome)
    base_name=output_file[:-4]
    if options.aligner == 'walt':
        cmd = align_walt(input_files, base_name, genome_file)
    else:
        #bismark uses ~5 time the number of cores as specified in --parallel
        threads = int(options.threads) // 5 if int(options.threads) > 5 else 1
        cmd=aligner_select(genome_file, str(threads), input_files, 
                           base_name, options.is_paired_end)
    logger.log(MESSAGE, timestamp(cmd))
    os.system(cmd)
    #name fix
    mv_cmd=bam_name_fix(base_name)
    os.system(mv_cmd)
    checkBam(options.sample_name)
    
    with logger_mutex:
        logger.log(MESSAGE,  timestamp("Alignment done"))

@transform(align_fastq, regex(r".bam"), ["_s.bam", "_s.bam.bai"], logger, logger_mutex)

def sort_bam(input_file, output_file, logger, logger_mutex):
     samtoolsPath = getFunctionPath("samtools")
     os.system("%s sort -m 2G -O bam -@ 12 -o %s %s " % (samtoolsPath, output_file[0], input_file))
     checkBam(options.sample_name + "_s")
     index_cmd = ("%s index %s" % (samtoolsPath, output_file[0]))
     os.system(index_cmd)
     checkOutput(output_file, "sort_bam")
     
     with logger_mutex:
          logger.log(MESSAGE,  timestamp("Sam sorted and converted to bam"))

#@transform(compress_sam, regex(r".bam"), "_d.bam")
@transform(sort_bam, formatter(".*_s.bam$"), ['{basename[0]}d.bam','{basename[0]}d_flagstat.txt'], logger, logger_mutex)
def mDuplicates(input_file, output_file, logger, logger_mutex):
    if not os.path.exists("./tmp"):
        os.mkdir("./tmp")
    picardPath = getFunctionPath("picard")
    samtoolsPath = getFunctionPath("samtools")
    cmd="%s -Xms8g MarkDuplicates I=%s O=%s M=%s_picard_MarkDuplicates_metrics.test ASSUME_SORT_ORDER=coordinate REMOVE_DUPLICATES=true TAGGING_POLICY=All CREATE_INDEX=true TMP_DIR=./tmp" % (picardPath, input_file[0], output_file[0], options.sample_name)
    logger.log(MESSAGE,  timestamp(cmd))
    os.system(cmd)
    checkBam(options.sample_name + "_sd")
    flagstat_cmd = '%s flagstat %s > %s' % (samtoolsPath, output_file[0], output_file[1])
    os.system(flagstat_cmd)
    checkOutput(output_file, "mDuplicates")
    
    with logger_mutex:
        logger.log(MESSAGE,  timestamp("Picard Completed"))

@transform(mDuplicates, formatter('.*_sd.bam$'), ['{basename[0]}_genomeCoverageBed.txt','{basename[0]}_coverage.txt'], logger, logger_mutex)
def calculateCoverage(input_file, output_file, logger, logger_mutex):
     samtoolsPath = getFunctionPath("samtools")
     covBedpath = getFunctionPath("genomeCoverageBed")
     cmd = "%s view -b -F 0x400 %s | %s -ibam - -g %s/%s/%s.genome > %s" % (samtoolsPath, input_file[0], covBedpath, options.genome_path, options.genome, options.genome, output_file[0])     
     os.system(cmd)
     covFile = pd.read_table(output_file[0], header=None, names=['chr','depth','base_count','chr_size_bp','fraction'])
     genomeCov = covFile[covFile['chr'] == 'genome']
     averageCov = sum(genomeCov['depth'] * genomeCov['base_count'])/genomeCov.loc[genomeCov.index[0],'chr_size_bp']
     pd.DataFrame(data={'sample_name':options.sample_name,'Genome':options.genome, 'Average_Coverage':averageCov}, index=['coveragedetails']).to_csv(path_or_buf=output_file[1], sep = '\t')
     checkOutput(output_file, "calculateCoverage")
     
     with logger_mutex:
          logger.log(MESSAGE, timestamp('Average genomic coverage = %sx' % averageCov))

@transform(mDuplicates, regex(r".bam$"), ["_CpG.bedGraph",'_OB.svg','_OT.svg'], logger, logger_mutex)
def call_meth(input_file, output_file, logger, logger_mutex):
     methPath = getFunctionPath("MethylDackel")
     bias_cmd=("%s mbias -@ %s %s/%s/%s.fa %s %s" % (methPath, options.threads, options.genome_path, options.genome, options.genome, input_file[0], re.sub(pattern = '\.bam$', repl='',string = input_file[0])))
     os.system(bias_cmd)
     cmd=("%s extract -@ %s --mergeContext %s%s/%s.fa %s" % (methPath, options.threads, options.genome_path, options.genome, options.genome, input_file[0]))
     os.system(cmd)
     checkOutput(output_file, "call_meth")
     
     with logger_mutex:
          logger.log(MESSAGE,  timestamp("MethylDackel Completed"))

@transform(call_meth, regex(r"_sd_CpG.bedGraph"), ["_PMD.bed", "_UMRLMR.bed", "_wPMD_UMRLMR.bed", "_sd_CpG.tdf"], logger, logger_mutex)
def methylseekrAndTDF(input_file, output_file, logger, logger_mutex):
     threads = int(options.threads) // 8 if int(options.threads) > 8 else 1 
     Rscript_cmd = "Rscript %s/Rscripts/CallMethylseekrRegions_and_convertMethCallsToTdf.R -g %s -i %s -t %s/data/TissueToEmbryoMap.csv -e %s -p %s" % (pipeline_path, options.genome, input_file[0], pipeline_path, options.genome_path, str(threads))
     print(Rscript_cmd)
     os.system(Rscript_cmd)
     checkOutput(output_file, "methylseekrAndTDF")
     
     with logger_mutex:
          logger.log(MESSAGE,  timestamp("MethylSeekR and toTDF Completed"))

cmdline.run(options)
