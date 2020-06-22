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

#define all paths to sofware (Popen will not use $PATH)
def getFunctionPath(target):
     if re.sub(".*\\.","",target) == "jar":
          proc = Popen(shlex.split('which %s' % target), stdout = PIPE, stderr = PIPE)
          functionPath = re.sub("\\n", "",proc.communicate()[0].decode("utf-8"))
     else:
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
parser.add_argument("--aligner_threads", help="Speed up alignment by increasing number of threads. Values depends on the aligner. Defaults to 5",
                          default=5)
parser.add_argument("--pbat", action='store_true',
                          help="Specify when aligning pbat library")
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
# Utility functions
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
    if options.pbat == 'True':
        align_cmd = align_cmd + ' ' + getPBAT()
        logger.log(MESSAGE,  timestamp('Aligning in PBAT mode'))
        logger.log(MESSAGE,  timestamp("Align Command - '%s'" % align_cmd))
    return align_cmd


def bamMerge(bam_files, file_prefix):
     thread_factor = 1
     if options.aligner == 'bismark':
          thread_factor = 5
     merge_cmd = '%s merge -@ %s -O BAM %s.bam %s' % (getFunctionPath("samtools"), options.aligner_threads*thread_factor, file_prefix, ' '.join(bam_files))
     logger.log(MESSAGE, timestamp('merging %s files into 1' % len(bam_files)))
     os.system(merge_cmd)


def bam_name_fix(base_name):
     samtoolsPath = getFunctionPath("samtools")
     if options.is_paired_end == "True":
          aligners={'bismark':'mv %s_r1_bismark_bt2_pe.bam %s.bam' % (base_name,base_name),
                        'walt':'%s view -b -@ %s -o %s.bam %s.sam' % (samtoolsPath, options.aligner_threads, base_name,base_name)}
     else:
          aligners={'bismark':'mv %s_r1_bismark_bt2.bam %s.bam' % (base_name,base_name),
                        'walt':'%s view -b -@ %s -o %s.bam %s.sam' % (samtoolsPath, options.aligner_threads, base_name,base_name)}
     fix_cmd = aligners.get(options.aligner)
     logger.log(MESSAGE,  timestamp("Correct filename - '%s'" % fix_cmd))
     return fix_cmd


def getPBAT():
     aligners={'bismark':'--pbat',
                  'walt':'-P'}
     return aligners.get(options.aligner)


def detect_input_type(input_files):
     if all([re.match(".+\\.sra", file) for file in input_files]):
          logger.log(MESSAGE, timestamp("Running form sra"))
          return True, False
     elif all([re.match(".+\\.f(ast)?q(.gz)?", file) for file in input_files]):
          logger.log(MESSAGE, timestamp("Running from fastq"))
          return False, True
     else:
          logger.log(MESSAGE, timestamp("Input files appear to be in incorrect or mixed formats"))
          sys.exit(1)
     


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
#     walt_cmd = aligner_select(target_genome, options.aligner_threads, [fq_files], base_name)
#     return walt_cmd


def getFastqPairs():
     target_files = []
     for file in os.listdir(options.data_dir):
          if file.endswith(".fastq.gz"):
                target_files.append(options.data_dir + file)
     target_files.sort()
     r1 = (target_files[::2])
     r2 = (target_files[1::2])
     return [(r1, options.sample_name + '_r1.fastq.gz'),
                (r2, options.sample_name + '_r2.fastq.gz')]


# Locate the files
sraFiles = []
for file in os.listdir(options.data_dir):
     if file.endswith('sra'):
          sraFiles.append(options.data_dir + file)
          
          
infiles=[os.path.join(options.data_dir, file) for file in os.listdir(options.data_dir)]
infiles.sort()
run_if_sra, run_if_fastq = detect_input_type(infiles)
#pair files and prepare for ruffus 
if options.is_paired_end == "True":
    infiles = [[r1,r2] for r1,r2 in zip(infiles[::2], infiles[1::2])]

#break if sample_name is not within convention
tmp_pd = pd.read_csv(pipeline_path + "/data/TissueToEmbryoMap.csv")
test = options.sample_name.split("_")
if len(test) >= 4 and test[0].lower() in tmp_pd.tissue.str.lower().tolist():
    logger.log(MESSAGE,timestamp("Sample name is good, proceeding to pipeline"))
else:
    logger.log(MESSAGE,timestamp("Check sample name before proceeding"))


#start the pipeline
@active_if(run_if_sra)
@transform(sraFiles, formatter(".*.sra$"), [os.getcwd() + "/{basename[0]}_1.fastq.gz", os.getcwd() + "/{basename[0]}_2.fastq.gz"] if options.is_paired_end == "True" else [os.getcwd() + "/{basename[0]}_1.fastq.gz"], logger, logger_mutex)
#extract from sra
def sraToFastq(input_file, output_files, logger, logger_mutex):
     cmd = ("%s --split-files --gzip %s" % (getFunctionPath("fastq-dump"),input_file))
     logger.log(MESSAGE,  timestamp("Running Command: %s" % cmd))
     exitcode, out, err = execute_cmd(cmd)
     checkOutput(output_files, "sraToFastq")
     
     with logger_mutex:
          logger.log(MESSAGE,"Extracted fastq from sra")


@active_if(run_if_sra)         
@collate(sraToFastq, regex("_[12].fastq.gz"), ["_1_val_1.fq.gz","_2_val_2.fq.gz"] if options.is_paired_end == "True" else ["_1_trimmed.fq.gz"], logger, logger_mutex)

def trim_fastq(input_files, output_paired_files, logger, logger_mutex):
    trimPath = getFunctionPath("trim_galore")
    if options.is_paired_end == "True":
        if all([len(files) == 2 for files in input_files]):
            cmd=('%s --fastqc --fastqc_args "--noextract" --gzip --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --paired %s %s' % tuple([trimPath] + input_files[0]))
            exitcode, out, err = execute_cmd(cmd)
        else:
            raise Exception("Unpaired files present. Check output from sraToFastq")
        
    else:
        cmd=('%s --fastqc --fastqc_args "--noextract" --gzip --clip_R1 10 --three_prime_clip_R1 10 %s' % tuple([trimPath] + input_files[0]))
        exitcode, out, err = execute_cmd(cmd)
    
    checkOutput(output_paired_files, "trim_fastq")
     
    with logger_mutex:
        logger.log(MESSAGE,  timestamp("trim_galore done"))


@active_if(run_if_fastq)         
@collate(infiles, formatter("(.*/)*(?P<file_details>.+_[rR]*)[12](?P<set_number>(_...)*)\.f(ast)*q(.gz)*"), 
                                     ['{file_details[0]}1{set_number[0]}_val_1.fq.gz',
                                      '{file_details[0]}2{set_number[0]}_val_2.fq.gz'], logger, logger_mutex)
def trim_fastq2(input_files, output_paired_files, logger, logger_mutex):
    trimPath = getFunctionPath("trim_galore")
    if options.is_paired_end == "True":
        if all([len(files) == 2 for files in input_files]):
            cmd=('%s --fastqc --fastqc_args "--noextract" --gzip --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --paired %s %s' % tuple([trimPath] + input_files[0]))
            exitcode, out, err = execute_cmd(cmd)
        else:
            raise Exception("Unpaired files in input.")
        
    else:
        raise Exception("Fastq processing only supports paired end sequencing")
    
    checkOutput(output_paired_files, "trim_fastq")
     
    with logger_mutex:
        logger.log(MESSAGE,  timestamp("trim_galore done"))

      
@merge([trim_fastq, trim_fastq2], [options.sample_name + "_r1.fq.gz", options.sample_name + "_r2.fq.gz"] if options.is_paired_end == "True" else options.sample_name + "_r1.fq.gz", logger, logger_mutex)

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
        cmd=aligner_select(genome_file, options.aligner_threads, input_files, 
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
    javaPath = getFunctionPath("java")
    picardPath = getFunctionPath("picard.jar")
    samtoolsPath = getFunctionPath("samtools")
    cmd="%s -Xms8g -jar %s MarkDuplicates I=%s O=%s M=%s_picard_MarkDuplicates_metrics.test ASSUME_SORT_ORDER=coordinate REMOVE_DUPLICATES=true TAGGING_POLICY=All CREATE_INDEX=true TMP_DIR=./tmp" % (javaPath, picardPath, input_file[0], output_file[0], options.sample_name)
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
     cmd = "%s view -b -F 0x400 %s | %s -ibam - -g %s%s/%s.genome > %s" % (samtoolsPath, input_file[0], covBedpath, options.genome_path, options.genome, options.genome, output_file[0])
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
     bias_cmd=("%s mbias -@ %s %s%s/%s.fa %s %s" % (methPath, options.aligner_threads, options.genome_path, options.genome, options.genome, input_file[0], re.sub(pattern = '\.bam$', repl='',string = input_file[0])))
     os.system(bias_cmd)
     cmd=("%s extract -@ %s --mergeContext %s%s/%s.fa %s" % (methPath, options.aligner_threads, options.genome_path, options.genome, options.genome, input_file[0]))
     os.system(cmd)
     checkOutput(output_file, "call_meth")
     
     with logger_mutex:
          logger.log(MESSAGE,  timestamp("MethylDackel Completed"))

@transform(call_meth, regex(r"_sd_CpG.bedGraph"), ["_PMD.bed", "_UMRLMR.bed", "_wPMD_UMRLMR.bed", "_sd_CpG.tdf"], logger, logger_mutex)
def methylseekrAndTDF(input_file, output_file, logger, logger_mutex):
     Rscript_cmd = "Rscript %s/Rscripts/CallMethylseekrRegions_and_convertMethCallsToTdf.R -g %s -i %s -t %s/data/TissueToEmbryoMap.csv -e %s -p %s" % (pipeline_path, options.genome, input_file[0], pipeline_path, options.genome_path, options.aligner_threads)
     os.system(Rscript_cmd)
     checkOutput(output_file, "methylseekrAndTDF")
     
     with logger_mutex:
          logger.log(MESSAGE,  timestamp("MethylSeekR and toTDF Completed"))

cmdline.run(options)