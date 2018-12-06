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

#define all paths to sofware (Popen will not use $PATH)
def getFunctionPath(target):
    if re.sub(".*\\.","",target) == "jar":
        proc = Popen(shlex.split('locate %s' % target), stdout = PIPE, stderr = PIPE)
        functionPath = re.sub("\\n", "",proc.communicate()[0].decode("utf-8"))
    else:
        proc = Popen(shlex.split('which %s' % target), stdout = PIPE, stderr = PIPE)
        functionPath = re.sub("\\n", "",proc.communicate()[0].decode("utf-8"))
    return functionPath


#softwarePaths = {'bismark': getFunctionPath('bismark'),
#                 'samtools': getFunctionPath('samtools'),
#                 'java': getFunctionPath('java'),
#                 'picard': getFunctionPath('picard.jar'),
#                 'trim_galore': getFunctionPath('trim_galore'),
#                 'fastq-dump': getFunctionPath('fastq-dump')}


# get important directory locations
pipeline_path = os.path.dirname(os.path.abspath(__file__))
genomes_path = pipeline_path + '/Genomes'

# Parse command line arguments
parser = cmdline.get_argparse(description='test.py - Automated WGBS processing pipeline')

# Add command line arguments
#removing multiple alingers for the time being
#parser.add_argument("--aligner", help="Select prefered aligner [bismark, bwameth]. Defaults to bismark",
#                    default="bismark")
parser.add_argument("--genome", help="Genome reads are aligned too. Check genome files for your prefered aligner are in " + genomes_path + ". Defaults to hg38 (hg38)",
                    default="hg38")
parser.add_argument("--data_dir", help="Directory for fastq files")
parser.add_argument("--sampleID", help="Sample name. Used for directory name and output file names")
parser.add_argument("--aligner_threads", help="Speed up alignment by increasing number of threads. Values depends on the aligner. Defaults to 5",
                    default=5)
parser.add_argument("--pbat", help="True if library method uses post-bis adapter tagging",
                    default=False)
parser.add_argument("--file_type", help="Starting file type (sra or fastq)")
parser.add_argument("--isPairedEnd",help = "IS the libarary paired end (defaults to True)",
                    default="True")
parser.add_argument("--bowenPath",help = "path to bowen volume (defaults to HPC path '/OSM/CBR/HB_FSP_TBI/work/'",
                    default="/OSM/CBR/HB_FSP_TBI/work/")
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


def check_valid_bam(file):
    logger.log(MESSAGE,  timestamp("Checking '%s' is a valid BAM" % file))
    exitcode, out, err = execute_cmd("samtools view -H %s" % file)
    if exitcode or err:
        raise JobFailException("'%s' appears to be a corrupt BAM!" % file)


def check_expected_output(files, taskname):
    for f in files:
        if not os.path.isfile(f):
            raise JobFailException('%s: Expected output %s is not present' % (taskname, f))


def genome_select():
    genomes={'hg19':'%spipeline_data/Genomes/hg19/' % options.bowenPath, 'hg38':'%spipeline_data/Genomes/hg38/' % options.bowenPath}
    if options.aligner == 'bismark':
        genome_file = genomes.get(options.genome)
    if options.aligner == 'walt':
        genome_file = "%s%s.dbindex" % (genomes.get(options.genome),options.genome)
    logger.log(MESSAGE, timestamp("Genome Path - '%s'" % genome_file))
    return genome_file


def choose_alignCommand(target_genome, nThreads, fq_files, base_name, isPaired):
    alignerPath = getFunctionPath(options.aligner)
    if isPaired == "True":
        aligners = {'bismark':"%s --genome %s --parallel %s -1 %s -2 %s" % tuple([alignerPath] + [target_genome] + [nThreads] + list(fq_files[0])),
                    'walt':'%s -t %s -i %s -1 %s -2 %s -o %s.sam' % tuple([alignerPath] + [nThreads] + [target_genome] + list(fq_files[0]) + [base_name])}
    else:
        aligners = {'bismark':"%s --genome %s --parallel %s --se %s" % tuple([alignerPath] + [target_genome] + [nThreads] + list(fq_files[0])),
                    'walt':'%s -t %s -i %s -r %s -o %s.sam' % tuple([alignerPath] + [nThreads] + [target_genome] + list(fq_files[0]) + [base_name])}
    return aligners
    

def aligner_select(target_genome, nThreads, fq_files, base_name, pairedReads):
    aligners= choose_alignCommand(target_genome, nThreads, fq_files, base_name, pairedReads)
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
    if options.isPairedEnd == "True":
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


#Walt code to be fixed in the future
#def align_walt(input_files, base_name, target_genome):
#    logger.log(MESSAGE, timestamp('Unzipping %s and %s for WALT' % tuple(input_files[0])))
#    cmd1 = 'gunzip --keep %s' % input_files[0][0]
#    cmd2 = 'gunzip --keep %s' % input_files[0][1]
#    logger.log(MESSAGE, timestamp('Unzip r1: %s' % cmd1))
#    os.system(cmd1)
#    logger.log(MESSAGE, timestamp('Unzip r2: %s' % cmd2))
#    os.system(cmd2)
#    fq_files = []
#    for file in os.listdir(os.getcwd()):
#        if file.endswith('fq'):
#            fq_files.append(os.getcwd() + '/' + file)
#    fq_files.sort()
#    walt_cmd = aligner_select(target_genome, options.aligner_threads, [fq_files], base_name)
#    return walt_cmd


def getFastqPairs():
    target_files = []
    for file in os.listdir(options.data_dir):
        if file.endswith(".fastq.gz"):
            target_files.append(options.data_dir + file)
    target_files.sort()
    r1 = (target_files[::2])
    r2 = (target_files[1::2])
    return [(r1, options.sampleID + '_r1.fastq.gz'),
            (r2, options.sampleID + '_r2.fastq.gz')]


# The actual pipeline
sraFiles = []
for file in os.listdir(options.data_dir):
    if file.endswith('sra'):
        sraFiles.append(options.data_dir + file)
        
        
infiles=[]
for file in os.listdir(options.data_dir):
    target_suffix = options.file_type
    if target_suffix == 'fastq':
        target_suffix = ('fastq.gz', 'fastq', 'fq.gz')
    if file.endswith(target_suffix):
        infiles.append(options.data_dir + file)
infiles.sort()
inType = options.file_type
run_if_sra = inType == 'sra'
run_if_fastq = inType == 'fastq'

def type_select():
    if run_if_sra == True:
        return('Running from sra')
    else:
        return('Running from fastq')

    
type_text = type_select() + ': ' + str(infiles)
logger.log(MESSAGE, timestamp(type_text))


@active_if(run_if_sra)
@transform(sraFiles, formatter(".*.sra$"), [os.getcwd() + "/{basename[0]}_1.fastq.gz", os.getcwd() + "/{basename[0]}_2.fastq.gz"] if options.isPairedEnd == "True" else [os.getcwd() + "/{basename[0]}_1.fastq.gz"], logger, logger_mutex)
#extract from sra
def sraToFastq(input_file, output_files, logger, logger_mutex):
    cmd = ("%s --split-files --gzip %s" % (getFunctionPath("fastq-dump"),input_file))
    logger.log(MESSAGE,  timestamp("Running Command: %s" % cmd))
    exitcode, out, err = execute_cmd(cmd)
    
    with logger_mutex:
        logger.log(MESSAGE,"Extracted fastq from sra")


@active_if(run_if_sra)       
@collate(sraToFastq, regex("_[12].fastq.gz"), ["_1_val_1.fq.gz","_2_val_2.fq.gz"] if options.isPairedEnd == "True" else ["_1_trimmed.fq.gz"], logger, logger_mutex)

def trim_fastq(input_files, output_paired_files, logger, logger_mutex):
    trimPath = getFunctionPath("trim_galore")
    if options.isPairedEnd == "True":
        if len(input_files[0]) < 2:
            raise Exception("One of read pairs %s missing" % (input_files[0]))
        if len(input_files[0]) == 2:
            cmd=('%s --fastqc --fastqc_args "--extract" --gzip --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --paired %s %s' % tuple([trimPath] + input_files[0]))
            exitcode, out, err = execute_cmd(cmd)
        if len(input_files[0]) > 2 and (len(input_files[0])/2).is_integer() == False:
            raise Exception("File missing in %s" % (input_files))
        if len(input_files[0]) > 2 and (len(input_files[0])/2).is_integer() == True:
            input_files.sort()
            target_files_1 = ' '.join(input_files[0][::2])
            target_files_2 = ' '.join(input_files[0][1::2])
            logger.log(MESSAGE, timestamp(target_files_1))
            logger.log(MESSAGE, timestamp(target_files_2))
            for r1, r2 in zip(target_files_1, target_files_2):
                cmd=('%s --fastqc --fastqc_args "--extract" --gzip --paired %s %s' % tuple(trimPath,r1, r2))
                exitcode, out, err = execute_cmd(cmd)
    else:
        cmd=('%s --fastqc --fastqc_args "--extract" --gzip %s' % tuple([trimPath] + input_files[0]))
        exitcode, out, err = execute_cmd(cmd)

    with logger_mutex:
        logger.log(MESSAGE,  timestamp("trim_galore done"))


@active_if(run_if_fastq)       
@collate(infiles, formatter("(.*/)*(?P<file_details>.+_[rR]*)[12](?P<set_number>(_...)*)\.f(ast)*q(.gz)*"), 
                            ['{file_details[0]}1{set_number[0]}_val_1.fq.gz',
                             '{file_details[0]}2{set_number[0]}_val_2.fq.gz'], logger, logger_mutex)
def trim_fastq2(input_files, output_paired_files, logger, logger_mutex):
    trimPath = getFunctionPath("trim_galore")
    if len(input_files) < 2:
        raise Exception("One of read pairs %s missing" % (input_files))
    if len(input_files) == 2:
        cmd_args = tuple([trimPath]) + input_files
        cmd=('%s --fastqc --fastqc_args "--extract" --gzip --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --paired %s %s' % cmd_args)
        exitcode, out, err = execute_cmd(cmd)
    if len(input_files) > 2 and (len(input_files)/2).is_integer() == False:
        raise Exception("File missing in %s" % (input_files))
    if len(input_files) > 2 and (len(input_files)/2).is_integer() == True:
        input_files.sort()
        target_files_1 = ' '.join(input_files[::2])
        target_files_2 = ' '.join(input_files[1::2])
        for r1, r2 in zip(target_files_1, target_files_2):
            cmd=('%s --fastqc --fastqc_args "--extract" --gzip --paired %s %s' % tuple(trimPath, r1, r2))
            exitcode, out, err = execute_cmd(cmd)

    with logger_mutex:
        logger.log(MESSAGE,  timestamp("trim_galore done"))

     
@collate([trim_fastq, trim_fastq2], formatter("R*[12].*_(?:val_[12]|trimmed).fq.gz$"), ["{path[0]}/" + options.sampleID + "_r1.fq.gz", "{path[0]}/" + options.sampleID + "_r2.fq.gz"] if options.isPairedEnd == "True" else ["{path[0]}/" + options.sampleID + "_r1.fq.gz"], logger, logger_mutex)

def merge_fastq(input_files, output_files, logger, logger_mutex):
    if len(input_files) == 1:
        input_flat = [item for sublist in input_files for item in sublist]
        logger.log(MESSAGE, timestamp('Only 1 fastq pair, no need to merge. Renaming to keep pipline flowing'))
        os.system('mv %s %s' % (input_flat[0], output_files[0]))
        os.system('mv %s %s' % (input_flat[1], output_files[1]))
    elif options.isPairedEnd == "True":
        r1 = []
        r2 = []
        for targets in input_files:
            r1.append(targets[0])
            r2.append(targets[1])
        
        logger.log(MESSAGE, timestamp('merging %s into %s' % (' '.join(r1), output_files[0])))
        cmd = 'cat %s > %s' % (' '.join(r1), output_files[0])
        os.system(cmd)
        logger.log(MESSAGE, timestamp('merging %s into %s' % (' '.join(r2), output_files[1])))
        cmd_next = 'cat %s > %s' % (' '.join(r2), output_files[1])
        os.system(cmd_next)
    else:
        input_flat = [item for sublist in input_files for item in sublist]
        logger.log(MESSAGE, timestamp('merging %s into %s' % (' '.join(input_flat), output_files[0])))
        cmd = 'cat %s > %s' % (' '.join(input_flat), output_files[0])
        os.system(cmd)
        
    with logger_mutex:
         logger.log(MESSAGE,  timestamp("merge_fastq done"))
    
    

@collate(merge_fastq, formatter('.*r[12].fq.gz'), "{path[0]}/" + options.sampleID + ".bam", logger, logger_mutex)
def align_fastq(input_files, output_file, logger, logger_mutex):
    genome_file=genome_select()
    base_name=output_file[:-4]
    if options.aligner == 'walt':
        cmd = align_walt(input_files, base_name, genome_file)
    else:
        cmd=aligner_select(genome_file, options.aligner_threads, input_files, base_name, options.isPairedEnd)
    logger.log(MESSAGE, timestamp(cmd))
    os.system(cmd)
    #name fix
    mv_cmd=bam_name_fix(base_name)
    os.system(mv_cmd)
    
    with logger_mutex:
        logger.log(MESSAGE,  timestamp("Alignment done"))

@transform(align_fastq, regex(r".bam"), ["_s.bam", "_s.bam.bai"], logger, logger_mutex)

def sort_bam(input_file, output_file, logger, logger_mutex):
    samtoolsPath = getFunctionPath("samtools")
    os.system("%s sort -m 2G -O bam -@ 12 -o %s %s " % (samtoolsPath, output_file[0], input_file))
    index_cmd = ("%s index %s" % (samtoolsPath, output_file[0]))
    os.system(index_cmd)
    with logger_mutex:
        logger.log(MESSAGE,  timestamp("Sam sorted and converted to bam"))

#@transform(compress_sam, regex(r".bam"), "_d.bam")
@transform(sort_bam, formatter(".*_s.bam$"), ['{basename[0]}d.bam','{basename[0]}d_flagstat.txt'], logger, logger_mutex)
def mDuplicates(input_file, output_file, logger, logger_mutex):
    os.mkdir("./tmp")
    javaPath = getFunctionPath("java")
    picardPath = getFunctionPath("picard.jar")
    cmd="%s -Xms8g -jar %s MarkDuplicates I=%s O=%s M=%s_picard_MarkDuplicates_metrics.test ASSUME_SORT_ORDER=coordinate REMOVE_DUPLICATES=false TAGGING_POLICY=All CREATE_INDEX=true TMP_DIR=./tmp" % (javaPath, picardPath, input_file[0], output_file[0], options.sampleID)
    logger.log(MESSAGE,  timestamp(cmd))
    os.system(cmd)
    flagstat_cmd = '/home/loc100/miniconda3/bin/samtools flagstat %s > %s' % (output_file[0], output_file[1])
    os.system(flagstat_cmd)
    with logger_mutex:
        logger.log(MESSAGE,  timestamp("Picard Complteted"))

@transform(mDuplicates, formatter('.*_sd.bam$'), ['{basename[0]}_genomeCoverageBed.txt','{basename[0]}_coverage.txt'], logger, logger_mutex)
def calculateCoverage(input_file, output_file, logger, logger_mutex):
    samtoolsPath = getFunctionPath("samtools")
    covBedpath = getFunctionPath("genomeCoverageBed")
    cmd = "%s view -b -F 0x400 %s | %s -ibam - -g /media/bowen_work/pipeline_data/Genomes/%s/%s.genome > %s" % (samtoolsPath, input_file[0], covBedpath, options.genome, options.genome, output_file[0])
    os.system(cmd)
    covFile = pd.read_table(output_file[0], header=None, names=['chr','depth','base_count','chr_size_bp','fraction'])
    genomeCov = covFile[covFile['chr'] == 'genome']
    averageCov = sum(genomeCov['depth'] * genomeCov['base_count'])/genomeCov.ix[genomeCov.index[0],'chr_size_bp']
    pd.DataFrame(data={'SampleID':options.sampleID,'Genome':options.genome, 'Average_Coverage':averageCov}, index=['coveragedetails']).to_csv(path_or_buf=output_file[1], sep = '\t')
    with logger_mutex:
        logger.log(MESSAGE, timestamp('Average genomic coverage = %sx' % averageCov))

@transform(mDuplicates, regex(r".bam$"), ["_CpG.bedGraph",'_OB.svg','_OT.svg'], logger, logger_mutex)
def call_meth(input_file, output_file, logger, logger_mutex):
    methPath = getFunctionPath("MethylDackel")
    bias_cmd=("%s mbias -@ %s /media/bowen_work/pipeline_data/Genomes/%s/%s.fa %s %s" % (methPath, options.aligner_threads, options.genome, options.genome, input_file[0], re.sub(pattern = '\.bam$', repl='',string = input_file[0])))
    os.system(bias_cmd)
    cmd=("%s extract -@ %s --mergeContext /media/bowen_work/pipeline_data/Genomes/%s/%s.fa %s" % (methPath, options.aligner_threads, options.genome, options.genome, input_file[0]))
    os.system(cmd)
    
    with logger_mutex:
        logger.log(MESSAGE,  timestamp("MethylDackel Completed"))
    
cmdline.run(options)