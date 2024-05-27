import os
import re
import pysam
import logging
import multiprocessing
from multiprocessing import Pool
from intervaltree import IntervalTree
from math import ceil
import time
import argparse
import subprocess

# Setup basic configuration for logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]

def write_cpg_sites_to_contig_beds(cpg_trees, bed_path_prefix):
    """
    Writes CpG sites to separate BED files for each contig.
    Args:
        cpg_trees (dict): Dictionary containing contig names as keys and IntervalTrees as values.
        bed_path_prefix (str): Path prefix for output BED files, including directory and filename prefix.
    """
    for contig, tree in cpg_trees.items():
        bed_file_path = f"{bed_path_prefix}{contig}.bed"
        with open(bed_file_path, 'w') as bed_file:
            sites = sorted(tree)
            for interval in sites:
                bed_file.write(f"{contig}\t{interval.begin + 1}\t{interval.end + 1}\n")
        logging.info(f"CpG sites for {contig} saved to {bed_file_path}")

def update_md_tag(read, ref_genome):
    """
    Update the MD:Z tag of a read after bases have been reverted.
    This method recalculates the MD tag based on the read's current sequence
    and its alignment to the reference.
    """
    aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)
    md_string = []
    num_matches = 0  # Counter for consecutive matches

    for qpos, rpos, rbase in aligned_pairs:
        if rpos is None:
            continue  # Skip insertions in the read
        if qpos is None:
            continue  # Skip deletions in the read
        ref_base = ref_genome.fetch(read.reference_name, rpos, rpos + 1).upper()
        read_base = read.query_sequence[qpos]
        if read_base == ref_base:
            num_matches += 1  # Increment match counter
        else:
            if num_matches > 0:
                md_string.append(str(num_matches))  # Record number of matches
                num_matches = 0
            md_string.append(ref_base)  # Record the mismatched reference base

    if num_matches > 0:
        md_string.append(str(num_matches))  # Ensure final matches are recorded

    read.set_tag('MD', ''.join(md_string))  # Set the new MD tag on the read

def process_bam_segment(args):
    (bam_path, bed_file_path, revert_path, masked_path, revert_to, check_base, 
     contig, start, end, ref_genome_path, progress_dict, process_id) = args
    
    cpg_sites = IntervalTree()
    try:
        with open(bed_file_path, 'r') as file:
            for line in file:
                parts = line.strip().split()
                site_start = int(parts[1]) - 1  # BED is zero-based, start inclusive
                site_end = int(parts[2]) - 1    # BED is zero-based, end exclusive

                # Check if the CpG site is beyond the segment of interest
                if site_start > end + 300:
                    break  # Exit the loop as subsequent sites will also be out of range
                
                # Only add CpG sites that are within the segment being processed
                if site_start >= start and site_end <= end:
                    cpg_sites.addi(site_start, site_end)

    except FileNotFoundError:
        logging.error(f"CpG BED file not found for {contig}")

    ref_genome = pysam.FastaFile(ref_genome_path)
    logging.info(f"Start processing {contig} from {start} to {end} in process {process_id}")
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        total_reads = bam.count(contig=contig, start=start, end=end)
        header = bam.header.copy()
        read_count = 0
        reverted_count = 0
        percent_step = 25
        next_percent_threshold = total_reads * (percent_step / 100)
        percent_increment = next_percent_threshold

        start_time = time.time()

        out_revert_bam = pysam.AlignmentFile(revert_path, "wb", header=header)
        out_masked_bam = pysam.AlignmentFile(masked_path, "wb", header=header)
        
        for read in bam.fetch(contig, start, end):
            if read.reference_start < start or read.reference_start >= end:
                continue

            revert_seq = list(read.query_sequence)
            masked_seq = list(read.query_sequence)
            logging.debug(f"seq: {revert_seq}")
            for idx, pos in enumerate(read.get_reference_positions(full_length=True)):
                logging.debug(f"seq_base: {revert_seq[idx]}")
                if pos is None:
                    continue
                ref_base = ref_genome.fetch(bam.references[read.reference_id], pos, pos + 1).upper()
                logging.debug(f"ref_base: {ref_base}")

                # if base is part of CpG site 
                overlapping_segments = cpg_sites.overlap(pos, pos+1)
                if overlapping_segments:                  
                    CpG_coords = next(iter(overlapping_segments))
                    CpG_start = CpG_coords.begin
                    CpG_end = CpG_coords.end
                    logging.debug(f"CpG overlap: {pos}")

                    # if the reads is C to T
                    if revert_to == 'C':
                        if pos == CpG_end:
                            # Mask the G base of the CpG site in the reversion BAM, 
                            # as this conflicts with the methylation calling on the opposite strand
                            revert_seq[idx] = 'N'
                        elif pos == CpG_start and masked_seq[idx].upper() in ['C', 'T']:
                            # Mask the C base of the CpG in the masked BAM if it is C or T,
                            # as deamination makes these bases uninformtive for genomic variant calling
                            masked_seq[idx] = 'N'
                    elif revert_to == 'G':
                        if pos == CpG_end:
                            revert_seq[idx] = 'N'
                        elif pos == start and masked_seq[idx].upper() in ['G', 'A']:
                            masked_seq[idx] = 'N'

                # Otherwise if a non-CpG site matches expected conversion by deamination,
                # revert this base to the reference, and mask it for genomic variant calling
                elif ref_base == revert_to:
                    if revert_seq[idx].upper() == check_base:
                        logging.debug(f"revert_to: {revert_to} and masked")
                        revert_seq[idx] = revert_to
                        reverted_count += 1
                        masked_seq[idx] = 'N'
                    elif revert_seq[idx].upper() == check_base:
                        logging.debug(f"masked")
                        masked_seq[idx] = 'N'

            read.query_sequence = ''.join(revert_seq)
            update_md_tag(read, ref_genome)  
            out_revert_bam.write(read)
            
            read.query_sequence = ''.join(masked_seq)
            update_md_tag(read, ref_genome)  
            out_masked_bam.write(read)
                
            read_count += 1

            if read_count >= next_percent_threshold:
                current_percent = (read_count / total_reads) * 100
                elapsed_time = time.time() - start_time
                reads_per_second = read_count / elapsed_time
                logging.info(f"Proc {process_id} {contig} {start}-{end}: {read_count}/{total_reads} ({current_percent:.0f}%) [{reads_per_second:.0f} reads/s]")
                next_percent_threshold += percent_increment
        
        out_revert_bam.close()
        out_masked_bam.close()

        logging.info(f"Finished processing {contig} from {start} to {end} in process {process_id}. Reads processed: {read_count}, reverted: {reverted_count}")
        progress_dict[process_id] = 100

def calculate_segments_by_length(bam_path, num_threads, target_reads_per_segment, chromosomes, CpG_search):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        total_length = sum(
            l for r, l in zip(bam.references, bam.lengths) 
            if chromosomes is None or r in chromosomes)
        segments = {}
        for contig, length in zip(bam.references, bam.lengths):
            if chromosomes is not None and contig not in chromosomes:
                continue
            contig_reads = bam.count(reference=contig)
            if contig_reads == 0:
                logging.info(f"Skipping {contig} due to zero reads.")
                continue  # Skip this contig entirely if no reads are present

            if CpG_search:
                proportion = length / total_length
                num_segments = num_threads
            else:
                total_reads = sum(
                    bam.count(reference=r) for r in (bam.references if chromosomes is None else chromosomes))
                total_segments = int(ceil(total_reads / target_reads_per_segment))
                num_segments = max(num_threads, total_segments)
                proportion = contig_reads / total_reads
            segments[contig] = max(1, int(ceil(proportion * num_segments)))
    return segments

def calculate_segments_by_reads(bam_path, num_threads, target_reads_per_segment, chromosome, CpG_search):
    segments = calculate_segments_by_length(
        bam_path, num_threads, target_reads_per_segment, chromosome, CpG_search)
    bam = pysam.AlignmentFile(bam_path, "rb")
    adjusted_segments = []

    for contig in segments:
        logging.info(f"Segmenting {contig}")
        contig_length = bam.get_reference_length(contig)
        contig_reads = bam.count(reference=contig)
        read_positions = sorted([r.reference_start for r in bam.fetch(contig)])

        if not read_positions:
            logging.warning(f"No reads available in {contig} for processing.")
            continue

        # Adjust segment calculation if total reads in the contig are less than target_reads_per_segment
        if contig_reads < target_reads_per_segment:
            logging.info(f"Contig {contig} has fewer reads ({contig_reads}) than target per segment ({target_reads_per_segment}). Processing as a single segment.")
            adjusted_segments.append((contig, read_positions[0], contig_length))
        else:
            contig_read_threshold = len(read_positions) // segments[contig]
            segment_start = 0
            for i in range(segments[contig]):
                if i < segments[contig] - 1:
                    segment_end = read_positions[min((i + 1) * contig_read_threshold - 1, len(read_positions) - 1)]
                else:
                    # Extend the last segment up to the length of the contig plus the maximum read length
                    segment_end = contig_length
                adjusted_segments.append((contig, segment_start, segment_end))
                segment_start = segment_end + 1

    return adjusted_segments

def split_bam_processing(bam_path, ref_genome_path, num_threads, revert_to, check_base, bed_file_paths,
                         target_reads_per_segment, chromosome_group, merged_reverted_path, merged_masked_path):
    logging.info(f"Segmenting BAM reads into contig segments for multicore processing using {num_threads} threads")
    segments = calculate_segments_by_reads(bam_path, num_threads, target_reads_per_segment, chromosome_group, False)
    jobs = []
    manager = multiprocessing.Manager()
    progress_dict = manager.dict()
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        header = bam.header
        job_id = 0

        for segment in segments:
            contig, start, end = segment 
            if contig in bed_file_paths:
                revert_path = f"{merged_reverted_path.replace('.bam', '')}_reverted_temp{job_id}.bam"
                masked_path = f"{merged_reverted_path.replace('.bam', '')}_masked_temp{job_id}.bam"
                job = (bam_path, bed_file_paths[contig], revert_path, masked_path, revert_to, check_base,
                       contig, start, end, ref_genome_path, progress_dict, job_id)
                logging.info(f"Process {job_id}. Segment of {contig}: {start}-{end}")
                jobs.append(job)
                job_id += 1

        with multiprocessing.Pool(processes=num_threads) as pool:
            pool.map(process_bam_segment, jobs)

        with pysam.AlignmentFile(merged_reverted_path, "wb", header=header) as out_reverted_bam,\
             pysam.AlignmentFile(merged_masked_path, "wb", header=header) as out_masked_bam:
              
            for job in jobs:
                with pysam.AlignmentFile(job[2], "rb") as temp_reverted_bam:
                    for read in temp_reverted_bam:
                        out_reverted_bam.write(read)

                os.remove(job[2])

                with pysam.AlignmentFile(job[3], "rb") as temp_masked_bam:
                    for read in temp_masked_bam:
                        out_masked_bam.write(read)

                os.remove(job[3])

def count_total_cpg_sites(cpg_trees):
    total_cpg_sites = 0
    for tree in cpg_trees.values():
        total_cpg_sites += len(tree)
    return total_cpg_sites

def load_cpg_sites_from_contig_bed(bed_file_path):
    logging.info(f"Loading CpG sites from {bed_file_path}")
    cpg_tree = IntervalTree()
    if os.path.exists(bed_file_path):
        with open(bed_file_path, 'r') as file:
            for line in file:
                parts = line.strip().split()
                start = int(parts[1]) - 1
                end = int(parts[2]) - 1
                cpg_tree.addi(start, end)
    logging.info(f"Completed loading {len(cpg_tree)} CpG sites from {bed_file_path}")
    return cpg_tree

def load_cpg_sites_parallel(ref_genome_path, cpg_bed_path_prefix, num_threads, chromosomes):
    logging.info("Loading CpG sites from BED files:")
    contigs = pysam.FastaFile(ref_genome_path).references if not chromosomes else chromosomes
     
    pool = multiprocessing.Pool(processes=num_threads)
    bed_file_paths = [f"{cpg_bed_path_prefix}{contig}.bed" for contig in contigs]
    trees = pool.map(load_cpg_sites_from_contig_bed, bed_file_paths)

    cpg_trees = {contig: tree for contig, tree in zip(contigs, trees)}
    pool.close()
    pool.join()
    total_cpg_sites = count_total_cpg_sites(cpg_trees)
    logging.info(f"Completed loading {total_cpg_sites} CpG sites across all chromosomes")
    return cpg_trees

def get_contigs(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        contigs = bam_file.references
    return list(contigs)

def load_cpg_sites_for_segment(args):
    contig, sequence, start, end, process_id = args
    cpg_sites = IntervalTree()
    segment_length = end - start
    total_positions = segment_length - 1  # Total positions to check for CpG sites
    progress_step = total_positions / 4  # Log progress every 25%
    next_progress_threshold = progress_step
    count_positions = 0

    logging.info(f"Proc {process_id}: Starting CpG detection for {contig} from {start} to {end}")
    pos = start
    while pos < end:  # Ensure it goes up to, but not including, `end` which should be the length of the segment
        if pos + 1 < start + len(sequence):  # Check there's a next base
            dinucleotide = sequence[pos-start:pos-start+2]
            logging.debug(f"Checking position {pos}: {dinucleotide}")
            if dinucleotide.upper() == 'CG':
                cpg_sites.addi(pos, pos + 2)  # Adjust if necessary
                logging.debug(f"Detected CpG site at position {pos}: {dinucleotide}")
                pos += 2
            else:
                pos += 1
        else:
            break
 
        count_positions += 1
        if count_positions >= next_progress_threshold:
            percent_done = (count_positions / total_positions) * 100
            logging.info(
                f"Proc {process_id} {contig} segment {start}-{end}: "
                f"Processed {count_positions}/{total_positions} positions ({percent_done:.0f}%)"
            )
            next_progress_threshold += progress_step

    logging.info(f"Proc {process_id}: Finished CpG detection for {contig} from {start} to {end}")
    return contig, start, end, cpg_sites

def load_cpg_sites_from_genome(ref_genome_path, bam_path, num_threads, target_reads_per_segment, chromosome):
    nchunks = calculate_segments_by_length(bam_path, num_threads, target_reads_per_segment, chromosome, True)
    segments = []
    process_id = 0
    ref_genome = pysam.FastaFile(ref_genome_path)

    for contig, chunks in nchunks.items():
        length = ref_genome.get_reference_length(contig)
        segment_length = length // chunks
        for i in range(chunks):
            start = i * segment_length
            end = (i + 1) * segment_length + 1 if i != chunks - 1 else length
            sequence = ref_genome.fetch(contig, start, end)
            logging.info(f"Fetched sequence for {contig} between {start}-{end}")
            segments.append((contig, sequence, start, end, process_id))
            process_id += 1

    with multiprocessing.Pool(processes=num_threads) as pool:
        results = pool.map(load_cpg_sites_for_segment, segments)

    cpg_trees = {}
    for contig, start, end, sites in results:
        if contig not in cpg_trees:
            cpg_trees[contig] = IntervalTree()
        
        # Debug log step 
        logging.info(f"Number of CpG sites for {contig} between {start}-{end}: {len(sites)}")
        
        # Update the existing tree instead of overwriting it.
        cpg_trees[contig].update(sites)

    return cpg_trees

def check_cpg_bed_files_exist(cpg_bed_path_prefix, ref_genome_path, chromosomes):
    contigs = pysam.FastaFile(ref_genome_path).references if not chromosomes else chromosomes
    for contig in contigs:
        bed_file_path = f"{cpg_bed_path_prefix}{contig}.bed"
        if not os.path.exists(bed_file_path):
            logging.error(f"Missing CpG BED file for contig {contig}: {bed_file_path}")
            return False
    return True

def get_reads_per_chromosome(bam_path, chromosomes):
    logging.info("Calculating number of reads per chromosome")
    reads_per_chromosome = {}
    total_reads = 0  # Initialize total reads counter

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom in chromosomes:
            read_count = bam.count(chrom)
            reads_per_chromosome[chrom] = read_count
            total_reads += read_count  # Accumulate total reads

    # Log the read counts for each chromosome and the total
    for chrom, count in reads_per_chromosome.items():
        logging.info(f"Chromosome {chrom}: {count} reads")

    logging.info(f"Total number of reads across all specified chromosomes: {total_reads}")
    return reads_per_chromosome

def group_chromosomes_by_reads(reads_per_chromosome, num_threads, target_reads_per_segment):
    # Sort chromosomes by read count for potentially more balanced groups
    sorted_chromosomes = sorted(reads_per_chromosome, key=reads_per_chromosome.get, reverse=True)
    groups = []
    current_group = []
    current_reads = 0

    # Grouping logic
    for chrom in sorted_chromosomes:
        if current_reads + reads_per_chromosome[chrom] > target_reads_per_segment * num_threads:
            if current_group:
                groups.append(current_group)
                # Log the group details when it's added
                logging.info(f"Formed group with {len(current_group)} chromosomes, totaling {current_reads} reads: {current_group}")
            current_group = [chrom]
            current_reads = reads_per_chromosome[chrom]
        else:
            current_group.append(chrom)
            current_reads += reads_per_chromosome[chrom]

    if current_group:
        groups.append(current_group)
        # Log the last group details
        logging.info(f"Formed group with {len(current_group)} chromosomes, totaling {current_reads} reads: {current_group}")

    return groups

def parse_args():
    parser = argparse.ArgumentParser(description="Process BAM files")
    parser.add_argument('--ct_bam_path', type=str, required=True,
                        help='Path for CT bam file')
    parser.add_argument('--ga_bam_path', type=str, required=True,
                        help='Path for GA bam file')
    parser.add_argument('--num_threads', type=int, default=1,
                        help='Number of threads to use')
    parser.add_argument('--target_reads_per_segment', type=int, default=10000,
                        help='Maximum number of reads per segment')
    parser.add_argument('--ref_genome_path', type=str, required=True,
                        help='Path for reference genome')
    parser.add_argument('--cpg_bed_path_prefix', type=str, default='',
                        help='Path for CpG bed file')
    parser.add_argument('--chromosomes', type=str, default=None,
                        help='A semicolon-separated list of chromosomes to analyze')
    parser.add_argument('--make_cpg_bed', action='store_true',
                        help='Create a CpG bed file if set, defaults to False')
    return parser.parse_args()

def main():
    args = parse_args()
    chromosomes = args.chromosomes.split(',') if args.chromosomes else None

    if "core_human" in chromosomes:
        chromosomes = ['chr'+str(i) for i in range(1, 23)] + ['chrX', 'chrY']

    # Dictionary to store BED file paths and list to track missing BED files
    bed_file_paths = {}
    missing_chromosomes = []

    for chrom in chromosomes:
        bed_path = f"{args.cpg_bed_path_prefix}{chrom}.bed"
        if os.path.exists(bed_path):
            bed_file_paths[chrom] = bed_path
        else:
            logging.error(f"Missing CpG BED file for chromosome {chrom}: {bed_path}")
            missing_chromosomes.append(chrom)

    # Check if there are missing BED files and process accordingly
    if missing_chromosomes:
        logging.info("Locating CpG sites from reference genome for missing chromosomes")
        cpg_trees = load_cpg_sites_from_genome(args.ref_genome_path, args.ct_bam_path, args.num_threads,
                                               args.target_reads_per_segment, missing_chromosomes)
        if args.make_cpg_bed:
            write_cpg_sites_to_contig_beds(cpg_trees, args.cpg_bed_path_prefix)
            # Update bed_file_paths with newly created BED files
            for chrom in missing_chromosomes:
                bed_path = f"{args.cpg_bed_path_prefix}{chrom}.bed"
                bed_file_paths[chrom] = bed_path

    bam_files = [("CT", args.ct_bam_path), ("GA", args.ga_bam_path)]

    for bam_type, bam_path in bam_files:
        logging.info(f"Starting {bam_type} processing")
        
        if not chromosomes:
            chromosomes = get_contigs(bam_path)
        
        chromosome_bam_files_reverted = []
        chromosome_bam_files_masked = []
        
        reads_per_chromosome = get_reads_per_chromosome(bam_path, chromosomes)
        chromosome_groups = group_chromosomes_by_reads(reads_per_chromosome, args.num_threads, args.target_reads_per_segment)

        for group_index, group in enumerate(chromosome_groups):
            revert_to = 'C' if bam_type == "CT" else 'G'
            check_base = 'T' if bam_type == "CT" else 'A'
            
            group_id = f"group{group_index}"
            merged_reverted_group_path = bam_path.replace('.bam', f'_{group_id}_reverted.bam')
            merged_masked_group_path = bam_path.replace('.bam', f'_{group_id}_masked.bam')

            split_bam_processing(bam_path, 
                                 args.ref_genome_path, 
                                 args.num_threads, 
                                 revert_to, 
                                 check_base, 
                                 bed_file_paths, 
                                 args.target_reads_per_segment, 
                                 group,
                                 merged_reverted_group_path,
                                 merged_masked_group_path)

            # Keep track of the created files
            chromosome_bam_files_reverted.append(merged_reverted_group_path)
            chromosome_bam_files_masked.append(merged_masked_group_path)

        # Merge into final bam files
        merged_reverted_path = bam_path.replace('.bam', '_reverted.bam')
        merged_masked_path = bam_path.replace('.bam', '_masked.bam')

        subprocess.run(['samtools', 'merge', '-f', '-@', str(args.num_threads), merged_reverted_path] 
                       + chromosome_bam_files_reverted, check=True)
        subprocess.run(['samtools', 'merge', '-f',  '-@', str(args.num_threads), merged_masked_path] 
                       + chromosome_bam_files_masked, check=True)

        for file_path in chromosome_bam_files_reverted:
            os.remove(file_path)

        for file_path in chromosome_bam_files_masked:
            os.remove(file_path)

    logging.info("All processing completed")

if __name__ == "__main__":
    main()
