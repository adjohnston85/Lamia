import os
import re
import pysam
import logging
import multiprocessing
from collections import defaultdict
from math import ceil
import subprocess
import argparse

# Setup basic configuration for logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    parser = argparse.ArgumentParser(description="Process BAM files")
    parser.add_argument('--bam_path', type=str, required=True,
                        help='Path for bam file')
    parser.add_argument('--num_threads', type=int, default=1,
                        help='Number of threads to use')
    parser.add_argument('--target_reads_per_segment', type=int, default=2000000,
                        help='Maximum number of reads per segment')
    parser.add_argument('--ref_genome_path', type=str, required=True,
                        help='Path for reference genome')
    parser.add_argument('--chromosomes', type=str, default=None,
                        help='A comma-separated list of chromosomes to analyze')
    return parser.parse_args()

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]

def check_methylation_status_using_xm(read):
    xm_tag = read.get_tag('XM')
    methylation_statuses = [char.isupper() for char in xm_tag if char in 'zZxX']
    return methylation_statuses

def clip_read_to_overlap(read, clip_start, clip_end):
    if read.query_sequence is None:
        logging.debug(f"query_sequence is None for read {read.query_name}")
        return read

    if read.query_qualities is None:
        logging.debug(f"query_qualities is None for read {read.query_name}")

    # Clip the sequences
    new_seq = read.query_sequence[clip_start:clip_end]
    new_qual = read.query_qualities[clip_start:clip_end] if read.query_qualities else None
    new_seq_len = len(new_seq)

    new_cigartuples = []
    if clip_start > 0:
        new_cigartuples.append((4, clip_start))
    if new_seq_len > 0:
        new_cigartuples.append((0, new_seq_len))
    if clip_end < len(read.query_sequence):
        new_cigartuples.append((4, len(read.query_sequence) - clip_end))

    read.query_sequence = new_seq
    read.query_qualities = new_qual
    read.cigartuples = new_cigartuples

    # Clip the XM tag if it exists
    if read.has_tag('XM'):
        xm_tag = read.get_tag('XM')
        new_xm_tag = xm_tag[clip_start:clip_end]
        read.set_tag('XM', new_xm_tag)

    return read

def mask_discrepant_cpg_sites(seq1, seq2, xm1, xm2, read1_is_reverse, read2_is_reverse, start1, end1, start2, end2):
    if seq1 is None or seq2 is None or xm1 is None or xm2 is None:
        logging.debug("One of the sequences or XM tags is None.")
        return None, None, None, None, None, None

    logging.info(f"Masking CpG sites: start1={start1}, end1={end1}, start2={start2}, end2={end2}")

    # Find the overlapping region
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    # Clip the reads to the overlapping region
    clip_start1 = overlap_start - start1
    clip_end1 = overlap_end - start1
    clip_start2 = overlap_start - start2
    clip_end2 = overlap_end - start2

    seq1 = seq1[clip_start1:clip_end1]
    seq2 = seq2[clip_start2:clip_end2]
    xm1 = xm1[clip_start1:clip_end1]
    xm2 = xm2[clip_start2:clip_end2]

    if not seq1 or not seq2:
        logging.debug("Sequences after clipping are empty.")
        logging.debug(f"Clipped Read1: {seq1}, Clipped Read2: {seq2}")
        return None, None, None, None, None, None

    masked_seq1 = list(seq1)
    masked_seq2 = list(seq2)
    masked_xm1 = list(xm1)
    masked_xm2 = list(xm2)

    len_seq = min(len(seq1), len(seq2))
    for i in range(len_seq - 1):
        # Identify CpG sites and their possible converted forms
        seq1_cpg = seq1[i:i+2] in ['CG', 'TG']
        seq2_cpg = seq2[i:i+2] in ['CG', 'CA']

        if seq1_cpg and seq2_cpg:
            # Check methylation status in XM tags
            xm1_cpg = xm1[i] in 'zZ'
            xm2_cpg = xm2[i+1] in 'zZ'

            if (xm1[i] == 'Z' and xm2[i+1] == 'Z') or (xm1[i] == 'z' and xm2[i+1] == 'z'):
                # Both reads agree on methylation status, do nothing
                continue
            else:
                # Mask the CpG sites if they do not agree
                masked_seq1[i] = 'N'
                masked_seq1[i+1] = 'N'
                masked_xm1[i] = '.'
                masked_xm1[i+1] = '.'

                masked_seq2[i] = 'N'
                masked_seq2[i+1] = 'N'
                masked_xm2[i] = '.'
                masked_xm2[i+1] = '.'

    return ''.join(masked_seq1), ''.join(masked_xm1), ''.join(masked_seq2), ''.join(masked_xm2), clip_start1, clip_end1, clip_start2, clip_end2

def process_single_end_reads(reads, base_difference_threshold=0):
    result_reads = []
    reads_by_umi = defaultdict(list)
    processed_reads = set()  # Set to keep track of processed reads

    # Group reads by UMI
    for read in reads:
        umi = read.query_name.split(':')[-1]
        reads_by_umi[umi].append(read)

    # Process reads within each UMI group
    for umi, umi_reads in reads_by_umi.items():
        # Sort reads by chromosome and coordinates
        umi_reads.sort(key=lambda r: (r.reference_name, r.reference_start))

        for i, read1 in enumerate(umi_reads):
            if read1.query_name in processed_reads:
                continue

            for j in range(i + 1, len(umi_reads)):
                read2 = umi_reads[j]
                if read2.query_name in processed_reads:
                    continue

                # Ensure reads are on the same chromosome
                if read1.reference_name == read2.reference_name:
                    # Check if starts and ends are within the base difference threshold
                    start1, end1 = read1.reference_start, read1.reference_end
                    start2, end2 = read2.reference_start, read2.reference_end

                    start_diff = abs(start1 - start2)
                    end_diff = abs(end1 - end2)

                    if start_diff <= base_difference_threshold and end_diff <= base_difference_threshold:
                        logging.info(f"\nread1: {read1.query_name}, read2: {read2.query_name}")
                        logging.info(f"\nread1: {read1.query_sequence}\nread2: {read2.query_sequence}")
                        logging.info(f"\nread1: {read1.get_tag('XM')}\nread2: {read2.get_tag('XM')}")

                        if (read1.get_tag('XG') == 'CT' and read2.get_tag('XG') == 'GA') or (read1.get_tag('XG') == 'GA' and read2.get_tag('XG') == 'CT'):
                            masked_seq1, masked_xm1, masked_seq2, masked_xm2, clip_start1, clip_end1, clip_start2, clip_end2 = mask_discrepant_cpg_sites(
                                read1.query_sequence, read2.query_sequence, 
                                read1.get_tag('XM'), read2.get_tag('XM'), 
                                read1.is_reverse, read2.is_reverse,
                                read1.reference_start, read1.reference_end,
                                read2.reference_start, read2.reference_end)

                            if masked_seq1 and masked_xm1 and masked_seq2 and masked_xm2:
                                read1.query_sequence = masked_seq1
                                read1.set_tag('XM', masked_xm1)
                                read1 = clip_read_to_overlap(read1, clip_start1, clip_end1)
                                
                                read2.query_sequence = masked_seq2
                                read2.set_tag('XM', masked_xm2)
                                read2 = clip_read_to_overlap(read2, clip_start2, clip_end2)

                                logging.info(f"\nmask1: {read1.query_sequence}\nmask2: {read2.query_sequence}")
                                logging.info(f"\nmask1: {read1.get_tag('XM')}\nmask2: {read2.get_tag('XM')}")

                                result_reads.append(read1)
                                result_reads.append(read2)
                                
                                processed_reads.add(read1.query_name)
                                processed_reads.add(read2.query_name)
                                break
                        else:
                            logging.debug(f"Failed to mask reads: {read1.query_name}, {read2.query_name}")
                            logging.debug(f"Read1: {read1.query_sequence} (start: {read1.reference_start}, end: {read1.reference_end}), Read2: {read2.query_sequence} (start: {read2.reference_start}, end: {read2.reference_end})")
                    else:
                        break  # Since reads are sorted, no further reads can overlap
                else:
                    break  # No further reads on different chromosomes

    return result_reads

def process_paired_end_reads(reads, base_difference_threshold=1):
    result_reads = []
    paired_reads_by_umi = defaultdict(list)
    processed_pairs = set()  # Set to keep track of processed pairs

    # Group paired-end reads by UMI
    for pair in reads:
        read1, read2 = pair
        umi = read1.query_name.split(':')[-1]
        paired_reads_by_umi[umi].append(pair)

    # Process paired-end reads within each UMI group
    for umi, umi_pairs in paired_reads_by_umi.items():
        # Sort pairs by chromosome and coordinates
        umi_pairs.sort(key=lambda pair: (pair[0].reference_name, pair[0].reference_start))

        for i, pair1 in enumerate(umi_pairs):
            read1a, read1b = pair1
            if read1a.query_name in processed_pairs or read1b.query_name in processed_pairs:
                continue

            for j in range(i + 1, len(umi_pairs)):
                pair2 = umi_pairs[j]
                read2a, read2b = pair2
                if read2a.query_name in processed_pairs or read2b.query_name in processed_pairs:
                    continue

                # Ensure reads are on the same chromosome
                if read1a.reference_name == read2a.reference_name:
                    # Identify the most upstream and downstream reads
                    if read1a.reference_start <= read1b.reference_start:
                        most_upstream1, most_downstream1 = read1a, read1b
                    else:
                        most_upstream1, most_downstream1 = read1b, read1a

                    if read2a.reference_start <= read2b.reference_start:
                        most_upstream2, most_downstream2 = read2a, read2b
                    else:
                        most_upstream2, most_downstream2 = read2b, read2a

                    # Check if starts of the most upstream reads and ends of the most downstream reads match
                    start_diff = abs(most_upstream1.reference_start - most_upstream2.reference_start)
                    end_diff = abs(most_downstream1.reference_end - most_downstream2.reference_end)

                    if start_diff <= base_difference_threshold and end_diff <= base_difference_threshold:
                        logging.info(f"\nread1a: {read1a.query_name}, read1b: {read1b.query_name}")
                        logging.info(f"\nread1a: {read1a.query_sequence}\nread1b: {read1b.query_sequence}")
                        logging.info(f"\nread1a: {read1a.get_tag('XM')}\nread1b: {read1b.get_tag('XM')}")
                        logging.info(f"\nread2a: {read2a.query_name}, read2b: {read2b.query_name}")
                        logging.info(f"\nread2a: {read2a.query_sequence}\nread2b: {read2b.query_sequence}")
                        logging.info(f"\nread2a: {read2a.get_tag('XM')}\nread2b: {read2b.get_tag('XM')}")

                        if (read1a.get_tag('XG') == 'CT' and read2a.get_tag('XG') == 'GA') or (read1a.get_tag('XG') == 'GA' and read2a.get_tag('XG') == 'CT'):
                            masked_seq1a, masked_xm1a, masked_seq2a, masked_xm2a, clip_start1a, clip_end1a, clip_start2a, clip_end2a = mask_discrepant_cpg_sites(
                                read1a.query_sequence, read2a.query_sequence, 
                                read1a.get_tag('XM'), read2a.get_tag('XM'), 
                                read1a.is_reverse, read2a.is_reverse,
                                read1a.reference_start, read1a.reference_end,
                                read2a.reference_start, read2a.reference_end)

                            masked_seq1b, masked_xm1b, masked_seq2b, masked_xm2b, clip_start1b, clip_end1b, clip_start2b, clip_end2b = mask_discrepant_cpg_sites(
                                read1b.query_sequence, read2b.query_sequence, 
                                read1b.get_tag('XM'), read2b.get_tag('XM'), 
                                read1b.is_reverse, read2b.is_reverse,
                                read1b.reference_start, read1b.reference_end,
                                read2b.reference_start, read2b.reference_end)

                            if masked_seq1a and masked_xm1a and masked_seq2a and masked_xm2a and masked_seq1b and masked_xm1b and masked_seq2b and masked_xm2b:
                                read1a.query_sequence = masked_seq1a
                                read1a.set_tag('XM', masked_xm1a)
                                read1a = clip_read_to_overlap(read1a, clip_start1a, clip_end1a)

                                read2a.query_sequence = masked_seq2a
                                read2a.set_tag('XM', masked_xm2a)
                                read2a = clip_read_to_overlap(read2a, clip_start2a, clip_end2a)

                                read1b.query_sequence = masked_seq1b
                                read1b.set_tag('XM', masked_xm1b)
                                read1b = clip_read_to_overlap(read1b, clip_start1b, clip_end1b)

                                read2b.query_sequence = masked_seq2b
                                read2b.set_tag('XM', masked_xm2b)
                                read2b = clip_read_to_overlap(read2b, clip_start2b, clip_end2b)

                                logging.info(f"\nmask1a: {read1a.query_sequence}\nmask2a: {read2a.query_sequence}")
                                logging.info(f"\nmask1a: {read1a.get_tag('XM')}\nmask2a: {read2a.get_tag('XM')}")
                                logging.info(f"\nmask1b: {read1b.query_sequence}\nmask2b: {read2b.query_sequence}")
                                logging.info(f"\nmask1b: {read1b.get_tag('XM')}\nmask2b: {read2b.get_tag('XM')}")

                                result_reads.extend([read1a, read2a, read1b, read2b])
                                
                                processed_pairs.update([read1a.query_name, read2a.query_name, read1b.query_name, read2b.query_name])
                                break
                        else:
                            logging.debug(f"Failed to mask reads: {read1a.query_name}, {read2a.query_name}, {read1b.query_name}, {read2b.query_name}")
                            logging.debug(f"Read1a: {read1a.query_sequence} (start: {read1a.reference_start}, end: {read1a.reference_end}), Read2a: {read2a.query_sequence} (start: {read2a.reference_start}, end: {read2a.reference_end})")
                            logging.debug(f"Read1b: {read1b.query_sequence} (start: {read1b.reference_start}, end: {read1b.reference_end}), Read2b: {read2b.query_sequence} (start: {read2b.reference_start}, end: {read2b.reference_end})")
                    else:
                        break  # Since reads are sorted, no further reads can overlap
                else:
                    break  # No further reads on different chromosomes

    return result_reads

def fetch_reads(bam_file, contig, start, end):
    single_end_reads = []
    paired_end_reads = defaultdict(list)

    for read in bam_file.fetch(contig, start, end):
        if read.is_unmapped:
            continue

        if read.is_secondary or read.is_supplementary:
            continue

        if read.is_paired:
            paired_end_reads[read.query_name].append(read)
        else:
            single_end_reads.append(read)

    paired_end_reads_list = [reads for reads in paired_end_reads.values() if len(reads) == 2]
    return single_end_reads, paired_end_reads_list

def process_paired_end_reads(reads, base_difference_threshold=1):
    result_reads = []

    for pair in reads:
        read1, read2 = pair

        # Ensure reads are on the same chromosome
        if read1.reference_name == read2.reference_name:
            # Check if starts and ends are within the base difference threshold
            start1, end1 = read1.reference_start, read1.reference_end
            start2, end2 = read2.reference_start, read2.reference_end

            start_diff = abs(start1 - start2)
            end_diff = abs(end1 - end2)

            if start_diff <= base_difference_threshold and end_diff <= base_difference_threshold:
                logging.info(f"read1: {read1.query_sequence}\n read2: {read2.query_sequence}")
                logging.info(f"read1: {read1.get_tag('XM')}\n read2: {read2.get_tag('XM')}")

                if (read1.get_tag('XG') == 'CT' and read2.get_tag('XG') == 'GA') or (read1.get_tag('XG') == 'GA' and read2.get_tag('XG') == 'CT'):
                    masked_seq1, masked_xm1, masked_seq2, masked_xm2, clip_start1, clip_end1, clip_start2, clip_end2 = mask_discrepant_cpg_sites(
                        read1.query_sequence, read2.query_sequence, 
                        read1.get_tag('XM'), read2.get_tag('XM'), 
                        read1.is_reverse, read2.is_reverse,
                        read1.reference_start, read1.reference_end,
                        read2.reference_start, read2.reference_end)

                    if masked_seq1 and masked_xm1 and masked_seq2 and masked_xm2:
                        read1.query_sequence = masked_seq1
                        read1.set_tag('XM', masked_xm1)
                        read1 = clip_read_to_overlap(read1, clip_start1, clip_end1)
                        
                        read2.query_sequence = masked_seq2
                        read2.set_tag('XM', masked_xm2)
                        read2 = clip_read_to_overlap(read2, clip_start2, clip_end2)
                        
                        logging.info(f"read1: {read1.query_sequence}\n read2: {read2.query_sequence}")
                        logging.info(f"read1: {read1.get_tag('XM')}\n read2: {read2.get_tag('XM')}")

                        result_reads.append(read1)
                        result_reads.append(read2)
                    else:
                        logging.debug(f"Failed to mask reads: {read1.query_name}, {read2.query_name}")
                        logging.debug(f"Read1: {read1.query_sequence} (start: {read1.reference_start}, end: {read1.reference_end}), Read2: {read2.query_sequence} (start: {read2.reference_start}, end: {read2.reference_end})")
                else:
                    logging.debug(f"Reads {read1.query_name} and {read2.query_name} do not have the correct XG tags for pairing.")
            else:
                logging.debug(f"Paired reads {read1.query_name} and {read2.query_name} have different start/end positions exceeding the base difference threshold.")
        else:
            logging.debug(f"Paired reads {read1.query_name} and {read2.query_name} are on different chromosomes.")
    return result_reads

def process_bam_segment(args):
    (bam_path, phase_path, contig, start, end, progress_dict, process_id) = args

    logging.info(f"Start processing {contig} from {start} to {end} in process {process_id}")
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    header = bam_file.header.copy()

    out_phase_bam = pysam.AlignmentFile(phase_path, "wb", header=header)

    single_end_reads, paired_end_reads = fetch_reads(bam_file, contig, start, end)

    processed_single_end_reads = process_single_end_reads(single_end_reads)
    processed_paired_end_reads = process_paired_end_reads(paired_end_reads)

    for read in processed_single_end_reads:
        out_phase_bam.write(read)

    for read in processed_paired_end_reads:
        out_phase_bam.write(read)

    out_phase_bam.close()
    bam_file.close()
    logging.info(f"Finished processing {contig} from {start} to {end} in process {process_id}.")
    progress_dict[process_id] = 100

def calculate_segments_by_length(bam_path, num_threads, target_reads_per_segment, chromosomes):
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
                continue

            proportion = length / total_length
            num_segments = max(num_threads, int(ceil(proportion * target_reads_per_segment)))
            segments[contig] = num_segments
    return segments

def calculate_segments_by_reads(bam_path, num_threads, target_reads_per_segment, chromosomes):
    segments = calculate_segments_by_length(
        bam_path, num_threads, target_reads_per_segment, chromosomes)
    bam = pysam.AlignmentFile(bam_path, "rb")
    adjusted_segments = []

    for contig in segments:
        logging.info(f"Segmenting {contig}")
        contig_length = bam.get_reference_length(contig)
        read_positions = sorted([r.reference_start for r in bam.fetch(contig) if not r.is_unmapped and not r.is_secondary and not r.is_supplementary])

        if not read_positions:
            logging.warning(f"No reads available in {contig} for processing.")
            continue

        if len(read_positions) < target_reads_per_segment:
            logging.info(f"Contig {contig} has fewer reads ({len(read_positions)}) than target per segment ({target_reads_per_segment}). Processing as a single segment.")
            adjusted_segments.append((contig, read_positions[0], contig_length))
        else:
            contig_read_threshold = len(read_positions) // segments[contig]
            segment_start = 0
            for i in range(segments[contig]):
                if i < segments[contig] - 1:
                    segment_end = read_positions[min((i + 1) * contig_read_threshold - 1, len(read_positions) - 1)]
                else:
                    segment_end = contig_length
                adjusted_segments.append((contig, segment_start, segment_end))
                segment_start = segment_end + 1

    return adjusted_segments

def split_bam_processing(bam_path, num_threads, target_reads_per_segment, chromosome_group, merged_duplex_path):
    logging.info(f"Segmenting BAM reads into contig segments for multicore processing using {num_threads} threads")
    segments = calculate_segments_by_reads(bam_path, num_threads, target_reads_per_segment, chromosome_group)
    jobs = []
    manager = multiprocessing.Manager()
    progress_dict = manager.dict()
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        header = bam.header
        job_id = 0

        for segment in segments:
            contig, start, end = segment 
            phase_path = f"{merged_duplex_path.replace('.bam', '')}_duplex_temp{job_id}.bam"
            job = (bam_path, phase_path, contig, start, end, progress_dict, job_id)
            logging.info(f"Process {job_id}. Segment of {contig}: {start}-{end}")
            jobs.append(job)
            job_id += 1

        with multiprocessing.Pool(processes=num_threads) as pool:
            pool.map(process_bam_segment, jobs)

        with pysam.AlignmentFile(merged_duplex_path, "wb", header=header) as out_duplex_bam:
            for job in jobs:
                try:
                    with pysam.AlignmentFile(job[1], "rb") as temp_duplex_bam:
                        for read in temp_duplex_bam:
                            out_duplex_bam.write(read)
                except Exception as e:
                    logging.error(f"Error reading temporary BAM file {job[1]}: {e}")
                finally:
                    if os.path.exists(job[1]):
                        os.remove(job[1])

def get_contigs(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        contigs = bam_file.references
    return list(contigs)

def get_reads_per_chromosome(bam_path, chromosomes):
    logging.info("Calculating number of reads per chromosome")
    reads_per_chromosome = {}
    total_reads = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom in chromosomes:
            read_count = bam.count(chrom)
            reads_per_chromosome[chrom] = read_count
            total_reads += read_count

    for chrom, count in reads_per_chromosome.items():
        logging.info(f"Chromosome {chrom}: {count} reads")

    logging.info(f"Total number of reads across all specified chromosomes: {total_reads}")
    return reads_per_chromosome

def group_chromosomes_by_reads(reads_per_chromosome, num_threads, target_reads_per_segment):
    sorted_chromosomes = sorted(reads_per_chromosome, key=reads_per_chromosome.get, reverse=True)
    groups = []
    current_group = []
    current_reads = 0

    for chrom in sorted_chromosomes:
        if current_reads + reads_per_chromosome[chrom] > target_reads_per_segment * num_threads:
            if current_group:
                groups.append(current_group)
                logging.info(f"Formed group with {len(current_group)} chromosomes, totaling {current_reads} reads: {current_group}")
            current_group = [chrom]
            current_reads = reads_per_chromosome[chrom]
        else:
            current_group.append(chrom)
            current_reads += reads_per_chromosome[chrom]

    if current_group:
        groups.append(current_group)
        logging.info(f"Formed group with {len(current_group)} chromosomes, totaling {current_reads} reads: {current_group}")

    return groups

def main():
    args = parse_args()
    chromosomes = args.chromosomes.split(',') if args.chromosomes else []

    if "core_human" in chromosomes:
        chromosomes = ['chr'+str(i) for i in range(1, 23)] + ['chrX', 'chrY']

    bam_files = [("CT", args.bam_path)]

    for bam_type, bam_path in bam_files:
        logging.info(f"Starting {bam_type} processing")
        
        if not chromosomes:
            chromosomes = get_contigs(bam_path)
        
        chromosome_bam_files_duplex = []
        
        reads_per_chromosome = get_reads_per_chromosome(bam_path, chromosomes)
        chromosome_groups = group_chromosomes_by_reads(reads_per_chromosome, args.num_threads, args.target_reads_per_segment)

        for group_index, group in enumerate(chromosome_groups):
            group_id = f"group{group_index}"
            merged_duplex_group_path = bam_path.replace('.bam', f'_{group_id}_duplex.bam')

            split_bam_processing(bam_path, 
                                 args.num_threads, 
                                 args.target_reads_per_segment, 
                                 group,
                                 merged_duplex_group_path)

            chromosome_bam_files_duplex.append(merged_duplex_group_path)

        merged_duplex_path = bam_path.replace('.bam', '_duplex.bam')

        subprocess.run(['samtools', 'merge', '-f', '-@', str(args.num_threads), merged_duplex_path] 
                       + chromosome_bam_files_duplex, check=True)

        for file_path in chromosome_bam_files_duplex:
            os.remove(file_path)

    logging.info("All processing completed")

if __name__ == "__main__":
    main()

