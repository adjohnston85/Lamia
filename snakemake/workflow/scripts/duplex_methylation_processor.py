#!/usr/bin/env python

'''
Title: duplex_methylation_processor.py
Date: 02/06/2024
Author: Andrew D Johnston
Description:
    This script processes BAM files containing methylation sequenced consensus reads 
    and conducts duplex matching. Clips these duplex read pairs to the overlap regions
    and masks discrepant CpG sites in bisulfite sequencing data. The script segments the 
    BAM file for parallel processing, handles both single-end and paired-end reads, 
    and writes the processed reads to a new BAM file. The script requires a BAM 
    file with bisulfite sequencing tags (XM and XG).

List of functions:
    main()
    parse_args()
    natural_sort_key()
    check_methylation_status_using_xm()
    clip_read_to_overlap()
    mask_discrepant_cpg_sites()
    process_single_end_reads()
    process_paired_end_reads()
    fetch_reads()
    process_bam_segment()
    calculate_segments_by_length()
    calculate_segments_by_reads()
    split_bam_processing()
    get_contigs()
    get_reads_per_chromosome()
    group_chromosomes_by_reads()

Procedure:
    1. Parse command-line arguments
    2. Get the list of chromosomes from the BAM file or specified arguments
    3. Calculate the number of reads per chromosome
    4. Group chromosomes for balanced parallel processing
    5. Segment the BAM file for multicore processing
    6. Process each segment to clip reads to overlaps and mask discrepant CpG sites
    7. Merge the processed segments into a final BAM file

Usage:
    ./process_bam_files.py [-h] --bam_path BAM_PATH [--num_threads NUM_THREADS] 
    [--target_reads_per_segment TARGET_READS_PER_SEGMENT] --ref_genome_path REF_GENOME_PATH 
    [--chromosomes CHROMOSOMES]

Example:
    ./process_bam_files.py --bam_path input.bam --num_threads 4 --target_reads_per_segment 2000000 
    --ref_genome_path reference.fa --chromosomes chr1,chr2,chrX
'''

import os
import re
import pysam  # A Python module for reading and manipulating BAM files
import logging  # For logging information
import multiprocessing  # For parallel processing
from collections import defaultdict  # For creating default dictionaries
from math import ceil  # For ceiling function
import subprocess  # For running subprocesses (e.g., samtools merge)
import argparse  # For command-line argument parsing

# Setup basic configuration for logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    """
    Parse command-line arguments.
    """
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
    """
    Key function for natural sorting of strings containing numbers.
    """
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]

def check_methylation_status_using_xm(read):
    """
    Check methylation status using the 'XM' tag in the read.
    Returns a list of boolean values indicating methylation status.
    """
    xm_tag = read.get_tag('XM')
    methylation_statuses = [char.isupper() for char in xm_tag if char in 'zZxX']
    return methylation_statuses

def clip_read_to_overlap(read, clip_start, clip_end):
    """
    Clip the read sequence and qualities to the specified start and end positions.
    """
    if read.query_sequence is None:
        logging.debug(f"query_sequence is None for read {read.query_name}")
        return read

    if read.query_qualities is None:
        logging.debug(f"query_qualities is None for read {read.query_name}")

    # Clip the sequences
    new_seq = read.query_sequence[clip_start:clip_end]
    new_qual = read.query_qualities[clip_start:clip_end] if read.query_qualities else None
    new_seq_len = len(new_seq)

    # Create new CIGAR tuples for clipped read
    new_cigartuples = []
    if clip_start > 0:
        new_cigartuples.append((4, clip_start))  # Hard clip at start
    if new_seq_len > 0:
        new_cigartuples.append((0, new_seq_len))  # Match/mismatch for the remaining sequence
    if clip_end < len(read.query_sequence):
        new_cigartuples.append((4, len(read.query_sequence) - clip_end))  # Hard clip at end

    # Update read sequence, qualities, and CIGAR
    read.query_sequence = new_seq
    read.query_qualities = new_qual
    read.cigartuples = new_cigartuples

    # Clip the XM tag if it exists
    if read.has_tag('XM'):
        xm_tag = read.get_tag('XM')
        new_xm_tag = xm_tag[clip_start:clip_end]
        read.set_tag('XM', new_xm_tag)

    return read

def parse_cigar(cigar):
    """
    Parse CIGAR string into operations and their lengths.
    """
    operations = []
    num = ""
    for char in cigar:
        if char.isdigit():
            num += char
        else:
            operations.append((char, int(num)))
            num = ""
    return operations

def update_cigar(cigar, diff):
    """
    Update CIGAR string based on alignment differences.
    """
    parsed_cigar = parse_cigar(cigar)
    updated_cigar = []

    for op, length in parsed_cigar:
        if op == 'M':
            updated_cigar.append((op, length - diff))
        elif op in ('I', 'D'):
            continue
        else:
            updated_cigar.append((op, length))

    return ''.join(f'{length}{op}' for op, length in updated_cigar if length > 0)

def apply_cigar_to_sequences(read1, read2):
    """
    Align two sequences based on their CIGAR strings.
    Remove insertions and add deletions to match both sequences.
    """
    seq1, qual1, xm1, cigar1 = read1.query_sequence, read1.query_qualities, read1.get_tag('XM'), read1.cigarstring
    seq2, qual2, xm2, cigar2 = read2.query_sequence, read2.query_qualities, read2.get_tag('XM'), read2.cigarstring

    parsed_cigar1 = parse_cigar(cigar1)
    parsed_cigar2 = parse_cigar(cigar2)

    aligned_seq1, aligned_qual1, aligned_xm1 = [], [], []
    aligned_seq2, aligned_qual2, aligned_xm2 = [], [], []

    i, j = 0, 0
    idx1, idx2 = 0, 0

    while i < len(parsed_cigar1) and j < len(parsed_cigar2):
        op1, len1 = parsed_cigar1[i]
        op2, len2 = parsed_cigar2[j]

        if op1 == 'M' and op2 == 'M':
            length = min(len1, len2)
            aligned_seq1.extend(seq1[idx1:idx1+length])
            aligned_qual1.extend(qual1[idx1:idx1+length])
            aligned_xm1.extend(xm1[idx1:idx1+length])
            aligned_seq2.extend(seq2[idx2:idx2+length])
            aligned_qual2.extend(qual2[idx2:idx2+length])
            aligned_xm2.extend(xm2[idx2:idx2+length])
            idx1 += length
            idx2 += length
            parsed_cigar1[i] = (op1, len1 - length)
            parsed_cigar2[j] = (op2, len2 - length)
            if parsed_cigar1[i][1] == 0:
                i += 1
            if parsed_cigar2[j][1] == 0:
                j += 1
        elif op1 == 'I' and op2 != 'I':
            idx1 += len1
            i += 1
        elif op2 == 'I' and op1 != 'I':
            idx2 += len2
            j += 1
        elif op1 == 'D' and op2 != 'D':
            aligned_seq1.extend(['N'] * len1)
            aligned_qual1.extend([0] * len1)
            aligned_xm1.extend(['.'] * len1)
            aligned_seq2.extend(seq2[idx2:idx2+len1])
            aligned_qual2.extend(qual2[idx2:idx2+len1])
            aligned_xm2.extend(xm2[idx2:idx2+len1])
            idx2 += len1
            i += 1
        elif op2 == 'D' and op1 != 'D':
            aligned_seq2.extend(['N'] * len2)
            aligned_qual2.extend([0] * len2)
            aligned_xm2.extend(['.'] * len2)
            aligned_seq1.extend(seq1[idx1:idx1+len2])
            aligned_qual1.extend(qual1[idx1:idx1+len2])
            aligned_xm1.extend(xm1[idx1:idx1+len2])
            idx1 += len2
            j += 1
        else:
            if op1 == 'M':
                aligned_seq1.extend(seq1[idx1:idx1+len1])
                aligned_qual1.extend(qual1[idx1:idx1+len1])
                aligned_xm1.extend(xm1[idx1:idx1+len1])
                idx1 += len1
            elif op1 == 'I':
                idx1 += len1
            elif op1 == 'D':
                aligned_seq1.extend(['N'] * len1)
                aligned_qual1.extend([0] * len1)
                aligned_xm1.extend(['.'] * len1)
            if op2 == 'M':
                aligned_seq2.extend(seq2[idx2:idx2+len2])
                aligned_qual2.extend(qual2[idx2:idx2+len2])
                aligned_xm2.extend(xm2[idx2:idx2+len2])
                idx2 += len2
            elif op2 == 'I':
                idx2 += len2
            elif op2 == 'D':
                aligned_seq2.extend(['N'] * len2)
                aligned_qual2.extend([0] * len2)
                aligned_xm2.extend(['.'] * len2)
            i += 1
            j += 1

    read1.query_sequence = ''.join(aligned_seq1)
    read1.query_qualities = aligned_qual1
    read1.set_tag('XM', ''.join(aligned_xm1))
    read2.query_sequence = ''.join(aligned_seq2)
    read2.query_qualities = aligned_qual2
    read2.set_tag('XM', ''.join(aligned_xm2))

    read1.cigarstring = update_cigar(cigar1, len(aligned_seq1) - len(seq1))
    read2.cigarstring = update_cigar(cigar2, len(aligned_seq2) - len(seq2))

    return read1, read2

def mask_discrepant_cpg_sites(seq1, seq2, xm1, xm2, start1, end1, start2, end2):
    """
    Mask CpG sites with discrepancies in methylation status between paired reads.
    """
    if seq1 is None or seq2 is None or xm1 is None or xm2 is None:
        logging.debug("One of the sequences or XM tags is None.")
        return None, None, None, None, None, None, None, None

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
        return None, None, None, None, None, None, None, None

    # Initialize masked sequences and XM tags
    masked_seq1 = list(seq1)
    masked_seq2 = list(seq2)
    masked_xm1 = list(xm1)
    masked_xm2 = list(xm2)

    # Iterate over sequences and mask discrepancies at CpG sites
    len_seq = min(len(seq1), len(seq2))
    for i in range(len_seq - 1):
        # Identify CpG sites and their possible converted forms
        seq1_cpg = seq1[i:i+2] in ['CG', 'TG']
        seq2_cpg = seq2[i:i+2] in ['CG', 'CA']

        if seq1_cpg and seq2_cpg:
            # Check methylation status in XM tags
            if (xm1[i] == 'Z' and xm2[i+1] == 'Z') or (xm1[i] == 'z' and xm2[i+1] == 'z'):
                continue  # Both reads agree on methylation status, do nothing
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
    """
    Process single-end reads to handle CpG masking and clipping.

    Args:
        reads (list): List of single-end reads.
        base_difference_threshold (int): Maximum base pair difference allowed for overlapping reads.

    Returns:
        list: List of processed reads with masked CpG sites.
    """
    result_reads = []  # List to store result reads
    reads_by_umi = defaultdict(list)  # Dictionary to group reads by UMI (Unique Molecular Identifier)
    processed_reads = set()  # Set to keep track of processed reads

    # Group reads by UMI
    for read in reads:
        umi = read.query_name.split(':')[-1]  # Extract UMI from the read name
        reads_by_umi[umi].append(read)  # Append read to the corresponding UMI group

    # Process reads within each UMI group
    for umi, umi_reads in reads_by_umi.items():
        umi_reads.sort(key=lambda r: (r.reference_name, r.reference_start))  # Sort reads by chromosome and coordinates

        for i, read1 in enumerate(umi_reads):
            if read1.query_name in processed_reads:  # Skip already processed reads
                continue

            for j in range(i + 1, len(umi_reads)):
                read2 = umi_reads[j]
                if read2.query_name in processed_reads:  # Skip already processed reads
                    continue

                if read1.reference_name == read2.reference_name:  # Ensure reads are on the same chromosome
                    start1, end1 = read1.reference_start, read1.reference_end
                    start2, end2 = read2.reference_start, read2.reference_end

                    start_diff = abs(start1 - start2)
                    end_diff = abs(end1 - end2)

                    if start_diff <= base_difference_threshold and end_diff <= base_difference_threshold:  # Check if starts and ends are within the base difference threshold
                        logging.debug(f"\nread1: {read1.query_name}, read2: {read2.query_name}")
                        logging.debug(f"\nread1: {read1.query_sequence}\nread2: {read2.query_sequence}")
                        logging.debug(f"\nread1: {read1.get_tag('XM')}\nread2: {read2.get_tag('XM')}")

                        # Check if reads have complementary XG tags
                        if (read1.get_tag('XG') == 'CT' and read2.get_tag('XG') == 'GA') or (read1.get_tag('XG') == 'GA' and read2.get_tag('XG') == 'CT'):
                            if read1.cigarstring != read2.cigarstring:
                                # Align reads using CIGAR strings
                                read1, read2 = apply_cigar_to_sequences(read1, read2)

                            masked_seq1, masked_xm1, masked_seq2, masked_xm2, clip_start1, clip_end1, clip_start2, clip_end2 = mask_discrepant_cpg_sites(
                                read1.query_sequence, read2.query_sequence, 
                                read1.get_tag('XM'), read2.get_tag('XM'), 
                                read1.reference_start, read1.reference_end,
                                read2.reference_start, read2.reference_end)

                            if masked_seq1 and masked_xm1 and masked_seq2 and masked_xm2:  # If masking was successful
                                # Update and clip the first read
                                read1.query_sequence = masked_seq1
                                read1.set_tag('XM', masked_xm1)
                                read1 = clip_read_to_overlap(read1, clip_start1, clip_end1)
                                
                                # Update and clip the second read
                                read2.query_sequence = masked_seq2
                                read2.set_tag('XM', masked_xm2)
                                read2 = clip_read_to_overlap(read2, clip_start2, clip_end2)

                                logging.debug(f"\nmask1: {read1.query_sequence}\nmask2: {read2.query_sequence}")
                                logging.debug(f"\nmask1: {read1.get_tag('XM')}\nmask2: {read2.get_tag('XM')}")

                                # Add processed reads to the result list and mark them as processed
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
    """
    Process paired-end reads to handle CpG masking and clipping.

    Args:
        reads (list of tuples): List of paired-end reads (tuples of read1 and read2).
        base_difference_threshold (int): Maximum base pair difference allowed for overlapping reads.

    Returns:
        list: List of processed paired-end reads with masked CpG sites.
    """
    result_reads = []  # List to store result reads
    paired_reads_by_umi = defaultdict(list)  # Dictionary to group paired-end reads by UMI
    processed_pairs = set()  # Set to keep track of processed pairs

    # Group paired-end reads by UMI
    for pair in reads:
        read1, read2 = pair
        umi = read1.query_name.split(':')[-1]  # Extract UMI from read1 name
        paired_reads_by_umi[umi].append(pair)  # Append read pair to the corresponding UMI group

    # Process paired-end reads within each UMI group
    for umi, umi_pairs in paired_reads_by_umi.items():
        umi_pairs.sort(key=lambda pair: (pair[0].reference_name, pair[0].reference_start))  # Sort pairs by chromosome and coordinates

        for i, pair1 in enumerate(umi_pairs):
            read1a, read1b = pair1
            if read1a.query_name in processed_pairs or read1b.query_name in processed_pairs:  # Skip already processed pairs
                continue

            for j in range(i + 1, len(umi_pairs)):
                pair2 = umi_pairs[j]
                read2a, read2b = pair2
                if read2a.query_name in processed_pairs or read2b.query_name in processed_pairs:  # Skip already processed pairs
                    continue

                if read1a.reference_name == read2a.reference_name:  # Ensure reads are on the same chromosome
                    # Identify the most upstream and downstream reads for both pairs
                    if read1a.reference_start <= read1b.reference_start:
                        most_upstream1, most_downstream1 = read1a, read1b
                    else:
                        most_upstream1, most_downstream1 = read1b, read1a

                    if read2a.reference_start <= read2b.reference_start:
                        most_upstream2, most_downstream2 = read2a, read2b
                    else:
                        most_upstream2, most_downstream2 = read2b, read2a

                    start_diff = abs(most_upstream1.reference_start - most_upstream2.reference_start)
                    end_diff = abs(most_downstream1.reference_end - most_downstream2.reference_end)

                    if start_diff <= base_difference_threshold and end_diff <= base_difference_threshold:  # Check if starts and ends are within the base difference threshold
                        logging.debug(f"\nread1a: {read1a.query_name}, read1b: {read1b.query_name}")
                        logging.debug(f"\nread1a: {read1a.query_sequence}\nread1b: {read1b.query_sequence}")
                        logging.debug(f"\nread1a: {read1a.get_tag('XM')}\nread1b: {read1b.get_tag('XM')}")
                        logging.debug(f"\nread2a: {read2a.query_name}, read2b: {read2b.query_name}")
                        logging.debug(f"\nread2a: {read2a.query_sequence}\nread2b: {read2b.query_sequence}")
                        logging.debug(f"\nread2a: {read2a.get_tag('XM')}\nread2b: {read2b.get_tag('XM')}")

                        # Check if reads have complementary XG tags
                        if (read1a.get_tag('XG') == 'CT' and read2a.get_tag('XG') == 'GA') or (read1a.get_tag('XG') == 'GA' and read2a.get_tag('XG') == 'CT'):
                            if read1a.cigarstring != read2a.cigarstring:
                                # Align reads using CIGAR strings
                                read1a, read2a = apply_cigar_to_sequences(read1a, read2a)
                            # Mask CpG sites for the first pair of reads
                            masked_seq1a, masked_xm1a, masked_seq2a, masked_xm2a, clip_start1a, clip_end1a, clip_start2a, clip_end2a = mask_discrepant_cpg_sites(
                                read1a.query_sequence, read2a.query_sequence, 
                                read1a.get_tag('XM'), read2a.get_tag('XM'), 
                                read1a.reference_start, read1a.reference_end,
                                read2a.reference_start, read2a.reference_end)

                            # Mask CpG sites for the second pair of reads
                            if read1b.cigarstring != read2b.cigarstring:
                                # Align reads using CIGAR strings
                                read1b, read2b = apply_cigar_to_sequences(read1b, read2b)
                            masked_seq1b, masked_xm1b, masked_seq2b, masked_xm2b, clip_start1b, clip_end1b, clip_start2b, clip_end2b = mask_discrepant_cpg_sites(
                                read1b.query_sequence, read2b.query_sequence, 
                                read1b.get_tag('XM'), read2b.get_tag('XM'), 
                                read1b.reference_start, read1b.reference_end,
                                read2b.reference_start, read2b.reference_end)

                            if masked_seq1a and masked_xm1a and masked_seq2a and masked_xm2a and masked_seq1b and masked_xm1b and masked_seq2b and masked_xm2b:  # If masking was successful
                                # Update and clip the first read in the first pair
                                read1a.query_sequence = masked_seq1a
                                read1a.set_tag('XM', masked_xm1a)
                                read1a = clip_read_to_overlap(read1a, clip_start1a, clip_end1a)

                                # Update and clip the second read in the first pair
                                read2a.query_sequence = masked_seq2a
                                read2a.set_tag('XM', masked_xm2a)
                                read2a = clip_read_to_overlap(read2a, clip_start2a, clip_end2a)

                                # Update and clip the first read in the second pair
                                read1b.query_sequence = masked_seq1b
                                read1b.set_tag('XM', masked_xm1b)
                                read1b = clip_read_to_overlap(read1b, clip_start1b, clip_end1b)

                                # Update and clip the second read in the second pair
                                read2b.query_sequence = masked_seq2b
                                read2b.set_tag('XM', masked_xm2b)
                                read2b = clip_read_to_overlap(read2b, clip_start2b, clip_end2b)

                                logging.debug(f"\nmask1a: {read1a.query_sequence}\nmask2a: {read2a.query_sequence}")
                                logging.debug(f"\nmask1a: {read1a.get_tag('XM')}\nmask2a: {read2a.get_tag('XM')}")
                                logging.debug(f"\nmask1b: {read1b.query_sequence}\nmask2b: {read2b.query_sequence}")
                                logging.debug(f"\nmask1b: {read1b.get_tag('XM')}\nmask2b: {read2b.get_tag('XM')}")

                                # Add processed reads to the result list and mark them as processed
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
    """
    Fetch reads from a BAM file within the specified region.

    Args:
        bam_file (pysam.AlignmentFile): BAM file to fetch reads from.
        contig (str): Contig (chromosome) name.
        start (int): Start position.
        end (int): End position.

    Returns:
        tuple: List of single-end reads and list of paired-end reads.
    """
    single_end_reads = []  # List to store single-end reads
    paired_end_reads = defaultdict(list)  # Dictionary to store paired-end reads grouped by read name

    # Iterate over reads in the specified region
    for read in bam_file.fetch(contig, start, end):
        if read.is_unmapped:
            continue  # Skip unmapped reads

        if read.is_secondary or read.is_supplementary:
            continue  # Skip secondary and supplementary alignments

        if read.is_paired:
            paired_end_reads[read.query_name].append(read)  # Append to paired-end reads
        else:
            single_end_reads.append(read)  # Append to single-end reads

    # Convert dictionary of paired-end reads to list of tuples (read1, read2)
    paired_end_reads_list = [reads for reads in paired_end_reads.values() if len(reads) == 2]
    return single_end_reads, paired_end_reads_list

def process_bam_segment(args):
    """
    Process a segment of a BAM file.

    Args:
        args (tuple): Arguments for processing the segment, including BAM file path, output path,
                      contig, start and end positions, progress dictionary, and process ID.
    """
    (bam_path, phase_path, contig, start, end, progress_dict, process_id) = args

    logging.info(f"Start processing {contig} from {start} to {end} in process {process_id}")

    # Open BAM file for reading
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    header = bam_file.header.copy()

    # Open output BAM file for writing
    out_phase_bam = pysam.AlignmentFile(phase_path, "wb", header=header)

    # Fetch reads from the specified region
    single_end_reads, paired_end_reads = fetch_reads(bam_file, contig, start, end)

    # Process single-end and paired-end reads
    processed_single_end_reads = process_single_end_reads(single_end_reads)
    processed_paired_end_reads = process_paired_end_reads(paired_end_reads)

    # Write processed reads to output BAM file
    for read in processed_single_end_reads:
        out_phase_bam.write(read)

    for read in processed_paired_end_reads:
        out_phase_bam.write(read)

    out_phase_bam.close()  # Close output BAM file
    bam_file.close()  # Close input BAM file
    logging.info(f"Finished processing {contig} from {start} to {end} in process {process_id}.")
    progress_dict[process_id] = 100  # Update progress dictionary

def calculate_segments_by_length(bam_path, num_threads, target_reads_per_segment, chromosomes):
    """
    Calculate segments by length for parallel processing.

    Args:
        bam_path (str): Path to the BAM file.
        num_threads (int): Number of threads to use.
        target_reads_per_segment (int): Target number of reads per segment.
        chromosomes (list): List of chromosomes to process.

    Returns:
        dict: Dictionary with contigs as keys and number of segments as values.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        total_length = sum(
            l for r, l in zip(bam.references, bam.lengths) 
            if chromosomes is None or r in chromosomes)  # Calculate total length of specified chromosomes
        segments = {}
        for contig, length in zip(bam.references, bam.lengths):
            if chromosomes is not None and contig not in chromosomes:
                continue
            contig_reads = bam.count(reference=contig)
            if contig_reads == 0:
                logging.info(f"Skipping {contig} due to zero reads.")
                continue

            proportion = length / total_length  # Calculate proportion of contig length to total length
            num_segments = max(num_threads, int(ceil(proportion * target_reads_per_segment)))  # Calculate number of segments
            segments[contig] = num_segments
    return segments

def calculate_segments_by_reads(bam_path, num_threads, target_reads_per_segment, chromosomes):
    """
    Calculate segments by read count for parallel processing.

    Args:
        bam_path (str): Path to the BAM file.
        num_threads (int): Number of threads to use.
        target_reads_per_segment (int): Target number of reads per segment.
        chromosomes (list): List of chromosomes to process.

    Returns:
        list: List of segments with contig, start, and end positions.
    """
    segments = calculate_segments_by_length(
        bam_path, num_threads, target_reads_per_segment, chromosomes)  # Calculate initial segments by length
    bam = pysam.AlignmentFile(bam_path, "rb")
    adjusted_segments = []

    # Adjust segments based on read positions
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
    """
    Split BAM processing into segments for multicore processing.

    Args:
        bam_path (str): Path to the BAM file.
        num_threads (int): Number of threads to use.
        target_reads_per_segment (int): Target number of reads per segment.
        chromosome_group (list): List of chromosomes to process.
        merged_duplex_path (str): Path to the output merged BAM file.

    Returns:
        None
    """
    logging.info(f"Segmenting BAM reads into contig segments for multicore processing using {num_threads} threads")
    segments = calculate_segments_by_reads(bam_path, num_threads, target_reads_per_segment, chromosome_group)  # Calculate segments by read count
    jobs = []
    manager = multiprocessing.Manager()
    progress_dict = manager.dict()
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        header = bam.header  # Get BAM file header
        job_id = 0

        for segment in segments:
            contig, start, end = segment 
            phase_path = f"{merged_duplex_path.replace('.bam', '')}_duplex_temp{job_id}.bam"
            job = (bam_path, phase_path, contig, start, end, progress_dict, job_id)
            logging.info(f"Process {job_id}. Segment of {contig}: {start}-{end}")
            jobs.append(job)
            job_id += 1

        # Process segments in parallel using multiprocessing
        with multiprocessing.Pool(processes=num_threads) as pool:
            pool.map(process_bam_segment, jobs)

        # Merge temporary BAM files into the final output BAM file
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
    """
    Get the list of contigs (chromosomes) from the BAM file.

    Args:
        bam_path (str): Path to the BAM file.

    Returns:
        list: List of contig names.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        contigs = bam_file.references  # Get contig names
    return list(contigs)

def get_reads_per_chromosome(bam_path, chromosomes):
    """
    Calculate the number of reads per chromosome.

    Args:
        bam_path (str): Path to the BAM file.
        chromosomes (list): List of chromosomes to process.

    Returns:
        dict: Dictionary with chromosome names as keys and read counts as values.
    """
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
    """
    Group chromosomes into segments for parallel processing based on read count.

    Args:
        reads_per_chromosome (dict): Dictionary with chromosome names as keys and read counts as values.
        num_threads (int): Number of threads to use.
        target_reads_per_segment (int): Target number of reads per segment.

    Returns:
        list: List of chromosome groups.
    """
    sorted_chromosomes = sorted(reads_per_chromosome, key=reads_per_chromosome.get, reverse=True)  # Sort chromosomes by read count
    groups = []
    current_group = []
    current_reads = 0

    # Group chromosomes based on the target reads per segment
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
    """
    Main function to orchestrate BAM file processing.
    """
    args = parse_args()  # Parse command-line arguments
    chromosomes = args.chromosomes.split(',') if args.chromosomes else []

    if "core_human" in chromosomes:
        chromosomes = ['chr'+str(i) for i in range(1, 23)] + ['chrX', 'chrY']

    bam_files = [("CT", args.bam_path)]

    for bam_type, bam_path in bam_files:
        logging.info(f"Starting {bam_type} processing")
        
        if not chromosomes:
            chromosomes = get_contigs(bam_path)  # Get list of contigs if not provided
        
        chromosome_bam_files_duplex = []
        
        reads_per_chromosome = get_reads_per_chromosome(bam_path, chromosomes)  # Get reads per chromosome
        chromosome_groups = group_chromosomes_by_reads(reads_per_chromosome, args.num_threads, args.target_reads_per_segment)  # Group chromosomes by reads

        for group_index, group in enumerate(chromosome_groups):
            group_id = f"group{group_index}"
            merged_duplex_group_path = bam_path.replace('.bam', f'_{group_id}_duplex.bam')

            split_bam_processing(bam_path, 
                                 args.num_threads, 
                                 args.target_reads_per_segment, 
                                 group,
                                 merged_duplex_group_path)  # Split BAM processing for each group

            chromosome_bam_files_duplex.append(merged_duplex_group_path)

        merged_duplex_path = bam_path.replace('.bam', '_duplex.bam')

        # Merge all temporary BAM files into the final output BAM file
        subprocess.run(['samtools', 'merge', '-f', '-@', str(args.num_threads), merged_duplex_path] 
                       + chromosome_bam_files_duplex, check=True)

        # Remove temporary files
        for file_path in chromosome_bam_files_duplex:
            os.remove(file_path)

    logging.info("All processing completed")

if __name__ == "__main__":
    main()
