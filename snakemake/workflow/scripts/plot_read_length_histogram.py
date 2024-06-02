import sys
import pysam
import matplotlib.pyplot as plt
import numpy as np

def plot_histogram(bam_file, histogram_file, length_counts_file):
    paired_end_lengths = []
    single_end_lengths = []
    
    # Dictionaries to store length counts
    single_end_counts = {}
    paired_end_counts = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_paired and read.is_read1:
                # Calculate insert size for paired-end reads (only for first read in pair to avoid duplicates)
                insert_size = abs(read.template_length)
                paired_end_lengths.append(insert_size)
                if insert_size in paired_end_counts:
                    paired_end_counts[insert_size] += 1
                else:
                    paired_end_counts[insert_size] = 1
            elif not read.is_paired:
                # Use read length for single-end reads
                read_length = read.query_length
                single_end_lengths.append(read_length)
                if read_length in single_end_counts:
                    single_end_counts[read_length] += 1
                else:
                    single_end_counts[read_length] = 1
    
    # Combine length lists and calculate the 99th percentile
    all_lengths = paired_end_lengths + single_end_lengths
    length_99th_percentile = np.percentile(all_lengths, 99.9) if all_lengths else 0

    # Determine the range for the bins
    max_length = max(paired_end_lengths + single_end_lengths) if paired_end_lengths or single_end_lengths else 0
    
    # Create the bins
    bins = np.arange(1, max_length + 2) - 0.5

    # Plot the histogram
    plt.figure(figsize=(16, 9))
    n, bins, patches = plt.hist([paired_end_lengths, single_end_lengths], bins=bins, stacked=True, color=["blue", "orange"], label=["Paired-End", "Single-End"])
    plt.title("Read Length Distribution")
    plt.xlabel("Read Length (bp)")
    plt.ylabel("Frequency")
    plt.legend()
    plt.grid(axis="y", alpha=0.75)
    
    # Set x-axis limit
    plt.xlim(0, length_99th_percentile)

    # Label the peaks
    peak_paired_idx = np.argmax(n[0])
    peak_single_idx = np.argmax(n[1])
    peak_paired = bins[peak_paired_idx]
    peak_single = bins[peak_single_idx]
    plt.text(peak_paired, n[0][peak_paired_idx] + 100, "Peak Paired-End: {}".format(int(peak_paired + 0.5)), color="blue", ha="center")
    plt.text(peak_single, n[1][peak_single_idx] + 100, "Peak Single-End: {}".format(int(peak_single + 0.5)), color="orange", ha="center")

    plt.savefig(histogram_file)
    plt.close()
    
    # Output length counts to a file
    with open(length_counts_file, "w") as f:
        f.write("Length\tSingle-End Count\tPaired-End Count\n")
        all_lengths = sorted(set(single_end_counts.keys()).union(paired_end_counts.keys()))
        for length in all_lengths:
            se_count = single_end_counts.get(length, 0)
            pe_count = paired_end_counts.get(length, 0)
            f.write("{}\t{}\t{}\n".format(length, se_count, pe_count))

if __name__ == "__main__":
    bam_file = sys.argv[1]
    histogram_file = sys.argv[2]
    length_counts_file = sys.argv[3]
    plot_histogram(bam_file, histogram_file, length_counts_file)
