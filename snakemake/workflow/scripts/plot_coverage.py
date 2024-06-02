import sys
import pandas as pd
import matplotlib.pyplot as plt

def read_coverage(file):
    """Reads coverage data from a genomeCoverageBed file."""
    data = pd.read_csv(file, sep='\t', header=None, names=["chrom", "pos", "coverage"])
    return data

def plot_coverage(s_file, sd_file, output_file, title):
    """Plots the coverage data with separate plots for each chromosome."""
    s_data = read_coverage(s_file)
    sd_data = read_coverage(sd_file)

    chromosomes = s_data['chrom'].unique()
    
    fig, axes = plt.subplots(len(chromosomes), 1, figsize=(20, len(chromosomes) * 2), sharex=True)
    if len(chromosomes) == 1:
        axes = [axes]

    for ax, chrom in zip(axes, chromosomes):
        s_chrom_data = s_data[s_data['chrom'] == chrom]
        sd_chrom_data = sd_data[sd_data['chrom'] == chrom]

        ax.fill_between(s_chrom_data['pos'], 0, s_chrom_data['coverage'], color='blue', step='mid', label='s', alpha=0.5)
        ax.fill_between(sd_chrom_data['pos'], 0, -sd_chrom_data['coverage'], color='orange', step='mid', label='sd', alpha=0.5)

        ax.set_title(f"{chrom}")
        ax.legend(loc='upper right')

    plt.suptitle(title)
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    s_file = sys.argv[1]
    sd_file = sys.argv[2]
    output_file = sys.argv[3]
    title = sys.argv[4]
    plot_coverage(s_file, sd_file, output_file, title)

