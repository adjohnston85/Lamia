rule plot_read_length_histograms:
    input:
        s_bam="{output_path}/{sample}/03_align_fastq/{sample}_merged_and_sorted.bam",  # Input sorted BAM file
        sd_bam=lambda wildcards: (
            "{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe_duplex.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
            if "umi_len" in D_sample_details[wildcards.sample]
            else "{output_path}/{sample}/04_deduplicate_bam/{sample}_dedup_merged_and_sorted_se_pe.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
        ),
    output:
        s_histogram="{output_path}/{sample}/07_calculate_statistics/{sample}_read_length_histogram_s.png",
        sd_histogram="{output_path}/{sample}/07_calculate_statistics/{sample}_read_length_histogram_sd.png",
    params:
        sdd_bam=lambda wildcards: (
            "{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe_duplex.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
            if "umi_len" in D_sample_details[wildcards.sample]
            else "None"
        ),
    log:
        "{output_path}/{sample}/logs/{sample}_plot_read_length_histogram.log"
    conda:
        "../envs/gencore.yaml"
    threads: 1
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 1, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        """
        python scripts/plot_read_length_histogram.py {input.s_bam} {output.s_histogram} &> {log}

        python scripts/plot_read_length_histogram.py {input.sd_bam} {output.sd_histogram} &>> {log}

        # Conditionally process duplex BAM
        if [ "{params.sdd_bam}" != "None" ]; then
            python scripts/plot_read_length_histogram.py {params.sdd_bam} "{wildcards.output_path}/{wildcards.sample}/07_calculate_statistics/{wildcards.sample}_read_length_histogram_sdd.png" &>> {log}   
        fi    
        """

rule calculate_coverage:
    input:
        # Input BAM and bedGraph files
        s_bam="{output_path}/{sample}/03_align_fastq/{sample}_merged_and_sorted.bam",  # Input sorted BAM file
        sd_bedGraph="{output_path}/{sample}/05_call_methylation/{sample}_sd_CpG.bedGraph",
        sd_bam=lambda wildcards: (
            "{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe_duplex.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
            if "umi_len" in D_sample_details[wildcards.sample]
            else "{output_path}/{sample}/04_deduplicate_bam/{sample}_dedup_merged_and_sorted_se_pe.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
        )
    output:
        # Output file containing genome coverage statistics
        s="{output_path}/{sample}/07_calculate_statistics/{sample}_s_genomeCoverageBed.txt",
        sd="{output_path}/{sample}/07_calculate_statistics/{sample}_sd_genomeCoverageBed.txt",
        s_roi="{output_path}/{sample}/07_calculate_statistics/{sample}_s_genomeCoverageBed_roi.txt",
        sd_roi="{output_path}/{sample}/07_calculate_statistics/{sample}_sd_genomeCoverageBed_roi.txt",
    params:
        # Define command to intersect with ROI bed if available
        intersect=lambda wcs: 'bedtools intersect -abam stdin -b ' +
                  D_sample_details[wcs.sample]["roi_bed"] + ' | '
                  if "roi_bed" in D_sample_details[wcs.sample] else '',
        # Genome file path
        genome_file=lambda wcs: D_sample_details[wcs.sample]["genome_path"] +
                    '/' + D_sample_details[wcs.sample]["genome"] + '.genome',
    log:
        "{output_path}/{sample}/logs/{sample}_calculate_coverage.log",
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    resources:
        time_min=lambda wildcards, input, threads: get_time_min(wildcards, input, 100, threads),
        cpus=1,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        # Shell command to filter and process BAM file, then calculate genome coverage using bedtools.
        "samtools view -h -b -F 0x400 {input.sd_bam} | bedtools genomecov -ibam - -g {params.genome_file} 1> {output.sd} 2>> {log}\n\n"
        "samtools view -h -b -F 0x400 {input.s_bam} | bedtools genomecov -ibam - -g {params.genome_file} 1> {output.s} 2>> {log}\n\n"
        "samtools view -h -b -F 0x400 {input.sd_bam} | {params.intersect}bedtools genomecov -ibam - -g {params.genome_file} 1> {output.sd_roi} 2>> {log}\n\n"
        "samtools view -h -b -F 0x400 {input.s_bam} | {params.intersect}bedtools genomecov -ibam - -g {params.genome_file} 1> {output.s_roi} 2>> {log}"

rule plot_coverage:
    input:
        s="{output_path}/{sample}/07_calculate_statistics/{sample}_s_genomeCoverageBed.txt",
        sd="{output_path}/{sample}/07_calculate_statistics/{sample}_sd_genomeCoverageBed.txt",
        s_roi="{output_path}/{sample}/07_calculate_statistics/{sample}_s_genomeCoverageBed_roi.txt",
        sd_roi="{output_path}/{sample}/07_calculate_statistics/{sample}_sd_genomeCoverageBed_roi.txt"
    output:
        plot="{output_path}/{sample}/07_calculate_statistics/{sample}_coverage_plot.png",
        plot_roi="{output_path}/{sample}/07_calculate_statistics/{sample}_coverage_plot_roi.png"
    log:
        "{output_path}/{sample}/logs/{sample}_plot_coverage.log",
    conda:
        "../envs/matplotlib.yaml"
    threads: 1
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "plot_coverage", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        """
        python scripts/plot_coverage.py {input.s} {input.sd} {output.plot} "Coverage Plot - {wildcards.sample}" &> {log}
        python scripts/plot_coverage.py {input.s_roi} {input.sd_roi} {output.plot_roi} "Coverage Plot ROI - {wildcards.sample}" &> {log}
        """

# Rule to calculate conversion statistics using the genome coverage data.
rule calculate_stats:
    input:
        # Input coverage data
        s_cov="{output_path}/{sample}/07_calculate_statistics/{sample}_s_genomeCoverageBed.txt",
        sd_cov="{output_path}/{sample}/07_calculate_statistics/{sample}_sd_genomeCoverageBed.txt",
        CHH="{output_path}/{sample}/05_call_methylation/{sample}_CHH_mbias.txt",
        CHG="{output_path}/{sample}/05_call_methylation/{sample}_CHG_mbias.txt",
        CpG="{output_path}/{sample}/05_call_methylation/{sample}_CpG_mbias.txt",
        bam=lambda wildcards: (
            "{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe_duplex.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
            if "umi_len" in D_sample_details[wildcards.sample]
            else "{output_path}/{sample}/04_deduplicate_bam/{sample}_dedup_merged_and_sorted_se_pe.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
        ),
    output:
        # Output file containing conversion and coverage statistics
        s_con="{output_path}/{sample}/07_calculate_statistics/{sample}_s_conversion_and_coverage.txt",
        sd_con="{output_path}/{sample}/07_calculate_statistics/{sample}_sd_conversion_and_coverage.txt",
        # split="{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_pe.deduplicated_splitting_report.txt",
        mbias="{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_pe.deduplicated.M-bias.txt",
        # dedup="{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_pe.deduplication_report.txt",
    threads: 1
    resources:
        time_min=30,
        cpus=1,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    run:
        D_conversion = {"CpG":[0,0], "CHH":[0,0], "CHG":[0,0]}
        # Read and process genome coverage data for both s_cov and sd_cov.
        for coverage_file, output_file in zip([input.s_cov, input.sd_cov], [output.s_con, output.sd_con]):
            covFile = pd.read_table(coverage_file, header=None, names=['chr','depth','base_count','chr_size_bp','fraction'])
            genomeCov = covFile[covFile['chr'] == 'genome']
            averageCov = sum(genomeCov['depth'] * genomeCov['base_count'])/genomeCov.loc[genomeCov.index[0],'chr_size_bp']

            # Write average coverage to the output file.
            pd.DataFrame(data={'sample_name':wildcards.sample,
                        'Genome':D_sample_details[wildcards.sample]["genome"], 'Average_Coverage':averageCov},
                         index=['coveragedetails']).to_csv(path_or_buf=output_file, sep = '\t')

            # Calculate conversion statistics for different contexts.
            with open(output_file, 'a') as F_coverage:
                w_option = 'w'
                F_coverage.write("\n")
                for context in D_conversion:
                    suffix = "_s" if "_sd_genomeCoverageBed.txt" not in coverage_file else "_sd"
                    L_stats, D_conversion = conversion_estimator(wildcards.sample, context, D_conversion, suffix)
                    if context == "CpG":
                        F_coverage.write('\t' + '\t'.join(L_stats[1]) + '\n')
                    F_coverage.write(context + '\t' + '\t'.join(L_stats[0]) + '\n')

                    transform_mbias(
                        "{0}/{1}/05_call_methylation/{1}_{2}_mbias.txt".format(
                            wildcards.output_path, wildcards.sample, context),
                        output.mbias, context, w_option
                    )
                    w_option = 'a'

        # process_dedup_json(input.json, output.dedup, input.og_bam)

        # transform_methylation(input.bam, output.split, D_conversion)

rule gencore:
    input:
        bam_pe="{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info.bam",
        bam_se="{output_path}/{sample}/03_align_fastq/{sample}_combined_bismark_bt2.bam",
    output:
        gencore_bam="{output_path}/{sample}/07_calculate_statistics/{sample}_{genome}_gencore.bam",
    params:
        umi_prefix= lambda wcs: "--umi_prefix={} ".format(D_sample_details[wcs.sample]["umi_prefix"]) if "umi_len" in D_sample_details[wcs.sample] else "",
        b=lambda wcs: '-b ' + D_sample_details[wcs.sample]["roi_bed"] + ' '
          if "roi_bed" in D_sample_details[wcs.sample] else "",
        r=lambda wcs: D_sample_details[wcs.sample]["genome_path"] + \
          "/Bisulfite_Genome/" + D_sample_details[wcs.sample]["genome"] + \
          "." + wcs.genome + "_conversion.fa",
        mG_bam="{output_path}/{sample}/07_calculate_statistics/{sample}_{genome}_genome.bam",
    log:
        "{output_path}/{sample}/logs/{sample}_gencore_{genome}.log",
    conda:
        "../envs/gencore.yaml"
    threads: lambda wcs: get_cpus(1,64)
    log:
        "{output_path}/{sample}/logs/{sample}_gencore.log",
    conda:
        "../envs/gencore.yaml"
    threads: lambda wcs: get_cpus(1,64)
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 500, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell: 
        'samtools merge -@ {threads} -f {output.gencore_bam}.unsorted {input.bam_se} {input.bam_pe} &> {log}\n\n'

        'samtools sort -@ {threads} -O bam -o {output.gencore_bam}.unprocessed {output.gencore_bam}.unsorted &>> {log}\n\n'

        'samtools view -h {output.gencore_bam}.unprocessed | grep -e "XG:Z:{wildcards.genome}" -e "@" '
        '| samtools view -bS - 1> {params.mG_bam} '
        '2>> {log} \n\n'

        'gencore {params.umi_prefix}-s 1 -d 2 '
        '-i {params.mG_bam} -o {output.gencore_bam} -r {params.r} {params.b}'
        '-j {wildcards.output_path}/{wildcards.sample}/07_calculate_statistics/gencore_{wildcards.genome}.json '
        '-h {wildcards.output_path}/{wildcards.sample}/07_calculate_statistics/gencore_{wildcards.genome}.html '
        '&>> {log}\n\n'

        'rm {params.mG_bam} {output.gencore_bam}.unsorted {output.gencore_bam}.unprocessed'

# Rule to compile multiple Bismark reports into a single HTML report for easier analysis and visualization.
# rule bismark2report:
#     input:
#         # Required report files from various Bismark processing steps
#         align_rep="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_PE_report.txt",
#         nuc_rep="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.nucleotide_stats.txt",
#         # dedup_rep="{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_pe.deduplication_report.txt",
#         split_rep="{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_pe.deduplicated_splitting_report.txt",
#         mbias_rep="{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_pe.deduplicated.M-bias.txt",
#     output:
#         "{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_PE_report.html",
#     log:
#         "{output_path}/{sample}/logs/{sample}_bismark2report.log",
#     conda:
#         "../envs/bismark.yaml",
#     threads: 1  # Single-threaded as report generation is not CPU-intensive
#     resources:
#         time_min=10,
#         mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
#         cpus=1,
#         account=lambda wcs: D_sample_details[wcs.sample]['account'],
#         email=lambda wcs: D_sample_details[wcs.sample]['email'],
#         partition=""
#     shell:
#         'bismark2report --alignment_report {input.align_rep} '
#         '--splitting_report {input.split_rep} --mbias_report {input.mbias_rep} '
#         '--nucleotide_report {input.nuc_rep} --dir {wildcards.output_path}/{wildcards.sample}/07_calculate_statistics &>> {log}'