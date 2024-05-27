# Rule to calculate coverage statistics using bedtools and samtools.
rule calculate_coverage:
    input:
        # Input BAM and bedGraph files
        s_bam="{output_path}/{sample}/03_align_fastq/{sample}_s.bam",
        sd_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_clipped.bam",
        sd_bedGraph="{output_path}/{sample}/05_call_methylation/{sample}_sd_CpG.bedGraph",
    output:
        # Output file containing genome coverage statistics
        s="{output_path}/{sample}/07_calculate_statistics/{sample}_s_genomeCoverageBed.txt",
        sd="{output_path}/{sample}/07_calculate_statistics/{sample}_sd_genomeCoverageBed.txt",
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
        time_min=lambda wildcards, input, threads: get_time_min(wildcards, input, "calculate_coverage", threads),
        cpus=1,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        # Shell command to filter and process BAM file, then calculate genome coverage using bedtools.
        "samtools view -h -b -F 0x400 {input.sd_bam} | {params.intersect}bedtools "
        "genomecov -ibam - -g {params.genome_file} 1> {output.sd} "
        "2>> {log}\n\n"
        "samtools view -h -b -F 0x400 {input.s_bam} | {params.intersect}bedtools "
        "genomecov -ibam - -g {params.genome_file} 1> {output.s} "
        "2>> {log}"


# Rule to calculate conversion statistics using the genome coverage data.
rule calculate_stats:
    input:
        # Input coverage data
        s_cov="{output_path}/{sample}/07_calculate_statistics/{sample}_s_genomeCoverageBed.txt",
        sd_cov="{output_path}/{sample}/07_calculate_statistics/{sample}_sd_genomeCoverageBed.txt",
        bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_clipped.bam",
        CHH="{output_path}/{sample}/05_call_methylation/{sample}_CHH_mbias.txt",
        CHG="{output_path}/{sample}/05_call_methylation/{sample}_CHG_mbias.txt",
        CpG="{output_path}/{sample}/05_call_methylation/{sample}_CpG_mbias.txt",
        # json="{output_path}/{sample}/04_deduplicate_bam/{sample}_gencore.json",
        og_bam="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.bam",
    output:
        # Output file containing conversion and coverage statistics
        s_con="{output_path}/{sample}/07_calculate_statistics/{sample}_s_conversion_and_coverage.txt",
        sd_con="{output_path}/{sample}/07_calculate_statistics/{sample}_sd_conversion_and_coverage.txt",
        split="{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_pe.deduplicated_splitting_report.txt",
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

        transform_methylation(input.bam, output.split, D_conversion)

# Rule to compile multiple Bismark reports into a single HTML report for easier analysis and visualization.
rule bismark2report:
    input:
        # Required report files from various Bismark processing steps
        align_rep="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_PE_report.txt",
        nuc_rep="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.nucleotide_stats.txt",
        # dedup_rep="{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_pe.deduplication_report.txt",
        split_rep="{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_pe.deduplicated_splitting_report.txt",
        mbias_rep="{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_pe.deduplicated.M-bias.txt",
    output:
        "{output_path}/{sample}/07_calculate_statistics/{sample}_r1_bismark_bt2_PE_report.html",
    log:
        "{output_path}/{sample}/logs/{sample}_bismark2report.log",
    conda:
        "../envs/bismark.yaml",
    threads: 1  # Single-threaded as report generation is not CPU-intensive
    resources:
        time_min=10,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        cpus=1,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        'bismark2report --alignment_report {input.align_rep} '
        '--splitting_report {input.split_rep} --mbias_report {input.mbias_rep} '
        '--nucleotide_report {input.nuc_rep} --dir {wildcards.output_path}/{wildcards.sample}/07_calculate_statistics &>> {log}'