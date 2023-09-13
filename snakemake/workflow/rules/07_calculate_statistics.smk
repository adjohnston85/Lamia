# Rule to calculate coverage statistics using bedtools and samtools.
rule calculate_coverage:
    input:
        # Input BAM and bedGraph files
        bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_sd.bam",
        bedGraph="{output_path}/{sample}/05_call_methylation/{sample}_sd_CpG.bedGraph",
    output:
        # Output file containing genome coverage statistics
        "{output_path}/{sample}/07_calculate_statistics/{sample}_genomeCoverageBed.txt",
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
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        # Shell command to filter and process BAM file, then calculate genome coverage using bedtools.
        "samtools view -h -b -F 0x400 {input.bam} | {params.intersect}bedtools "
        "genomecov -ibam - -g {params.genome_file} 1> {output} "
        "2>> {log}"


# Rule to calculate conversion statistics using the genome coverage data.
rule calculate_conversion:
    input:
        # Input coverage data
        cov="{output_path}/{sample}/07_calculate_statistics/{sample}_genomeCoverageBed.txt",
    output:
        # Output file containing conversion and coverage statistics
        con="{output_path}/{sample}/07_calculate_statistics/{sample}_sd_conversion_and_coverage.txt",
    threads: 1
    resources:
        time_min=30,
        cpus=1,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    run:
        # Read and process genome coverage data.
        covFile = pd.read_table(input.cov, header=None, names=['chr','depth','base_count','chr_size_bp','fraction'])
        genomeCov = covFile[covFile['chr'] == 'genome']
        averageCov = sum(genomeCov['depth'] * genomeCov['base_count'])/genomeCov.loc[genomeCov.index[0],'chr_size_bp']
        
        # Write average coverage to the output file.
        pd.DataFrame(data={'sample_name':wildcards.sample,
                    'Genome':D_sample_details[wildcards.sample]["genome"], 'Average_Coverage':averageCov},
                     index=['coveragedetails']).to_csv(path_or_buf=output.con, sep = '\t')

        # Calculate conversion statistics for different contexts.
        L_contexts = ['CHH','CHG','CpG']
        with open(output.con, 'a') as F_coverage:
            F_coverage.write("\n")
            for context in L_contexts:
                L_stats = conversion_estimator(wildcards.sample, context)
                if context == 'CHH':
                    F_coverage.write('\t' + '\t'.join(L_stats[1]) + '\n')
                F_coverage.write(context + '\t' + '\t'.join(L_stats[0]) + '\n')