# Rule deduplicates aligned and paired BAM files using the Bismark deduplication tool.
rule bismark_deduplicate:
    input:
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.bam",
    output:
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.deduplicated.bam",
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.deduplication_report.txt",
    # Log files capture standard output and standard error..
    log:
        out="{output_path}/{sample}/logs/{sample}_bismark_deduplicate.stdout",
        err="{output_path}/{sample}/logs/{sample}_bismark_deduplicate.stderr",
    conda:
        "../envs/bismark.yaml",
    # Dynamic CPU allocation with a maximum limit of 64.
    threads: lambda wcs: get_cpus(1,64),
    # Resource constraints such as max time in minutes, CPU, memory, and other user/account specific details.
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "bismark_deduplicate", threads),
        mem_mb=get_mem_mb,
        cpus=lambda wcs, threads: threads,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        'deduplicate_bismark --bam --paired {input} --output_dir '
        '{wildcards.output_path}/{wildcards.sample}/03_align_fastq >> {log.out} 2>> {log.err}'


# Rule to extract methylation data from deduplicated BAM files using Bismark, generating various output files.
rule bismark_methylation:
    input:
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.deduplicated.bam",
    output:
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.deduplicated.bedGraph.gz",
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.deduplicated_splitting_report.txt",
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.deduplicated.M-bias.txt",
    log:
        out="{output_path}/{sample}/logs/{sample}_bismark_methylation.stdout",
        err="{output_path}/{sample}/logs/{sample}_bismark_methylation.stderr",
    params:
        # '--ignore_r2 2' is applied for 'bs-seq' library type based on Bismark documentation.
        ignore='--ignore_r2 2 ' if (lambda wcs: D_sample_details[wcs.sample]['library_type']) == 'bs-seq' else '',
        # Multicore parallelism is set to 1/3 of the allocated cores due to high resource consumption.
        parallel=lambda wcs, threads: get_parallel(wcs, 3, threads),
    conda:
        "../envs/bismark.yaml",
    threads:
        lambda wcs: get_cpus(1,64),
    resources:
        time_min=lambda wildcards, input, threads: get_time_min(wildcards, input, "bismark_methylation", threads),
        mem_mb=get_mem_mb,
        cpus=lambda wcs, threads: threads,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        # Runs bismark_methylation_extractor with specified options and parameters
        'bismark_methylation_extractor {input} --bedGraph {params.ignore} '
        '--output {wildcards.output_path}/{wildcards.sample}/03_align_fastq --multicore {params.parallel} '
        '--scaffolds --gzip >> {log.out} 2>> {log.err} \n\n'
        # Deletes large, unnecessary intermediate files to save storage space.
        'rm -f {wildcards.output_path}/{wildcards.sample}/03_align_fastq/C??_O?_{wildcards.sample}_r1_bismark_bt2_pe.deduplicated.txt*'


# Rule to compile multiple Bismark reports into a single HTML report for easier analysis and visualization.
rule bismark2report:
    input:
        # Required report files from various Bismark processing steps
        align_rep="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_PE_report.txt",
        dedup_rep="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.deduplication_report.txt",
        split_rep="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.deduplicated_splitting_report.txt",
        mbias_rep="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.deduplicated.M-bias.txt",
        nuc_rep="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.nucleotide_stats.txt",
    output:
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_PE_report.html",
    log:
        "{output_path}/{sample}/logs/{sample}_bismark2report.log",
    conda:
        "../envs/bismark.yaml",
    threads: 1  # Single-threaded as report generation is not CPU-intensive
    resources:
        time_min=10,
        mem_mb=get_mem_mb,
        cpus=1,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        'bismark2report --alignment_report {input.align_rep} --dedup_report {input.dedup_rep} '
        '--splitting_report {input.split_rep} --mbias_report {input.mbias_rep} '
        '--nucleotide_report {input.nuc_rep} --dir {wildcards.output_path}/{wildcards.sample}/03_align_fastq &>> {log}'

