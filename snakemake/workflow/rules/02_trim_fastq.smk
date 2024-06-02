# Rule merges multiple trimmed FASTQ files for each read (R1 and R2) into single files for each sample.
rule merge_fastq:
    input:
        r1=lambda wcs: expand(wcs.output_path + "/" + wcs.sample + "/01_sequence_files/{initial_fq}", \
                       initial_fq=D_sample_details[wcs.sample]["initial_fqs"]["R1"]),
        r2=lambda wcs: expand(wcs.output_path + "/" + wcs.sample + "/01_sequence_files/{initial_fq}", \
                       initial_fq=D_sample_details[wcs.sample]["initial_fqs"]["R2"]),
    output:
        r1="{output_path}/{sample}/02_trim_fastq/{sample}_r1.fq.gz",
        r2="{output_path}/{sample}/02_trim_fastq/{sample}_r2.fq.gz",
    log:
        "{output_path}/{sample}/logs/{sample}_merge_fastq.log",
    threads: 1
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "merge_fastq", threads),
        cpus=1,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        "cat {input.r1} 1> {output.r1} 2>> {log}\n"
        "cat {input.r2} 1> {output.r2} 2>> {log}"

# Rule performs 3' trimming using Trim-galore!, moves Unique Molecular Identifiers (UMIs) to the read name using fastp, and merges the reads.
rule trim_and_merge_fastq:
    input:
        r1="{output_path}/{sample}/02_trim_fastq/{sample}_r1.fq.gz",
        r2="{output_path}/{sample}/02_trim_fastq/{sample}_r2.fq.gz",
    output:
        tr1="{output_path}/{sample}/02_trim_fastq/{sample}_r1_val_1.fq.gz",
        tr2="{output_path}/{sample}/02_trim_fastq/{sample}_r2_val_2.fq.gz",
        combined="{output_path}/{sample}/02_trim_fastq/{sample}_combined.fq.gz",
        uncombined_r1="{output_path}/{sample}/02_trim_fastq/{sample}_r1_uncombined.fq.gz",
        uncombined_r2="{output_path}/{sample}/02_trim_fastq/{sample}_r2_uncombined.fq.gz",
    params:
        umi_info=lambda wcs: get_umi_info(wcs.sample),
        trim_lens=lambda wcs: D_sample_details[wcs.sample]["fq_trim_lengths"],
        parallel=lambda wcs, threads: get_parallel(wcs, 4, threads),
    conda:
        "../envs/fastp.yaml",
    threads: lambda wildcards: get_cpus(1, 32),
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "trim_and_merge", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        'mkdir -p {wildcards.output_path}/{wildcards.sample}/02_trim_fastq\n\n'

        # Step 1: Trim 3' ends using Trim Galore!
        'trim_galore --fastqc --fastqc_args "--noextract" --gzip '
        '--cores {params.parallel} {params.trim_lens} --output_dir '
        '{wildcards.output_path}/{wildcards.sample}/02_trim_fastq --paired {input.r1} {input.r2} '
        '&> {wildcards.output_path}/{wildcards.sample}/logs/{wildcards.sample}_{rule}_trimming.log\n\n'

        # Step 2: Merge trimmed reads using fastp
        'fastp -i {output.tr1} -I {output.tr2} --disable_quality_filtering '
        '--merge --merged_out {output.combined} '
        '--out1 {output.uncombined_r1} --out2 {output.uncombined_r2} '
        '--thread={threads} {params.umi_info} '
        '--json={wildcards.output_path}/{wildcards.sample}/02_trim_fastq/{sample}.json '
        '--html={wildcards.output_path}/{wildcards.sample}/02_trim_fastq/{sample}.html '
        '&> {wildcards.output_path}/{wildcards.sample}/logs/{wildcards.sample}_{rule}_fastp.log\n\n'