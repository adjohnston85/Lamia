# Rule moves Unique Molecular Identifiers (UMIs) to the read name using the fastp tool and trims using Trim-galore!.
rule trim_fastq:
    input:
        r1="{output_path}/{sample}/01_sequence_files/{prefix}_{read}1{suffix}.fastq.gz",
        r2="{output_path}/{sample}/01_sequence_files/{prefix}_{read}2{suffix}.fastq.gz",
    output:
        tr1="{output_path}/{sample}/02_trim_fastq/{prefix}_{read}1{suffix}_fastp_1_val_1.fq.gz",
        tr2="{output_path}/{sample}/02_trim_fastq/{prefix}_{read}2{suffix}_fastp_2_val_2.fq.gz",
    wildcard_constraints:
        read = '[rR]?',
        suffix = '_?.*',
    params:
        fr1="{output_path}/{sample}/02_trim_fastq/{prefix}_{read}1{suffix}_fastp_1.fq.gz",
        fr2="{output_path}/{sample}/02_trim_fastq/{prefix}_{read}2{suffix}_fastp_2.fq.gz",
        umi_info=lambda wcs: get_umi_info(wcs.sample),
        trim_lens=lambda wcs: D_sample_details[wcs.sample]["trim_lengths"],
        parallel=lambda wcs, threads: get_parallel(wcs, 4, threads),
    conda:
        "../envs/fastp.yaml",
    threads: lambda wildcards: get_cpus(1,32),
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "move_umis", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        'mkdir -p {wildcards.output_path}/{wildcards.sample}/02_trim_fastq\n\n'

        'fastp -i {input.r1} -I {input.r2} --disable_quality_filtering '
        '-o {params.fr1} -O {params.fr2} --thread={threads} {params.umi_info}'
        '--json={wildcards.output_path}/{wildcards.sample}/02_trim_fastq/{wildcards.prefix}{wildcards.suffix}.json '
        '--html={wildcards.output_path}/{wildcards.sample}/02_trim_fastq/{wildcards.prefix}{wildcards.suffix}.html '
        '&> {wildcards.output_path}/{wildcards.sample}/logs/{wildcards.sample}_{rule}.log\n\n'

        'trim_galore --fastqc --fastqc_args "--noextract" --gzip '
        '--cores {params.parallel} {params.trim_lens}--output_dir '
        '{wildcards.output_path}/{wildcards.sample}/02_trim_fastq --paired {params.fr1} {params.fr2} '
        '&> {wildcards.output_path}/{wildcards.sample}/logs/{wildcards.sample}_{rule}.log\n\n'
        
        'rm {params.fr1} {params.fr2}'

# Rule merges multiple trimmed FASTQ files for each read (R1 and R2) into single files for each sample.
rule merge_fastq:
    input:
        r1=lambda wcs: expand(wcs.output_path + "/" + wcs.sample + "/02_trim_fastq/{trimmed_fq}", \
                       trimmed_fq=D_sample_details[wcs.sample]["trimmed_fqs"]["_fastp_1_val_1"]),
        r2=lambda wcs: expand(wcs.output_path + "/" + wcs.sample + "/02_trim_fastq/{trimmed_fq}", \
                       trimmed_fq=D_sample_details[wcs.sample]["trimmed_fqs"]["_fastp_2_val_2"]),
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
