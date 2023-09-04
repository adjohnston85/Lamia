# Rule moves Unique Molecular Identifiers (UMIs) to the read name using the fastp tool.
rule move_umis:
    input:
        # Input files come from a previous step located in the '01_sequence_files' directory.
        # Uses a variety of wildcard constraints for flexibility.
        r1="{output_path}/{sample}/01_sequence_files/{prefix}_{read}1{suffix}.fastq.gz",
        r2="{output_path}/{sample}/01_sequence_files/{prefix}_{read}2{suffix}.fastq.gz",

    output:
        # Output files are stored in '02_trim_fastq' with a '_fastp_1' and '_fastp_2' tag to denote processed files.
        r1="{output_path}/{sample}/02_trim_fastq/{prefix}_{read}1{suffix}_fastp_1.fq.gz",
        r2="{output_path}/{sample}/02_trim_fastq/{prefix}_{read}2{suffix}_fastp_2.fq.gz",

    # Using constrained wildcards for the 'read' and 'suffix' to avoid ambiguity.
    wildcard_constraints:
        read = '[rR]?',
        suffix = '_?.*',

    # UMI information is dynamically obtained via a custom function called 'get_umi_info'.
    params:
        umi_info=lambda wcs: get_umi_info(wcs.sample),

    conda:
        "../envs/fastp.yaml",

    # Dynamic CPU allocation. Max limit is 16.
    threads: lambda wildcards: get_cpus(1,16),
    
    # Resource constraints and user/account specific parameters.
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "move_umis", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    
    # Shell command to run fastp with various flags and arguments, including capturing JSON and HTML reports.
    shell:
        'fastp -i {input.r1} -I {input.r2} '
        '-o {output.r1} -O {output.r2} --thread={threads} -Q {params.umi_info} '
        '--json={wildcards.output_path}/{wildcards.sample}/02_trim_fastq/{wildcards.prefix}{wildcards.suffix}.json '
        '--html={wildcards.output_path}/{wildcards.sample}/02_trim_fastq/{wildcards.prefix}{wildcards.suffix}.html '
        '&> {wildcards.output_path}/{wildcards.sample}/logs/{wildcards.sample}_{rule}.log'


# Rule performs quality and adapter trimming on paired-end FASTQ files using trim_galore.
rule trim_fastq:
    input:
        r1="{output_path}/{sample}/02_trim_fastq/{prefix}_{read}1{suffix}_fastp_1.fq.gz",
        r2="{output_path}/{sample}/02_trim_fastq/{prefix}_{read}2{suffix}_fastp_2.fq.gz",
    output:
        r1="{output_path}/{sample}/02_trim_fastq/{prefix}_{read}1{suffix}_fastp_1_val_1.fq.gz",
        r2="{output_path}/{sample}/02_trim_fastq/{prefix}_{read}2{suffix}_fastp_2_val_2.fq.gz",
    wildcard_constraints:
        read = '[rR]?',
        suffix = '_?.*',
    params:
        trim_lens=lambda wcs: D_sample_details[wcs.sample]["trim_lengths"],
        parallel=lambda wcs, threads: get_parallel(wcs, 4, threads),
    conda:
        "../envs/trim_galore.yaml"
    threads: lambda wildcards: get_cpus(1,32) # Dynamic CPU assignment. Max limit is 32.
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "trim_fastq", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        'trim_galore --fastqc --fastqc_args "--noextract" --gzip '
        '--cores {params.parallel} {params.trim_lens}--output_dir '
        '{wildcards.output_path}/{wildcards.sample}/02_trim_fastq --paired {input.r1} {input.r2} '
        '&> {wildcards.output_path}/{wildcards.sample}/logs/{wildcards.sample}_{rule}.log'


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
    threads: 1 # Single-threaded operation as it's just a file merge.
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "merge_fastq", threads),
        cpus=1,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        "cat {input.r1} 1> {output.r1} 2>> {log} \n"
        "cat {input.r2} 1> {output.r2} 2>> {log}"

