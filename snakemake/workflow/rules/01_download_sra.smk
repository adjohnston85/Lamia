# Download SRA files using aria2c and srapath
rule sra_download:
    output:
        r1="{output_path}/{sample}/01_sequence_files/{r1}/{r1}.sra",
    conda:
        "../envs/sra-tools.yaml"
    threads: 1
    resources:
        # Calculate time limit based on a custom function in master.smk
        time_min=lambda wcs, threads: get_time_min(wcs, wcs, "sra_download", threads),
        mem_mb=512,
        cpus=1,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="--partition=io"
    retries: 5
    shell:
        "mkdir -p {wildcards.output_path}/{wildcards.sample}/01_sequence_files \n\n"
        "aria2c -c --auto-file-renaming=false -x 16 --dir=/ -o {output.r1} $(srapath "
        "{wildcards.r1}) &>> {wildcards.output_path}/{wildcards.sample}/logs/{wildcards.sample}_{rule}.log"

# Convert SRA files to FASTQ format using parallel-fastq-dump
rule sra_to_fastq:
    input:
        r1="{output_path}/{sample}/01_sequence_files/{r1}/{r1}.sra",
    output:
        r1="{output_path}/{sample}/01_sequence_files/{r1}/{r1}_1.fastq.gz",
        r2="{output_path}/{sample}/01_sequence_files/{r1}/{r1}_2.fastq.gz",
    conda:
        "../envs/sra-tools.yaml"
    # Dynamic assignment of threads based on a custom function in master.smk
    threads: lambda wcs: get_cpus(1,16)
    resources:
        # Calculate time limit based on a custom function in master.smk
        time_min=lambda wcs, input, threads: get_time_min(wcs, wcs, "sra_to_fastq", threads),
        # Fetch memory requirements from a custom function in master.smk
        mem_mb=lambda wcs, threads:get_mem_mb(wcs, threads, 2048),
        cpus=lambda wcs, threads: threads,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="--partition=io"
    shell:
        "parallel-fastq-dump --split-files --threads {resources.cpus} --gzip "
        "--outdir {wildcards.output_path}/{wildcards.sample}/01_sequence_files/{wildcards.r1} --sra-id {input} "
        "--tmpdir {wildcards.output_path}/{wildcards.sample}/01_sequence_files/{wildcards.r1} "
        "&>> {wildcards.output_path}/{wildcards.sample}/logs/{wildcards.sample}_{rule}.log \n\n"
