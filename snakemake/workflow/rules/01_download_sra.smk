rule sra_download:
    output:
        r1="{output_path}/{sample}/01_sequence_files/{r1}/{r1}.sra",
    conda:
        "../envs/sra-tools.yaml"
    threads: 1
    resources:
        time_min=lambda wcs, threads: get_time_min(wcs, wcs, "sra_download", threads),
        mem_mb=512,
        cpus=1,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="--partition=io"
    retries: 5
    shell:
        "aria2c -c --auto-file-renaming=false -x 16 --dir=/ -o {output.r1} $(srapath "
        "{wildcards.r1}) &>> {wildcards.output_path}/{wildcards.sample}/logs/{wildcards.sample}_{rule}.log"


rule sra_to_fastq:
    input:
        "{output_path}/{sample}/01_sequence_files/{run_accession}/{run_accession}.sra",
    output:
        "{output_path}/{sample}/01_sequence_files/{run_accession}_1.fastq.gz",
        "{output_path}/{sample}/01_sequence_files/{run_accession}_2.fastq.gz",
    conda:
        "../envs/sra-tools.yaml"
    threads: lambda wcs: get_cpus(1,16)
    resources:
        time_min=lambda wcs, threads: get_time_min(wcs, wcs, "sra_to_fastq", threads),
        mem_mb=get_mem_mb,
        cpus=lambda wcs, threads: threads,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="--partition=io"
    shell:
        "parallel-fastq-dump --split-files --threads {resources.cpus} "
        "--gzip --outdir {wildcards.output_path}/{wildcards.sample}/01_sequence_files/ --sra-id {input} "
        "--tmpdir {wildcards.output_path}/{wildcards.sample}/01_sequence_files/ "
        "&>> {wildcards.output_path}/{wildcards.sample}/logs/{wildcards.sample}_{rule}.log"
