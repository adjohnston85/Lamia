rule softlink_fastq:
    input:
        lambda wcs: D_sample_details[wcs.sample]["run_files"][wcs.fastq_name],
    output:
        "{output_path}/{sample}/01_sequence_files/{fastq_name}",
    threads: 1
    resources:
        time_min=5,
        mem_mb=512,
        cpus=1,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="--partition=io",
    shell:
        "ln -sf {input} {output} &>> {wildcards.output_path}/{wildcards.sample}/logs/{wildcards.sample}_{rule}.log"
