rule transfer_ref_genome:
    input:
        lambda wcs: D_genomes[wcs.sample],
    output:
        "{output_path}/{sample}/{sample}.fa",
        "{output_path}/{sample}/{sample}.fa.fai",
        "{output_path}/{sample}/{sample}.genome",
        "{output_path}/{sample}/CHH_ROI.bed",
        "{output_path}/{sample}/Bisulfite_Genome/{sample}.CT_conversion.fa",
        "{output_path}/{sample}/Bisulfite_Genome/{sample}.GA_conversion.fa",
        "{output_path}/{sample}/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        "{output_path}/{sample}/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
    log:
        "{output_path}/{sample}/logs/{sample}_transfer_ref_genome.log",
    conda:
        "../envs/rclone.yaml"
    threads: 8
    resources:
        time_min=30,
        mem_mb=2048,
        cpus=8,
        email=config["email"],
        account=config["account"],
        partition="--partition=io"
    shell:
        "rclone sync {input} {wildcards.output_path}/{wildcards.sample}/ --progress "
        "--transfers={threads} --log-file={log} --log-level INFO"
