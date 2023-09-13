# copy the reference genome to output directory using rclone, allowing slurm jobs to start faster
rule transfer_ref_genome:
    input:
        lambda wcs: D_genomes[wcs.sample],  # Dynamic input based on the sample name
    output:
        # The following output files are produced, including bisulfite-converted genome sequences
        "{output_path}/{sample}/{sample}.fa",
        "{output_path}/{sample}/{sample}.fa.fai",
        "{output_path}/{sample}/{sample}.genome",
        "{output_path}/{sample}/CHH_ROI.bed",
        "{output_path}/{sample}/Bisulfite_Genome/{sample}.CT_conversion.fa",
        "{output_path}/{sample}/Bisulfite_Genome/{sample}.GA_conversion.fa",
        "{output_path}/{sample}/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        "{output_path}/{sample}/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
    log:
        "{output_path}/{sample}/logs/{sample}_transfer_ref_genome.log",  # Log file stored under logs directory for each sample
    conda:
        "../envs/rclone.yaml",  # Specifies the conda environment needed for rclone
    threads: 1,
    resources:
        time_min=30,  # Sets the maximum allowed time for the job in minutes
        partition="--partition=io",  # Specifies a specific partition, in this case 'io'
    shell:
        "rclone sync {input} {wildcards.output_path}/{wildcards.sample}/ --progress "
        "--transfers=8 --log-file={log} --log-level INFO"
