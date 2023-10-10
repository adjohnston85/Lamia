# Rule to deduplicate BAM files and generate genome-specific BAM and deduplicated BAM files.
rule deduplicate_bam:
    input:
        "{output_path}/{sample}/03_align_fastq/{sample}_s.bam",  # Source BAM file
    output:
        bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_sd.bam",  # Merged and sorted BAM
        bai="{output_path}/{sample}/04_deduplicate_bam/{sample}_sd.bam.bai",  # BAM index
        flag="{output_path}/{sample}/04_deduplicate_bam/{sample}_sd_flagstat.txt",  # Flagstat report
        json="{output_path}/{sample}/04_deduplicate_bam/gencore.json",
        html="{output_path}/{sample}/04_deduplicate_bam/gencore.html",
    params:
        # Conditionally add BED file for region of interest if available
        b=lambda wcs: '-b ' + D_sample_details[wcs.sample]["roi_bed"] + ' '
            if "roi_bed" in D_sample_details[wcs.sample] else "",
        # UMI prefix dynamically obtained based on the sample
        umi_prefix=lambda wcs: get_umi_prefix(wcs.sample),
        r=lambda wcs: D_sample_details[wcs.sample]["genome_path"] + \
            "/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
    log:
        "{output_path}/{sample}/logs/{sample}_deduplicate_bam.log",
    conda:
        "../envs/gencore.yaml",
    threads:
        lambda wcs: get_cpus(1,64),  # Dynamic CPU allocation. Max limit is 64.
    resources:
        # Dynamic resource allocation and user/account specific parameters
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "deduplicate_bam", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        "mkdir -p {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam \n\n"

        # Run gencore for deduplication
        'gencore {params.umi_prefix}--supporting_reads=1 '
        '--umi_diff_threshold=4 -i {input} -r {params.r} {params.b}'
        '-j {output.json} '
        '-h {output.html} '
        '2>> {log} | samtools sort - -o {output.bam} \n\n'

        # Generate BAM index file
        "samtools index -@ {threads} {output.bam} &>> {log} \n\n"

        # Generate Flagstat report for the final BAM
        "samtools flagstat {output.bam} 1> {output.flag} 2>> {log}"
