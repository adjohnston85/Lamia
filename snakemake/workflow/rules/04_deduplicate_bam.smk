# Rule to deduplicate BAM files and generate genome-specific BAM and deduplicated BAM files.
rule deduplicate_bam:
    input:
        "{output_path}/{sample}/03_align_fastq/{sample}_s.bam",  # Source BAM file
    output:
        mG_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_{genome}_genome.bam",  # Genome-specific BAM
        mG_dedup_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_{genome}_genome_sd.bam",  # Deduplicated BAM
    params:
        # Conditionally add BED file for region of interest if available
        b=lambda wcs: '-b ' + D_sample_details[wcs.sample]["roi_bed"] + ' '
          if "roi_bed" in D_sample_details[wcs.sample] else "",
        # Path to reference genome's conversion FASTA
        r=lambda wcs: D_sample_details[wcs.sample]["genome_path"] + \
          "/Bisulfite_Genome/" + D_sample_details[wcs.sample]["genome"] + \
          "." + wcs.genome + "_conversion.fa",
        # UMI prefix dynamically obtained based on the sample
        umi_prefix=lambda wcs: get_umi_prefix(wcs.sample),
    log:
        "{output_path}/{sample}/logs/{sample}_deduplicate_{genome}_bam.log",
    conda:
        "../envs/gencore.yaml",
    threads:
        lambda wcs: get_cpus(1,64),  # Dynamic CPU allocation. Max limit is 64.
    resources:
        # Dynamic resource allocation and user/account specific parameters
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "deduplicate_bam", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        # Extract only the records for the CT or GA genome and create a new BAM file
        'samtools view -h {input} | grep -e "XG:Z:{wildcards.genome}" -e "@" '
        '| samtools view -bS - 1> {output.mG_bam} '
        '2>> {log} \n\n'

        # Run gencore for deduplication
        'gencore {params.umi_prefix}-s 1 -d 2 '
        '-i {output.mG_bam} -o {output.mG_dedup_bam} -r {params.r} {params.b}'
        '-j {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam/gencore_{wildcards.genome}.json '
        '-h {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam/gencore_{wildcards.genome}.html '
        '&>> {log}'


# Rule to merge CT and GA deduplicated BAM files, sort them, create index, and generate flagstat report.
rule merge_deduplicate_bams:
    input:
        CT="{output_path}/{sample}/04_deduplicate_bam/{sample}_CT_genome_sd.bam",  # CT genome-specific deduplicated BAM
        GA="{output_path}/{sample}/04_deduplicate_bam/{sample}_GA_genome_sd.bam",  # GA genome-specific deduplicated BAM
    output:
        bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_sd.bam",  # Merged and sorted BAM
        bai="{output_path}/{sample}/04_deduplicate_bam/{sample}_sd.bam.bai",  # BAM index
        flag="{output_path}/{sample}/04_deduplicate_bam/{sample}_sd_flagstat.txt",  # Flagstat report
    log:
        "{output_path}/{sample}/logs/{sample}_merge_deduplicate_bams.log",
    conda:
        "../envs/samtools.yaml",
    threads:
        lambda wcs: get_cpus(1,64),  # Dynamic CPU allocation with max limit of 64
    resources:
        # Dynamic resource allocation and user/account specific parameters
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "merge_deduplicate_bam", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        # Merge CT and GA BAMs, sort the merged BAM and save it.
        "samtools merge -@ {threads} - {input.CT} {input.GA} | "
        "samtools sort -m 2G -O bam -@ 12 -o {output.bam} - &>> {log} \n\n"

        # Generate BAM index file
        "samtools index -@ {threads} {output.bam} &>> {log} \n\n"

        # Generate Flagstat report for the final BAM
        "samtools flagstat {output.bam} 1> {output.flag} 2>> {log}"
