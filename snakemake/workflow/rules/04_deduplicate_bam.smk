rule deduplicate_bam:
    input:
        "{output_path}/{sample}/03_align_fastq/{sample}_s.bam",
    output:
        mG_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_{genome}_genome.bam",
        mG_dedup_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_{genome}_genome_sd.bam",
    params:
        b=lambda wcs: '-b ' + D_sample_details[wcs.sample]["roi_bed"] + ' '
          if "roi_bed" in D_sample_details[wcs.sample] else "",
        r=lambda wcs: D_sample_details[wcs.sample]["genome_path"] + \
          "/Bisulfite_Genome/" + D_sample_details[wcs.sample]["genome"] + \
          "." + wcs.genome + "_conversion.fa",
        umi_info= lambda wcs: get_umi_info(wcs.sample),
    log:
        "{output_path}/{sample}/logs/{sample}_deduplicate_{genome}_bam.log",
    conda:
        "../envs/gencore.yaml"
    threads: lambda wcs: get_cpus(1,64)
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "deduplicate_bam", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        'samtools view -h {input} | grep -e "XG:Z:{wildcards.genome}" -e "@" '
        '| samtools view -bS - 1> {output.mG_bam} '
        '2>> {log} \n\n'

        'gencore {params.umi_info}-s 1 -d 2 '
        '-i {output.mG_bam} -o {output.mG_dedup_bam} -r {params.r} {params.b}'
        '-j {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam/gencore_{wildcards.genome}.json '
        '-h {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam/gencore_{wildcards.genome}.html '
        '&>> {log}'


rule merge_deduplicate_bams:
    input:
        CT="{output_path}/{sample}/04_deduplicate_bam/{sample}_CT_genome_sd.bam",
        GA="{output_path}/{sample}/04_deduplicate_bam/{sample}_GA_genome_sd.bam",
    output:
        bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_sd.bam",
        bai="{output_path}/{sample}/04_deduplicate_bam/{sample}_sd.bam.bai",
        flag="{output_path}/{sample}/04_deduplicate_bam/{sample}_sd_flagstat.txt",
    log:
        "{output_path}/{sample}/logs/{sample}_merge_deduplicate_bams.log",
    conda:
        "../envs/samtools.yaml"
    threads: lambda wcs: get_cpus(1,64)
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "merge_deduplicate_bam", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        "samtools merge -@ {threads} - {input.CT} {input.GA} | "
        "samtools sort -m 2G -O bam -@ 12 -o {output.bam} - &>> {log} \n\n"

        "samtools index -@ {threads} {output.bam} &>> {log} \n\n"

        "samtools flagstat {output.bam} 1> {output.flag} 2>> {log}"
