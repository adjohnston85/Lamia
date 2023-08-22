rule bismark_align:
    input:
        r1="{output_path}/{sample}/02_trim_fastq/{sample}_r1.fq.gz",
        r2="{output_path}/{sample}/02_trim_fastq/{sample}_r2.fq.gz",
    output:
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.bam",
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_PE_report.txt",
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.nucleotide_stats.txt",
    params:
        genome_path = lambda wcs: D_sample_details[wcs.sample]["genome_path"],
        parallel = lambda wcs, threads: get_parallel(wcs, 5, threads),
        non_directional = lambda wcs: "--non_directional " if D_sample_details[wcs.sample]["non_directional"] else ""
    log:
        "{output_path}/{sample}/logs/{sample}_bismark_align.log",
    conda:
        "../envs/bismark.yaml"
    threads: lambda wcs: get_cpus(16,64)
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "bismark_align", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        "bismark --genome_folder {params.genome_path} --multicore {params.parallel} --bowtie2 {params.non_directional}"
        "--temp_dir {wildcards.output_path}/{wildcards.sample}/03_align_fastq -1 {input.r1} -2 {input.r2} "
        "--nucleotide_coverage --output_dir {wildcards.output_path}/{wildcards.sample}/03_align_fastq 1>> /dev/null 2>> {log}"


rule sort_bam:
    input:
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.bam",
    output:
        bam="{output_path}/{sample}/03_align_fastq/{sample}_s.bam",
        bai="{output_path}/{sample}/03_align_fastq/{sample}_s.bam.bai",
        flag="{output_path}/{sample}/03_align_fastq/{sample}_s_flagstat.txt",
    log:
        "{output_path}/{sample}/logs/{sample}_sort_bam.log",
    conda:
        "../envs/samtools.yaml"
    threads: lambda wcs: get_cpus(1,16),
    resources:
        time_min=lambda wcs, input, threads : get_time_min(wcs, input, "sort_bam", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        "samtools sort -@ {threads} -O bam -o {output.bam} {input} &>> {log} \n\n"

        "samtools index -@ {threads} {output.bam} &>> {log} \n\n"

        "samtools flagstat -@ {threads} {output.bam} 1> {output.flag} 2>> {log}"
