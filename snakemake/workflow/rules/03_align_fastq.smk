# Rule for bismark alignment of bisulfite or enzymatically converted DNA reads
rule bismark_align:
    input:
        # Input FASTQ files for forward and reverse reads
        r1="{output_path}/{sample}/02_trim_fastq/{sample}_r1.fq.gz",
        r2="{output_path}/{sample}/02_trim_fastq/{sample}_r2.fq.gz",
        # The path to the reference genome
        genome_path = lambda wcs: D_sample_details[wcs.sample]["genome_path"],
    output:
        # Aligned BAM and Bismark report files
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.bam",
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_PE_report.txt",
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.nucleotide_stats.txt",
    params:
        # Dynamic core allocation using threads divided by five due to how Bismark parallel works
        parallel = lambda wcs, threads: get_parallel(wcs, 5, threads),
        # Flag for non-directional libraries
        non_directional = lambda wcs: "--non_directional " if D_sample_details[wcs.sample]["non_directional"] else "",
        # Max insert size - Bismark default is 500, Majel automatically sets to 1000 for em-seq library_tpe
        maxins = lambda wcs: D_sample_details[wcs.sample]["maxins"],
    log:
        # Log file to capture the stdout and stderr
        "{output_path}/{sample}/logs/{sample}_bismark_align.log",
    conda:
        "../envs/bismark.yaml",
    threads: lambda wcs: get_cpus(16,64)
    resources:
        # Resource constraints, dynamically calculated
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "bismark_align", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 4096),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        "mkdir -p {wildcards.output_path}/{wildcards.sample}/03_align_fastq \n\n"
        "bismark --genome_folder {input.genome_path} --multicore {params.parallel} --bowtie2 {params.non_directional}"
        "--temp_dir {wildcards.output_path}/{wildcards.sample}/03_align_fastq -1 {input.r1} -2 {input.r2} "
        "--maxins {params.maxins} "
        "--nucleotide_coverage --output_dir {wildcards.output_path}/{wildcards.sample}/03_align_fastq 1>> /dev/null 2>> {log}"


# Rule for sorting BAM files and generating index and flagstat
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
        "../envs/samtools.yaml",
    threads: lambda wcs: get_cpus(1,16)
    resources:
        time_min=lambda wcs, input, threads : get_time_min(wcs, input, "sort_bam", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        "samtools sort -@ {threads} -O bam -o {output.bam} {input} &>> {log} \n\n"
        "samtools index -@ {threads} {output.bam} &>> {log} \n\n"
        "samtools flagstat -@ {threads} {output.bam} 1> {output.flag} 2>> {log}"


# Rule to collect Picard metrics such as insert size and hybrid selection metrics.
rule picard_metrics:
    input:
        bam="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.bam",
    output:
        histogram="{output_path}/{sample}/03_align_fastq/{sample}_insert_size_histogram.pdf",
        insert_size="{output_path}/{sample}/03_align_fastq/{sample}_insert_size_metrics.txt",
        hybrid_selection="{output_path}/{sample}/03_align_fastq/{sample}_hs_metrics.txt",
    params:
        roi_bed=lambda wcs: D_sample_details[wcs.sample].get("roi_bed",""),
        genome_path = lambda wcs: D_sample_details[wcs.sample]["genome_path"],
        genome = lambda wcs: D_sample_details[wcs.sample]["genome"],
        interval_list = lambda wcs: "{}/{}/03_align_fastq/ROI.interval_list".format(wcs.output_path, wcs.sample),
    log:
        "{output_path}/{sample}/logs/{sample}_picard_metrics.log",
    conda:
        "../envs/picard.yaml",
    threads: 1
    resources:
        time_min=lambda wcs, input, threads : get_time_min(wcs, input, "picard_metrics", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        """
        picard CollectInsertSizeMetrics -I {input.bam} -O {output.insert_size} -H {output.histogram} &>> {log}
        
        if [ -n '{params.roi_bed}' ]; then
            picard BedToIntervalList -I {params.roi_bed} -O {params.interval_list} \
            -SD {params.genome_path}/{params.genome}.dict &>> {log}
            
            picard CollectHsMetrics -I {input.bam} -O {output.hybrid_selection} \
            -BI {params.interval_list} -TI {params.interval_list} &>> {log}
        else
            echo 'no regions of interest specified' > {output.hybrid_selection}
        fi
        """
