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
    threads: lambda wcs: get_cpus(16,64)
    resources:
        # Resource constraints, dynamically calculated
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "bismark_align", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        # Bismark alignment command
        "bismark --genome_folder {input.genome_path} --multicore {params.parallel} --bowtie2 {params.non_directional}"
        "--temp_dir {wildcards.output_path}/{wildcards.sample}/03_align_fastq -1 {input.r1} -2 {input.r2} "
        "--maxins {params.maxins} "
        "--nucleotide_coverage --output_dir {wildcards.output_path}/{wildcards.sample}/03_align_fastq 1>> /dev/null 2>> {log}"


# Rule for sorting BAM files and generating index and flagstat
rule sort_bam:
    input:
        # Input BAM file from bismark_align rule
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.bam",
    output:
        # Sorted BAM file and its index
        bam="{output_path}/{sample}/03_align_fastq/{sample}_s.bam",
        bai="{output_path}/{sample}/03_align_fastq/{sample}_s.bam.bai",
        # SAMtools flagstat report
        flag="{output_path}/{sample}/03_align_fastq/{sample}_s_flagstat.txt",
    log:
        # Log file to capture the stdout and stderr
        "{output_path}/{sample}/logs/{sample}_sort_bam.log",
    threads: lambda wcs: get_cpus(1,16)
    resources:
        # Resource constraints, dynamically calculated
        time_min=lambda wcs, input, threads : get_time_min(wcs, input, "sort_bam", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        # SAMtools commands for sorting, indexing and generating flagstats
        "samtools sort -@ {threads} -O bam -o {output.bam} {input} &>> {log} \n\n"
        "samtools index -@ {threads} {output.bam} &>> {log} \n\n"
        "samtools flagstat -@ {threads} {output.bam} 1> {output.flag} 2>> {log}"


# Rule to collect Picard metrics such as insert size and hybrid selection metrics.
rule picard_metrics:
    input:
        # Input BAM file from the bismark_align rule
        bam="{output_path}/{sample}/03_align_fastq/{sample}_r1_bismark_bt2_pe.bam",
    output:
        # Output histogram PDF for insert size
        histogram="{output_path}/{sample}/03_align_fastq/{sample}_insert_size_histogram.pdf",
        # Output text file containing insert size metrics
        insert_size="{output_path}/{sample}/03_align_fastq/{sample}_insert_size_metrics.txt",
        # Output text file for hybrid selection metrics
        hybrid_selection="{output_path}/{sample}/03_align_fastq/{sample}_hs_metrics.txt",
    params:
        # Path to the regions of interest (ROI) BED file, defaults to an empty string
        roi_bed=lambda wcs: D_sample_details[wcs.sample].get("roi_bed",""),
        # Path to the reference genome
        genome_path = lambda wcs: D_sample_details[wcs.sample]["genome_path"],
        # Reference genome name
        genome = lambda wcs: D_sample_details[wcs.sample]["genome"],
        # Interval list for ROI, constructed dynamically from roid_bed based on output path and sample name
        interval_list = lambda wcs: "{wcs.output_path}/{wcs.sample}/03_align_fastq/ROI.interval_list",
    log:
        # Log file to capture the stdout and stderr of the Picard tool
        "{output_path}/{sample}/logs/{sample}_picard_metrics.log",
    threads: 1  # Single-threaded as Picard tools often are
    resources:
        # Resource constraints dynamically calculated
        time_min=lambda wcs, input, threads : get_time_min(wcs, input, "picard_metrics", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        """
        # Collect Insert Size Metrics
        picard CollectInsertSizeMetrics I={input.bam} O={output.insert_size} H={output.histogram}
        
        # If a regions of interest (ROI) BED file is provided
        if [ -n {params.roi_bed} ]; then
            # Convert BED to Interval List
            picard BedToIntervalList I={params.roi_bed} O={params.interval_list} SD={params.genome_path}/{params.genome}.dict
            
            # Collect Hybrid Selection Metrics
            picard CollectHsMetrics I={input.bam} O={output.hybrid_selection} \
            BAIT_INTERVALS={params.interval_list} TARGET_INTERVALS={params.interval_list}
        fi
        """
