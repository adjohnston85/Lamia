rule bismark_align:
    input:
        # Input FASTQ files for forward and reverse reads
        uncombined_r1="{output_path}/{sample}/02_trim_fastq/{sample}_r1_uncombined.fq.gz",
        uncombined_r2="{output_path}/{sample}/02_trim_fastq/{sample}_r2_uncombined.fq.gz",
        # Input combined FASTQ file
        combined="{output_path}/{sample}/02_trim_fastq/{sample}_combined.fq.gz",
        # The path to the reference genome
        genome_path=lambda wcs: D_sample_details[wcs.sample]["genome_path"],
    output:
        # Aligned BAM and Bismark report files for paired-end reads
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_uncombined_bismark_bt2_PE_report.txt",
        "{output_path}/{sample}/03_align_fastq/{sample}_r1_uncombined_bismark_bt2_pe.nucleotide_stats.txt",
        # Aligned BAM and Bismark report files for single-end reads
        "{output_path}/{sample}/03_align_fastq/{sample}_combined_bismark_bt2_SE_report.txt",
        "{output_path}/{sample}/03_align_fastq/{sample}_combined_bismark_bt2.nucleotide_stats.txt",
        bam_se="{output_path}/{sample}/03_align_fastq/{sample}_combined_bismark_bt2.bam",
        bam_pe="{output_path}/{sample}/03_align_fastq/{sample}_r1_uncombined_bismark_bt2_pe.bam",
    params:
        # Dynamic core allocation using threads divided by five due to how Bismark parallel works
        parallel=lambda wcs, threads: get_parallel(wcs, 5, threads),
        non_directional=lambda wcs: "--non_directional " if D_sample_details[wcs.sample]["non_directional"] else "",
        # Max insert size - Bismark default is 500, Majel automatically sets to 1000 for em-seq library_type
        maxins=lambda wcs: D_sample_details[wcs.sample]["maxins"],
    log:
        # Log file to capture the stdout and stderr
        "{output_path}/{sample}/logs/{sample}_bismark_align.log",
    conda:
        "../envs/bismark.yaml",
    threads: lambda wcs: get_cpus(16, 64),
    resources:
        # Resource constraints, dynamically calculated
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "bismark_align", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs,threads: get_mem_mb(wcs, threads, 4096),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        'mkdir -p {wildcards.output_path}/{wildcards.sample}/03_align_fastq \n\n'

        'bismark '
            '--genome_folder {input.genome_path} '
            '--multicore {params.parallel} '
            '--bowtie2 {params.non_directional} '
            '--temp_dir {wildcards.output_path}/{wildcards.sample}/03_align_fastq '
            '-1 {input.uncombined_r1} '
            '-2 {input.uncombined_r2} '
            '--maxins {params.maxins} '
            '--nucleotide_coverage '
            '--output_dir {wildcards.output_path}/{wildcards.sample}/03_align_fastq '
            '1>> /dev/null 2>> {log}\n\n'

        'bismark '
            '--genome_folder {input.genome_path} '
            '--multicore {params.parallel} '
            '--bowtie2 {params.non_directional} '
            '--temp_dir {wildcards.output_path}/{wildcards.sample}/03_align_fastq '
            '{input.combined} '
            '--nucleotide_coverage '
            '--output_dir {wildcards.output_path}/{wildcards.sample}/03_align_fastq '
            '1>> /dev/null 2>> {log}\n\n'

# Rule to fix mate information in BAM files using Picard
rule fix_mate_information:
    input:
        bam_pe="{output_path}/{sample}/03_align_fastq/{sample}_r1_uncombined_bismark_bt2_pe.bam",
    output:
        fixed_bam="{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info.bam",
    log:
        "{output_path}/{sample}/logs/{sample}_fix_mate_information.log",
    conda:
        "../envs/picard.yaml",
    threads: 1,
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 10, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        'picard FixMateInformation '
            '-INPUT {input.bam_pe} '
            '-OUTPUT {output.fixed_bam} '
            '-TMP_DIR {wildcards.output_path}/{wildcards.sample}/tmp '
            '-VALIDATION_STRINGENCY SILENT '
            '1> {log} 2>&1'

rule umi_extraction_to_rx:
    input:
        bam_se="{output_path}/{sample}/03_align_fastq/{sample}_combined_bismark_bt2.bam",
        bam_pe="{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info.bam"
    output:
        bam_with_umi_se="{output_path}/{sample}/03_align_fastq/{sample}_combined_umi.bam",
        bam_with_umi_pe="{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info_umi.bam"
    log:
        "{output_path}/{sample}/logs/{sample}_umi_extraction_to_rx.log",
    conda:
        "../envs/samtools.yaml",
    threads: 1,
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 40, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        """
        samtools view -h {input.bam_se} | \
        awk 'BEGIN {{FS=OFS="\t"}} {{
            if ($1 ~ /^@/) {{print; next}}  # Pass headers directly to output
            if (match($1, /UMI_([A-Z]+)_([A-Z]+)_/, arr)) {{
                umi = arr[1] "-" arr[2];
                sub(/:UMI_[A-Z]+_[A-Z]+_/, ":", $1);
                print $0, "RX:Z:" umi;  # Append the new RX tag
            }} else {{
                print;  # Output lines without UMI in the header unchanged
            }}
        }}' | samtools view -b -o {output.bam_with_umi_se}.temp

        samtools view -h {input.bam_pe} | \
        awk 'BEGIN {{FS=OFS="\t"}} {{
            if ($1 ~ /^@/) {{print; next}}  # Pass headers directly to output
            if (match($1, /UMI_([A-Z]+)_([A-Z]+)_/, arr)) {{
                umi = arr[1] "-" arr[2];
                sub(/:UMI_[A-Z]+_[A-Z]+_/, ":", $1);
                print $0, "RX:Z:" umi;  # Append the new RX tag
            }} else {{
                print;  # Output lines without UMI in the header unchanged
            }}
        }}' | samtools view -b -o {output.bam_with_umi_pe}.temp

        mv {output.bam_with_umi_se}.temp {output.bam_with_umi_se}
        mv {output.bam_with_umi_pe}.temp {output.bam_with_umi_pe}
        """

rule filter_non_conversion:
    input:
        se_bam=lambda wildcards: (
            "{output_path}/{sample}/03_align_fastq/{sample}_combined_umi.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
            if "umi_len" in D_sample_details[wildcards.sample]
            else "{output_path}/{sample}/03_align_fastq/{sample}_combined_bismark_bt2.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
        ),
        pe_bam=lambda wildcards: (
            "{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info_umi.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
            if "umi_len" in D_sample_details[wildcards.sample]
            else "{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
        )
    output:
        nonCG_filtered_se="{output_path}/{sample}/03_align_fastq/{sample}_combined_nonCG_filtered.bam",
        nonCG_filtered_pe="{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info_nonCG_filtered.bam"
    log:
        "{output_path}/{sample}/logs/{sample}_filter_non_conversion.log",
    conda:
        "../envs/bismark.yaml",
    threads: 1,
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 20, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        """
        filter_non_conversion --percentage_cutoff 10 {input.se_bam} &>> {log}

        if [ "{input.se_bam}" = "{wildcards.output_path}/{wildcards.sample}/03_align_fastq/{wildcards.sample}_combined_umi.bam" ]; then
            mv {wildcards.output_path}/{wildcards.sample}/03_align_fastq/{wildcards.sample}_combined_umi.nonCG_filtered.bam {output.nonCG_filtered_se}
        else
            mv {wildcards.output_path}/{wildcards.sample}/03_align_fastq/{wildcards.sample}_combined_bismark_bt2.nonCG_filtered.bam {output.nonCG_filtered_se}
        fi

        filter_non_conversion -p --percentage_cutoff 10 {input.pe_bam} &>> {log}

        if [ "{input.pe_bam}" = "{wildcards.output_path}/{wildcards.sample}/03_align_fastq/{wildcards.sample}_fixed_mate_info_umi.bam" ]; then
            mv {wildcards.output_path}/{wildcards.sample}/03_align_fastq/{wildcards.sample}_fixed_mate_info_umi.nonCG_filtered.bam {output.nonCG_filtered_pe}
        else
            mv {wildcards.output_path}/{wildcards.sample}/03_align_fastq/{wildcards.sample}_fixed_mate_info.nonCG_filtered.bam {output.nonCG_filtered_pe}
        fi
        """

rule merge_and_sort_bam:
    input:
        nonCG_filtered_se="{output_path}/{sample}/03_align_fastq/{sample}_combined_nonCG_filtered.bam",
        nonCG_filtered_pe="{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info_nonCG_filtered.bam"
    output:
        merged_bam="{output_path}/{sample}/03_align_fastq/{sample}_merged_and_sorted.bam",
        bai="{output_path}/{sample}/03_align_fastq/{sample}_merged_and_sorted.bam.bai",
        flag="{output_path}/{sample}/03_align_fastq/{sample}_merged_and_sorted_flagstat.txt",
        bam_se="{output_path}/{sample}/03_align_fastq/{sample}_combined_s.bam",
        bai_se="{output_path}/{sample}/03_align_fastq/{sample}_combined_s.bam.bai",
        flag_se="{output_path}/{sample}/03_align_fastq/{sample}_combined_s_flagstat.txt",
        bam_pe="{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info_s.bam",
        bai_pe="{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info_s.bam.bai",
        flag_pe="{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info_s_flagstat.txt",
    log:
        "{output_path}/{sample}/logs/{sample}_sort_bam.log",
    conda:
        "../envs/samtools.yaml",
    threads: lambda wcs: get_cpus(1, 16),
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 200, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        """
        # Merge BAM files
        samtools merge -f -@ {threads} -O bam -o {output.merged_bam}.unsorted {input.nonCG_filtered_se} {input.nonCG_filtered_pe} &>> {log}

        # Sort the merged BAM file
        samtools sort -@ {threads} -O bam -o {output.merged_bam} {output.merged_bam}.unsorted &>> {log}

        # Index the sorted BAM file
        samtools index -@ {threads} {output.merged_bam} &>> {log}

        # Generate flagstat report
        samtools flagstat -@ {threads} {output.merged_bam} 1> {output.flag} 2>> {log}

        rm {output.merged_bam}.unsorted

        samtools sort -@ {threads} -O bam -o {output.bam_se} {input.nonCG_filtered_se} &>> {log}
        samtools index -@ {threads} {output.bam_se} &>> {log}
        samtools flagstat -@ {threads} {output.bam_se} 1> {output.flag_se} 2>> {log}

        samtools sort -@ {threads} -O bam -o {output.bam_pe} {input.nonCG_filtered_pe} &>> {log}
        samtools index -@ {threads} {output.bam_pe} &>> {log}
        samtools flagstat -@ {threads} {output.bam_pe} 1> {output.flag_pe} 2>> {log}
        """

rule picard_metrics:
    input:
        bam="{output_path}/{sample}/03_align_fastq/{sample}_merged_and_sorted.bam"
    output:
        histogram="{output_path}/{sample}/03_align_fastq/{sample}_merged_insert_size_histogram.pdf",
        insert_size="{output_path}/{sample}/03_align_fastq/{sample}_merged_insert_size_metrics.txt",
        hybrid_selection="{output_path}/{sample}/03_align_fastq/{sample}_merged_hs_metrics.txt"
    params:
        roi_bed=lambda wcs: D_sample_details[wcs.sample].get("roi_bed", ""),
        genome_path=lambda wcs: D_sample_details[wcs.sample]["genome_path"],
        genome=lambda wcs: D_sample_details[wcs.sample]["genome"],
        interval_list=lambda wcs: "{}/{}/03_align_fastq/ROI.interval_list".format(wcs.output_path, wcs.sample),
    log:
        "{output_path}/{sample}/logs/{sample}_picard_metrics.log",
    conda:
        "../envs/picard.yaml",
    threads: 1,
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 5, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        """
        picard CollectInsertSizeMetrics -I {input.bam} -O {output.insert_size} -H {output.histogram} &>> {log}
        
        if [ -n '{params.roi_bed}' ]; then
            picard BedToIntervalList -I {params.roi_bed} -O {params.interval_list} -SD {params.genome_path}/{params.genome}.dict &>> {log}
            picard CollectHsMetrics -I {input.bam} -O {output.hybrid_selection} -BI {params.interval_list} -TI {params.interval_list} &>> {log}
        else
            echo 'no regions of interest specified' > {output.hybrid_selection}
        fi
        """