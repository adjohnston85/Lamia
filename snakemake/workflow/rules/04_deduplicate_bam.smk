rule picard_deduplication:
    input:
        bam_se="{output_path}/{sample}/03_align_fastq/{sample}_combined_s.bam",
        bam_pe="{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info_s.bam",
    output:
        dedup_bam_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_dedup_se.bam",
        dedup_bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_dedup_pe.bam",
        dedup_bam_merged="{output_path}/{sample}/04_deduplicate_bam/{sample}_dedup_merged_and_sorted_se_pe.bam",
        dedup_bam_merged_bai="{output_path}/{sample}/04_deduplicate_bam/{sample}_dedup_merged_and_sorted_se_pe.bam.bai",
    log:
        "{output_path}/{sample}/logs/{sample}_picard_deduplication.log",
    conda:
        "../envs/picard.yaml",
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
        picard MarkDuplicates I={input.bam_se} O={output.dedup_bam_se} M={output.dedup_bam_se}.metrics.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true 1>> {log} 2>&1

        picard MarkDuplicates I={input.bam_pe} O={output.dedup_bam_pe} M={output.dedup_bam_pe}.metrics.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true 1>> {log} 2>&1

        samtools merge -@ {threads} -f {output.dedup_bam_merged}.unsorted {output.dedup_bam_se} {output.dedup_bam_pe} &>> {log}

        # Sort the merged BAM file
        samtools sort -@ {threads} -O bam -o {output.dedup_bam_merged} {output.dedup_bam_merged}.unsorted &>> {log}

        # Index the sorted BAM file
        samtools index -@ {threads} {output.dedup_bam_merged} &>> {log}

        rm {output.dedup_bam_merged}.unsorted
        """

rule correct_umis:
    input:
        bam_se="{output_path}/{sample}/03_align_fastq/{sample}_combined_s.bam",
        bam_pe="{output_path}/{sample}/03_align_fastq/{sample}_fixed_mate_info_s.bam"
    output:
        corrected_bam_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_corrected_se.bam",
        corrected_bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_corrected_pe.bam",
        metrics_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_metrics_se.txt",
        metrics_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_metrics_pe.txt",
        rejected_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_rejected_se.bam",
        rejected_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_rejected_pe.bam"
    params:
        umi_list=lambda wcs: D_sample_details[wcs.sample]["umi_list"],
    log:
        "{output_path}/{sample}/logs/{sample}_correct_umis.log",
    conda:
        "../envs/fgbio.yaml",
    threads: 1,
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 100, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        """
        fgbio CorrectUmis -i {input.bam_se} -o {output.corrected_bam_se} -M {output.metrics_se} -r {output.rejected_se} -t RX -u {params.umi_list} --max-mismatches=3 --min-distance=1

        fgbio CorrectUmis -i {input.bam_pe} -o {output.corrected_bam_pe} -M {output.metrics_pe} -r {output.rejected_pe} -t RX -u {params.umi_list} --max-mismatches=3 --min-distance=1
        """

rule group_reads_by_umi:
    input:
        corrected_bam_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_corrected_se.bam",
        corrected_bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_corrected_pe.bam"
    output:
        grouped_bam_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_grouped_se.bam",
        grouped_bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_grouped_pe.bam"
    log:
        "{output_path}/{sample}/logs/{sample}_group_reads_by_umi.log",
    conda:
        "../envs/fgbio.yaml",
    threads: 1,
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 4, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        """
        fgbio GroupReadsByUmi --input={input.corrected_bam_se} --output={output.grouped_bam_se} --strategy=adjacency --edits=1 --min-map-q=20

        fgbio GroupReadsByUmi --input={input.corrected_bam_pe} --output={output.grouped_bam_pe} --strategy=paired --edits=0 --min-map-q=20
        """

rule call_molecular_consensus_reads:
    input:
        grouped_bam_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_grouped_se.bam",
        grouped_bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_grouped_pe.bam"
    output:
        consensus_bam_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_se.bam",
        consensus_bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_pe.bam"
    params:
        tmp_dir="{output_path}/{sample}/tmp"
    log:
        "{output_path}/{sample}/logs/{sample}_call_consensus_reads.log",
    conda:
        "../envs/fgbio.yaml",
    threads: 1,
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 10, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb= 65536,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        """
        fgbio -Xmx512g -Xms4G -Djava.io.tmpdir={params.tmp_dir} CallMolecularConsensusReads --input={input.grouped_bam_se} --output={output.consensus_bam_se}.preheader --min-reads=1 &> {log}

        samtools view -H {output.consensus_bam_se}.preheader | \
            awk 'BEGIN{{OFS="\\t"}} /^@RG/{{if($0 !~ /SM:/) $0=$0 "\\tSM:SampleName"; else $0=gensub(/(SM:)[^\\t]+/, "\\1SampleName", "g")}} {{print $0}}' | \
            samtools reheader - {output.consensus_bam_se}.preheader > {output.consensus_bam_se}

        rm {output.consensus_bam_se}.preheader

        fgbio -Xmx512g -Xms4G -Djava.io.tmpdir={params.tmp_dir} CallMolecularConsensusReads --input={input.grouped_bam_pe} --output={output.consensus_bam_pe}.preheader --min-reads=1 &> {log}

        samtools view -H {output.consensus_bam_pe}.preheader | \
            awk 'BEGIN{{OFS="\\t"}} /^@RG/{{if($0 !~ /SM:/) $0=$0 "\\tSM:SampleName"; else $0=gensub(/(SM:)[^\\t]+/, "\\1SampleName", "g")}} {{print $0}}' | \
            samtools reheader - {output.consensus_bam_pe}.preheader > {output.consensus_bam_pe}

        rm {output.consensus_bam_pe}.preheader
        """

rule filter_consensus_reads:
    input:
        bam_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_se.bam",
        bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_pe.bam"
    output:
        filtered_bam_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_filtered_se.bam",
        filtered_bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_filtered_pe.bam"
    params:
        genome_path=lambda wcs: D_sample_details[wcs.sample]["genome_path"],
        genome=lambda wcs: D_sample_details[wcs.sample]["genome"],
        tmp_dir="{output_path}/{sample}/tmp"
    log:
        "{output_path}/{sample}/logs/{sample}_filter_consensus_reads.log",
    conda:
        "../envs/fgbio.yaml",
    threads: 1,
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 1, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        """
        fgbio -Xmx512g -Xms4G -Djava.io.tmpdir={params.tmp_dir} FilterConsensusReads --input={input.bam_se} --output={output.filtered_bam_se} --ref={params.genome_path}/{params.genome}.fa \
            --min-reads=3 --max-read-error-rate=0.05 --max-base-error-rate=0.1 --min-base-quality=40 --max-no-call-fraction=0.1 &> {log}

        fgbio -Xmx512g -Xms4G -Djava.io.tmpdir={params.tmp_dir} FilterConsensusReads --input={input.bam_pe} --output={output.filtered_bam_pe} --ref={params.genome_path}/{params.genome}.fa \
            --min-reads=3 --max-read-error-rate=0.05 --max-base-error-rate=0.1 --min-base-quality=40 --max-no-call-fraction=0.1 &> {log}
        """

rule convert_bam_to_fastq:
    input:
        bam_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_filtered_se.bam",
        bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_filtered_pe.bam"
    output:
        fastq_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_se.fastq",
        fastq1_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_pe.fastq",
        fastq2_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r2_pe.fastq"
    log:
        "{output_path}/{sample}/logs/{sample}_convert_bam_to_fastq.log"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 1, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        """
        # Convert BAM to FASTQ for single-end reads with UMI
        samtools view -h {input.bam_se} | \
        awk 'BEGIN {{FS=OFS="\t"}} 
            /^@/ {{next}} 
            {{umi=""; 
            for(i=12;i<=NF;i++) 
                if ($i ~ /^RX:Z:/) umi=$i; 
            gsub(/RX:Z:/,"",umi); 
            print "@"$1":"umi"\\n"$10"\\n+\\n"$11}}' > {output.fastq_se}

        # Convert BAM to FASTQ for paired-end reads (read 1) with UMI
        samtools view -h -f 64 {input.bam_pe} | \
        awk 'BEGIN {{FS=OFS="\t"}} 
            /^@/ {{next}} 
            {{umi=""; 
            for(i=12;i<=NF;i++) 
                if ($i ~ /^RX:Z:/) umi=$i; 
            gsub(/RX:Z:/,"",umi); 
            print "@"$1":"umi"\\n"$10"\\n+\\n"$11}}' > {output.fastq1_pe}

        # Convert BAM to FASTQ for paired-end reads (read 2) with UMI
        samtools view -h -f 128 {input.bam_pe} | \
        awk 'BEGIN {{FS=OFS="\t"}} 
            /^@/ {{next}} 
            {{umi=""; 
            for(i=12;i<=NF;i++) 
                if ($i ~ /^RX:Z:/) umi=$i; 
            gsub(/RX:Z:/,"",umi); 
            print "@"$1":"umi"\\n"$10"\\n+\\n"$11}}' > {output.fastq2_pe}
        """

rule bismark_realign:
    input:
        r1_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_se.fastq",
        r1_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_pe.fastq",
        r2_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r2_pe.fastq",
        genome_path=lambda wcs: D_sample_details[wcs.sample]["genome_path"],
    output:
        report_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_se_bismark_bt2_SE_report.txt",
        nucleotide_stats_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_se_bismark_bt2.nucleotide_stats.txt",
        bam_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_se_bismark_bt2.bam", 
        report_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_pe_bismark_bt2_PE_report.txt",
        nucleotide_stats_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_pe_bismark_bt2_pe.nucleotide_stats.txt",
        bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_pe_bismark_bt2_pe.bam",
    params:
        parallel=lambda wcs, threads: get_parallel(wcs, 5, threads),
        non_directional=lambda wcs: "--non_directional " if D_sample_details[wcs.sample]["non_directional"] else "",
        maxins=lambda wcs: D_sample_details[wcs.sample]["maxins"],
    log:
        "{output_path}/{sample}/logs/{sample}_bismark_realign.log",
    conda:
        "../envs/bismark.yaml",
    threads: lambda wcs: get_cpus(16, 64),
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 500, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 4096),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        """
        bismark --genome_folder {input.genome_path} --multicore {params.parallel} --bowtie2 {params.non_directional} --temp_dir {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam \
            {input.r1_se} --nucleotide_coverage --output_dir {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam 1>> /dev/null 2>> {log}

        bismark --genome_folder {input.genome_path} --multicore {params.parallel} --bowtie2 {params.non_directional} --temp_dir {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam \
            -1 {input.r1_pe} -2 {input.r2_pe} --maxins {params.maxins} --nucleotide_coverage --output_dir {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam 1>> /dev/null 2>> {log}
        """

rule clip_bam:
    input:
        bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_pe_bismark_bt2_pe.bam"
    output:
        clipped_bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_clipped_pe.bam"
    params:
        genome_path=lambda wcs: D_sample_details[wcs.sample]["genome_path"],
        genome=lambda wcs: D_sample_details[wcs.sample]["genome"],
    log:
        "{output_path}/{sample}/logs/{sample}_clip_bam.log",
    conda:
        "../envs/fgbio.yaml",
    threads: 1,
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 1, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        """
        fgbio ClipBam --input={input.bam_pe} --output={output.clipped_bam_pe} --ref={params.genome_path}/{params.genome}.fa --clip-overlapping-reads=true &> {log}
        """

rule merge_se_pe_bams:
    input:
        bam_se="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_se_bismark_bt2.bam",
        bam_pe="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_clipped_pe.bam",
    output:
        merged_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe.bam",
        merged_bai="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe.bam.bai",
    log:
        "{output_path}/{sample}/logs/{sample}_merge_se_pe_bams.log",
    conda:
        "../envs/samtools.yaml",
    threads: 1,
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 1, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="",
    shell:
        """
        samtools merge -@ {threads} -f {output.merged_bam}.unsorted {input.bam_se} {input.bam_pe} &> {log}

        # Sort the merged BAM file
        samtools sort -@ {threads} -O bam -o {output.merged_bam} {output.merged_bam}.unsorted &>> {log}

        # Index the sorted BAM file
        samtools index -@ {threads} {output.merged_bam} &>> {log}

        rm {output.merged_bam}.unsorted
        """

rule duplex_methylation_processing:
    input:
        merged_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe.bam"
    output:
        duplex_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe_duplex.bam",
        duplex_bai="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe_duplex.bam.bai",
    params:
        genome_path=lambda wcs: D_sample_details[wcs.sample]["genome_path"],
        genome=lambda wcs: D_sample_details[wcs.sample]["genome"]
    log:
        "{output_path}/{sample}/logs/{sample}_duplex_methylation_processing.log"
    conda:
        "../envs/gencore.yaml"
    threads: get_cpus(1,64),
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, 1, threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        """
        python scripts/duplex_methylation_processor.py --bam_path {input.merged_bam} --num_threads {threads} --ref_genome_path {params.genome_path}/{params.genome}.fa &> {log}
        
        # Sort the duplex BAM file
        samtools sort -@ {threads} -O bam -o {output.duplex_bam}.sorted {output.duplex_bam} &>> {log}

        mv {output.duplex_bam}.sorted {output.duplex_bam}
        
        # Index the sorted duplex BAM file
        samtools index -@ {threads} {output.duplex_bam} &>> {log}
        """

# rule reject_mixed_methylation:
#     input:
#         clipped_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_clipped.bam"
#     output:
#         phased_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_clipped_phased.bam",
#     params:
#         genome_path = lambda wcs: D_sample_details[wcs.sample]["genome_path"],
#         genome = lambda wcs: D_sample_details[wcs.sample]["genome"],
#     log:
#         "{output_path}/{sample}/logs/{sample}_deamination_reversion.log"
#     conda:
#         "../envs/gencore.yaml",
#     threads:
#         lambda wcs: get_cpus(1,64),  # Dynamic CPU allocation. Max limit is 64.
#     resources:
#         # Dynamic resource allocation and user/account specific parameters
#         time_min=lambda wcs, input, threads: get_time_min(wcs, input, "revert_deamination", threads) + 15,
#         cpus=lambda wcs, threads: threads,
#         mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 4096),
#         account=lambda wcs: D_sample_details[wcs.sample]['account'],
#         email=lambda wcs: D_sample_details[wcs.sample]['email'],
#         partition=""
#     shell:
#         'python scripts/CpG_phaser.py --bam_path {input.clipped_bam} --num_threads {threads} '
#         '--cpg_bed_path_prefix data/CpGs_hg38_ --make_cpg_bed --ref_genome_path {params.genome_path}/{params.genome}.fa '
#         '--chromosomes core_human &>> {log} \n\n'