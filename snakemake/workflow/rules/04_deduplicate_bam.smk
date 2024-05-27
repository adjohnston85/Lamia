# rule split_deaminated_genomes:
#     input:
#         bam="{output_path}/{sample}/03_align_fastq/{sample}_s.bam"
#     output:
#         ct_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_CT_genome.bam",
#         ga_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_GA_genome.bam",
#         ct_bam_index="{output_path}/{sample}/04_deduplicate_bam/{sample}_CT_genome.bam.bai",
#         ga_bam_index="{output_path}/{sample}/04_deduplicate_bam/{sample}_GA_genome.bam.bai"
#     log:
#         "{output_path}/{sample}/logs/{sample}_split_genomes.log"
#     conda:
#         "../envs/gencore.yaml",
#     threads:
#         1
#     resources:
#         # Dynamic resource allocation and user/account specific parameters
#         time_min=lambda wcs, input, threads: get_time_min(wcs, input, "split_deaminated_genomes", threads),
#         cpus=lambda wcs, threads: threads,
#         mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
#         account=lambda wcs: D_sample_details[wcs.sample]['account'],
#         email=lambda wcs: D_sample_details[wcs.sample]['email'],
#         partition=""
#     shell:
#         'mkdir -p {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam \n\n' 
#         'samtools view -h {input.bam} | grep -e "XG:Z:CT" -e "@" '
#         '| samtools view -bS - 1> {output.ct_bam} 2>> {log} \n'
#         'samtools index {output.ct_bam} 2>> {log} \n\n'
#         'samtools view -h {input.bam} | grep -e "XG:Z:GA" -e "@" '
#         '| samtools view -bS - 1> {output.ga_bam} 2>> {log} \n'
#         'samtools index {output.ga_bam} 2>> {log} \n\n'

# rule revert_deamination:
#     input:
#         ct_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_CT_genome.bam",
#         ga_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_GA_genome.bam"
#     output:
#         ct_reverted_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_CT_genome_reverted.bam",
#         ga_reverted_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_GA_genome_reverted.bam",
#         ct_reverted_bam_index="{output_path}/{sample}/04_deduplicate_bam/{sample}_CT_genome_reverted.bam.bai",
#         ga_reverted_bam_index="{output_path}/{sample}/04_deduplicate_bam/{sample}_GA_genome_reverted.bam.bai",
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
#         'python scripts/deamination_reverter.py --ct_bam_path {input.ct_bam} --ga_bam_path {input.ga_bam} --num_threads {threads} '
#         '--cpg_bed_path_prefix data/CpGs_hg38_ --make_cpg_bed --ref_genome_path {params.genome_path}/{params.genome}.fa '
#         '--chromosomes core_human &>> {log} \n\n'

#         'samtools index {output.ct_reverted_bam} 2>> {log} \n'
#         'samtools index {output.ga_reverted_bam} 2>> {log} \n'

# rule merge_sort_reverted:
#     input:
#         ct_reverted_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_CT_genome_reverted.bam",
#         ga_reverted_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_GA_genome_reverted.bam"
#     output:
#         sorted_reverted_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_reverted.bam",
#         sorted_reverted_bam_index="{output_path}/{sample}/04_deduplicate_bam/{sample}_reverted.bam.bai"
#     log:
#         "{output_path}/{sample}/logs/{sample}_merge_sort_reverted.log"
#     conda:
#         "../envs/gencore.yaml",
#     threads:
#         lambda wcs: get_cpus(1,64),  # Dynamic CPU allocation. Max limit is 64.
#     resources:
#         # Dynamic resource allocation and user/account specific parameters
#         time_min=lambda wcs, input, threads: get_time_min(wcs, input, "merge_sort_reverted", threads),
#         cpus=lambda wcs, threads: threads,
#         mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
#         account=lambda wcs: D_sample_details[wcs.sample]['account'],
#         email=lambda wcs: D_sample_details[wcs.sample]['email'],
#         partition=""
#     shell:
#         'samtools merge -@ {threads} {output.sorted_reverted_bam}.unsorted {input.ct_reverted_bam} {input.ga_reverted_bam} 2>> {log} \n'
#         'samtools sort -m 2G -@ {threads} -o {output.sorted_reverted_bam} {output.sorted_reverted_bam}.unsorted 2>> {log} \n'
#         'samtools index {output.sorted_reverted_bam} 2>> {log} \n\n'

rule correct_umis:
    input:
        bam="{output_path}/{sample}/03_align_fastq/{sample}_s.bam"
    output:
        corrected_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_corrected.bam",
        metrics="{output_path}/{sample}/04_deduplicate_bam/{sample}_metrics.txt",
        rejected="{output_path}/{sample}/04_deduplicate_bam/{sample}_rejected.bam"
    params:
        umi_list = lambda wcs: D_sample_details[wcs.sample]["umi_list"],
    log:
        "{output_path}/{sample}/logs/{sample}_correct_umis.log"
    conda:
        "../envs/fgbio.yaml"
    threads: 1
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "correct_umis", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        'fgbio CorrectUmis '
            '-i {input.bam} '
            '-o {output.corrected_bam} '
            '-M {output.metrics} '
            '-r {output.rejected} '
            '-t RX -u {params.umi_list} '
            '--max-mismatches=3 '
            '--min-distance=1 '

rule group_reads_by_umi:
    input:
        corrected_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_corrected.bam"
    output:
        grouped_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_grouped.bam"
    log:
        "{output_path}/{sample}/logs/{sample}_group_reads_by_umi.log"
    conda:
        "../envs/fgbio.yaml"
    threads: 1
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "group_reads_by_umi", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        'fgbio GroupReadsByUmi '
            '--input={input.corrected_bam} '
            '--output={output.grouped_bam} '
            '--strategy=paired '
            '--edits=0 '
            '--min-map-q=20'

rule call_molecular_consensus_reads:
    input:
        grouped_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_grouped.bam"
    output:
        consensus_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus.bam",
    log:
        "{output_path}/{sample}/logs/{sample}_call_consensus_reads.log"
    conda:
        "../envs/fgbio.yaml"
    threads: 1
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "call_molecular_consensus_reads", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        """
        fgbio CallMolecularConsensusReads \
            --input={input.grouped_bam} \
            --output={output.consensus_bam}.preheader \
            --min-reads=1 &> {log}

        samtools view -H {output.consensus_bam}.preheader | \
            awk 'BEGIN{{OFS="\\t"}} /^@RG/{{if($0 !~ /SM:/) $0=$0 "\\tSM:SampleName"; else $0=gensub(/(SM:)[^\\t]+/, "\\1SampleName", "g")}} {{print $0}}' | \
            samtools reheader - {output.consensus_bam}.preheader > {output.consensus_bam}

        rm {output.consensus_bam}.preheader
        """

# rule reject_mixed_methylation:
#     input:
#         bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_bismark_bt2_pe.bam",
#     output:
#         phased_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_bismark_bt2_pe_phased.bam",
#         phased_bam_index="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_bismark_bt2_pe_phased.bam.bai",
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
#         'python scripts/CpG_phaser.py --bam_path {input.bam} --num_threads {threads} '
#         '--cpg_bed_path_prefix data/CpGs_hg38_ --make_cpg_bed --ref_genome_path {params.genome_path}/{params.genome}.fa '
#         '--chromosomes core_human &>> {log} \n\n'

#         'samtools index {output.phased_bam} 2>> {log} \n'

rule filter_consensus_reads:
    input:
        bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus.bam"
    output:
        filtered_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_filtered.bam"
    params:
        genome_path = lambda wcs: D_sample_details[wcs.sample]["genome_path"],
        genome = lambda wcs: D_sample_details[wcs.sample]["genome"],
    log:
        "{output_path}/{sample}/logs/{sample}_filter_consensus_reads.log"
    conda:
        "../envs/fgbio.yaml"  # Ensure this environment includes fgbio
    threads: 1
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "filter_consensus_reads", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        'fgbio FilterConsensusReads '
            '--input={input.bam} '
            '--output={output.filtered_bam} '
            '--ref={params.genome_path}/{params.genome}.fa '
            '--min-reads=1 '
            '--max-read-error-rate=0.05 '
            '--max-base-error-rate=0.1 '
            '--min-base-quality=40 '
            '--max-no-call-fraction=0.1 '
            '&> {log}'

rule convert_bam_to_fastq_pair:
    input:
        bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus.bam"
    output:
        fastq1="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1.fastq",
        fastq2="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r2.fastq"
    log:
        "{output_path}/{sample}/logs/{sample}_convert_bam_to_fastq_pair.log"
    conda:
        "../envs/picard.yaml"
    threads: 1
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "convert_bam_to_fastq_pair", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        """
        picard SamToFastq \
            -I {input.bam} \
            -F {output.fastq1} \
            -F2 {output.fastq2} \
            -VALIDATION_STRINGENCY LENIENT \
            -TMP_DIR {wildcards.output_path}/{wildcards.sample}/tmp &> {log}
        """

rule bismark_realign:
    input:
        r1="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1.fastq",
        r2="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r2.fastq",
        genome_path = lambda wcs: D_sample_details[wcs.sample]["genome_path"],
    output:
        report="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_bismark_bt2_PE_report.txt",
        nucleotide_stats="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_bismark_bt2_pe.nucleotide_stats.txt",
        bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_bismark_bt2_pe.bam",
    params:
        parallel = lambda wcs, threads: get_parallel(wcs, 5, threads),
        non_directional = lambda wcs: "--non_directional " if D_sample_details[wcs.sample]["non_directional"] else "",
        maxins = lambda wcs: D_sample_details[wcs.sample]["maxins"],
    log:
        "{output_path}/{sample}/logs/{sample}_bismark_align.log",
    conda:
        "../envs/bismark.yaml",
    threads: lambda wcs: get_cpus(16,64),
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "bismark_align", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 4096),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        'bismark '
            '--genome_folder {input.genome_path} '
            '--multicore {params.parallel} '
            '--bowtie2 {params.non_directional}'
            '--temp_dir {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam '
            '-1 {input.r1} '
            '-2 {input.r2} '
            '--maxins {params.maxins} '
            '--nucleotide_coverage '
            '--output_dir {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam '
            '1>> /dev/null 2>> {log}'

rule clip_bam:
    input:
        bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_r1_bismark_bt2_pe.bam",
    output:
        clipped_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_clipped.bam"
    params:
        genome_path = lambda wcs: D_sample_details[wcs.sample]["genome_path"],
        genome = lambda wcs: D_sample_details[wcs.sample]["genome"],
    log:
        "{output_path}/{sample}/logs/{sample}_clip_bam.log"
    conda:
        "../envs/fgbio.yaml"  # Ensure this environment includes fgbio
    threads: 1
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "clip_bam", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 32768),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        'fgbio ClipBam '
            '--input={input.bam} '
            '--output={output.clipped_bam} '
            '--ref={params.genome_path}/{params.genome}.fa '
            '--clip-overlapping-reads=true '
            '&> {log}'


# rule deduplicate_bam:
#     input:
#         sorted_reverted_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_reverted_s.bam"
#     output:
#         dedup_reverted_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_reverted_sd.bam",
#         dedup_reverted_bam_index="{output_path}/{sample}/04_deduplicate_bam/{sample}_reverted_sd.bam.bai",
#         json="{output_path}/{sample}/04_deduplicate_bam/{sample}_gencore.json",
#         html="{output_path}/{sample}/04_deduplicate_bam/{sample}_gencore.html",
#         flag="{output_path}/{sample}/04_deduplicate_bam/{sample}_sd_flagstat.txt"
#     params:
#         # Conditionally add BED file for region of interest if available
#         b=lambda wcs: '-b ' + D_sample_details[wcs.sample]["roi_bed"] + ' '
#         if "roi_bed" in D_sample_details[wcs.sample] else "",
#         # UMI prefix dynamically obtained based on the sample
#         umi_prefix=lambda wcs: get_umi_prefix(wcs.sample),
#         genome_path = lambda wcs: D_sample_details[wcs.sample]["genome_path"],
#         genome = lambda wcs: D_sample_details[wcs.sample]["genome"],
#     log:
#         "{output_path}/{sample}/logs/{sample}_deduplication.log"
#     conda:
#         "../envs/gencore.yaml",
#     threads:
#         lambda wcs: get_cpus(1,64),  # Dynamic CPU allocation. Max limit is 64.
#     resources:
#         # Dynamic resource allocation and user/account specific parameters
#         time_min=lambda wcs, input, threads: get_time_min(wcs, input, "deduplicate_bam", threads),
#         cpus=lambda wcs, threads: threads,
#         mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 4096),
#         account=lambda wcs: D_sample_details[wcs.sample]['account'],
#         email=lambda wcs: D_sample_details[wcs.sample]['email'],
#         partition=""
#     params:
#         umi_prefix=lambda wcs: get_umi_prefix(wcs.sample),
#         genome_path=lambda wcs: D_sample_details[wcs.sample]["genome_path"],
#         genome=lambda wcs: D_sample_details[wcs.sample]["genome"],
#         b=lambda wcs: '-b ' + D_sample_details[wcs.sample]["roi_bed"] + ' ' if "roi_bed" in D_sample_details[wcs.sample] else ""
#     shell:
#         'gencore {params.umi_prefix}-s 1 -d 2 -i {input.sorted_reverted_bam} -o {output.dedup_reverted_bam} '
#         '-r {params.genome_path}/{params.genome}.fa {params.b} -j {output.json} -h {output.html} &>> {log} \n\n'

#         'samtools index {output.dedup_reverted_bam} 2>> {log} \n'
#         'samtools flagstat {output.dedup_reverted_bam} 1> {output.flag} 2>> {log} \n\n'
