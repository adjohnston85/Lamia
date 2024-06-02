# Rule for the call_methylation step, which generates methylation information from both sorted and deduplicated BAM files
rule call_methylation:
    input:
        s_bam="{output_path}/{sample}/03_align_fastq/{sample}_merged_and_sorted.bam",  # Input sorted BAM file
        sd_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe.bam",  # Input deduplicated BAM file
        sdd_bam=lambda wildcards: (
            "{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe_duplex.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
            if "umi_len" in D_sample_details[wildcards.sample]
            else "None"
        )
    output:
        # SVG and bedGraph files for methylation in different contexts (CHH, CHG, CpG)
        # Outputs for sorted BAM
        "{output_path}/{sample}/05_call_methylation/{sample}_s_ROI_Conversion_CHH.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_s_ROI_Conversion_CHG.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_s_ROI_Conversion_CpG.bedGraph",
        # Outputs for deduplicated BAM
        "{output_path}/{sample}/05_call_methylation/{sample}_CHH_OT.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_CHH_OB.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_CHG_OT.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_CHG_OB.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_CpG_OT.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_CpG_OB.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_sd_ROI_Conversion_CHH.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_sd_ROI_Conversion_CHG.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_sd_ROI_Conversion_CpG.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_sd_CpG.bedGraph",
        CHH="{output_path}/{sample}/05_call_methylation/{sample}_CHH_mbias.txt",
        CHG="{output_path}/{sample}/05_call_methylation/{sample}_CHG_mbias.txt",
        CpG="{output_path}/{sample}/05_call_methylation/{sample}_CpG_mbias.txt",
    params:
        # Use regions of interest (ROI) if available, otherwise default to pre-defined CHH ROI
        ROI_path=lambda wcs: D_sample_details[wcs.sample]["roi_bed"]
        if "roi_bed" in D_sample_details[wcs.sample]
        else D_sample_details[wcs.sample]["genome_path"] + "/CHH_ROI.bed",
        # Locate the relevant genome fasta file
        genome_fa=lambda wcs: D_sample_details[wcs.sample]["genome_path"] + '/' + D_sample_details[wcs.sample]["genome"] + ".fa",
    log:
        "{output_path}/{sample}/logs/{sample}_call_methylation.log",
    conda:
        "../envs/methyldackel.yaml",
    threads: lambda wcs: get_cpus(1, 64),
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "call_methylation", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        """
        mkdir -p {wildcards.output_path}/{wildcards.sample}/05_call_methylation

        # Process sorted BAM
        MethylDackel extract --CHH --CHG --minOppositeDepth 5 --maxVariantFrac 0.2 \
        --opref {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_ROI_Conversion \
        -@ {threads} -l {params.ROI_path} {params.genome_fa} {input.s_bam} \
        &>> {log}

        # Process deduplicated BAM
        MethylDackel extract --CHH --CHG --minOppositeDepth 5 --maxVariantFrac 0.2 \
        --opref {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_ROI_Conversion \
        -@ {threads} -l {params.ROI_path} {params.genome_fa} {input.sd_bam} \
        &>> {log}

        # Generate mbias plots for CHH context
        MethylDackel mbias --CHH --noCpG --txt -@ {threads} {params.genome_fa} \
        {input.sd_bam} {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_CHH \
        > {output.CHH} 2>> {log}

        # Generate mbias plots for CHG context
        MethylDackel mbias --CHG --noCpG --txt -@ {threads} {params.genome_fa} \
        {input.sd_bam} {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_CHG \
        > {output.CHG} 2>> {log}

        # Generate mbias plots for CpG context
        MethylDackel mbias --txt -@ {threads} {params.genome_fa} {input.sd_bam} \
        {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_CpG \
        > {output.CpG} 2>> {log}

        # Generate overall methylation extract without context separation for deduplicated BAM
        MethylDackel extract -@ {threads} --mergeContext {params.genome_fa} {input.sd_bam} \
        -o {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd &>> {log}

        # Conditionally process duplex BAM
        if [ "{input.sdd_bam}" != "None" ]; then
            MethylDackel extract -@ {threads} --mergeContext {params.genome_fa} {input.sdd_bam} \
            -o {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_duplex &>> {log}
        fi
        """

# rule call_methylation_bismark:
#     input:
#         sd_bam="{output_path}/{sample}/04_deduplicate_bam/{sample}_merged_se_pe.bam",
#     output:
#         report="{output_path}/{sample}/05_call_methylation/{sample}_merged_se_pe_CpG_report.txt.gz",
#         bedGraph="{output_path}/{sample}/05_call_methylation/{sample}_merged_se_pe_CpG.bedGraph.gz",
#         summary="{output_path}/{sample}/05_call_methylation/{sample}_methylation_summary.txt",
#     params:
#         genome_folder=lambda wcs: D_sample_details[wcs.sample]["genome_path"],
#         trim_lens=lambda wcs: D_sample_details[wcs.sample]["bam_trim_lengths"],
#     log:
#         "{output_path}/{sample}/logs/{sample}_call_methylation_bismark.log",
#     conda:
#         "../envs/bismark.yaml",
#     threads: lambda wcs: get_cpus(1,64),
#     resources:
#         time_min=lambda wcs, input, threads: get_time_min(wcs, input, "call_methylation", threads),
#         cpus=lambda wcs, threads: threads,
#         mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
#         account=lambda wcs: D_sample_details[wcs.sample]['account'],
#         email=lambda wcs: D_sample_details[wcs.sample]['email'],
#         partition=""
#     shell:
#         """
#         mkdir -p {wildcards.output_path}/{wildcards.sample}/05_call_methylation

#         # Process duplex BAM with Bismark Methylation Extractor
#         bismark_methylation_extractor --comprehensive --bedGraph --gzip {params.trim_lens}\
#         --output {wildcards.output_path}/{wildcards.sample}/05_call_methylation \
#         --genome_folder {params.genome_folder} {input.duplex_bam} &>> {log}

#         # Summarize methylation extraction results
#         bismark2summary {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_consensus_clipped_duplex_CpG_report.txt.gz &>> {log}
#         """