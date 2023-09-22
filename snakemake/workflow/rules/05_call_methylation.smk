# Rule for the call_methylation step, which generates methylation information from deduplicated BAM files
rule call_methylation:
    input:
        "{output_path}/{sample}/04_deduplicate_bam/{sample}_sd.bam",  # Input deduplicated BAM file
    output:
        # SVG and bedGraph files for methylation in different contexts (CHH, CHG, CpG)
        "{output_path}/{sample}/05_call_methylation/{sample}_CHH_OT.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_CHH_OB.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_CHG_OT.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_CHG_OB.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_CpG_OT.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_CpG_OB.svg",
        "{output_path}/{sample}/05_call_methylation/{sample}_ROI_Conversion_CHH.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_ROI_Conversion_CHG.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_ROI_Conversion_CpG.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_sd_CpG.bedGraph",
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
    threads: lambda wcs: get_cpus(1,64),
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "call_methylation", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        "mkdir -p {wildcards.output_path}/{wildcards.sample}/05_call_methylation \n\n"
        # MethylDackel is used for methylation extraction and metrics calculation
        "MethylDackel extract --CHH --CHG --minOppositeDepth 5 --maxVariantFrac 0.2 "
        "--opref {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_ROI_Conversion "
        "-@ {threads} -l {params.ROI_path} {params.genome_fa} {input} "
        "&>> {log} \n\n"

        # Generate mbias plots for CHH context
        "MethylDackel mbias --CHH --noCpG --txt -@ {threads} {params.genome_fa} "
        "{input} {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_CHH "
        "&>> {log} \n\n"

        # Generate mbias plots for CHG context
        "MethylDackel mbias --CHG --noCpG --txt -@ {threads} {params.genome_fa} "
        "{input} {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_CHG "
        "&>> {log} \n\n"

        # Generate mbias plots for CpG context
        "MethylDackel mbias --txt -@ {threads} {params.genome_fa} {input} "
        "{wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_CpG "
        "&>> {log} \n\n"

        # Generate overall methylation extract without context separation
        "MethylDackel extract -@ {threads} --mergeContext {params.genome_fa} {input} "
        "-o {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd &>> {log}"
