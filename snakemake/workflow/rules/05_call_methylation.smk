rule call_methylation:
    input:
        "{output_path}/{sample}/04_deduplicate_bam/{sample}_sd.bam",
    output:
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
        ROI_path=lambda wcs: D_sample_details[wcs.sample]["roi_bed"] \
        if "roi_bed" in D_sample_details[wcs.sample] \
        else D_sample_details[wcs.sample]["genome_path"] + "/CHH_ROI.bed",
        genome_fa=lambda wcs: D_sample_details[wcs.sample]["genome_path"] + \
        '/' + D_sample_details[wcs.sample]["genome"] + ".fa"
    log:
        "{output_path}/{sample}/logs/{sample}_call_methylation.log",
    conda:
        "../envs/methyldackel.yaml"
    threads: lambda wcs: get_cpus(1,64),
    resources:
        time_min=lambda wcs, input, threads: get_time_min(wcs, input, "call_methylation", threads),
        cpus=lambda wcs, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        "MethylDackel extract --CHH --CHG --minOppositeDepth 5 --maxVariantFrac 0.2 "
        "--opref {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_ROI_Conversion "
        "-@ {threads} -l {params.ROI_path} {params.genome_fa} {input} "
        "&>> {log} \n\n"

        "MethylDackel mbias --CHH --noCpG --txt -@ {threads} {params.genome_fa} "
        "{input} {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_CHH "
        "&>> {log} \n\n"

        "MethylDackel mbias --CHG --noCpG --txt -@ {threads} {params.genome_fa} "
        "{input} {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_CHG "
        "&>> {log} \n\n"

        "MethylDackel mbias --txt -@ {threads} {params.genome_fa} {input} "
        "{wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_CpG "
        "&>> {log} \n\n"

        "MethylDackel extract -@ {threads} --mergeContext {params.genome_fa} {input} "
        "-o {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd &>> {log}"
