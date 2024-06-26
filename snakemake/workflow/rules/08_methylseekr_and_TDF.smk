# Rule to perform methylseekr analysis and convert methylation calls to TDF format.
rule methylseekr_and_TDF:
    input:
        # Input bedGraph and conversion statistics files
        bedgraph="{output_path}/{sample}/05_call_methylation/{sample}_sd_CpG.bedGraph",
        cov="{output_path}/{sample}/07_calculate_statistics/{sample}_sd_conversion_and_coverage.txt",
    output:
        # Output files for methylseekr results and TDF conversion
        "{output_path}/{sample}/08_methylseekr_and_TDF/{sample}_PMD.bed",
        "{output_path}/{sample}/08_methylseekr_and_TDF/{sample}_UMRLMR.bed",
        "{output_path}/{sample}/08_methylseekr_and_TDF/{sample}_wPMD_UMRLMR.bed",
        "{output_path}/{sample}/08_methylseekr_and_TDF/{sample}_sd_CpG.tdf",
    params:
        # Define genome and genome directory parameters
        g=lambda wcs: D_sample_details[wcs.sample]["genome"],
        e=lambda wcs: D_sample_details[wcs.sample]["genome_dir"],
    log:
        "{output_path}/{sample}/logs/{sample}_methylseekr_and_TDF.log",
    conda:
       "../envs/methylseekr.yaml"
    threads: lambda wildcards: get_cpus(1,64)
    resources:
        time_min=lambda wildcards, input, threads: get_time_min(wildcards, input, "methylseekr_and_TDF", threads),
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wcs, threads: get_mem_mb(wcs, threads, 2048),
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        # Shell command to execute R script for methylseekr analysis and TDF conversion
        'Rscript {workflow.basedir}/scripts/CallMethylseekrRegions_and_convertMethCallsToTdf.R '
        '-g {params.g} -i {input.bedgraph} -t {workflow.basedir}/data/TissueToEmbryoMap.csv '
        '-e {params.e} -p {resources.cpus} -o {wildcards.output_path}/{wildcards.sample}/08_methylseekr_and_TDF '
        '&>> {log}'

