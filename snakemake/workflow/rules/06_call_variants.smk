rule mask_converted_bases:
    input:
        "{output_path}/{sample}/04_deduplicate_bam/{sample}_sd.bam",
    output:
        calmd_bam="{output_path}/{sample}/06_call_variants/{sample}_sd_calmd.bam",
        calmd_bai="{output_path}/{sample}/06_call_variants/{sample}_sd_calmd.bam.bai",
        masked_bam="{output_path}/{sample}/06_call_variants/{sample}_sd_calmd_masked.bam",
        masked_bai="{output_path}/{sample}/06_call_variants/{sample}_sd_calmd_masked.bam.bai",
    params:
        genome_fa=lambda wcs: D_sample_details[wcs.sample]["genome_path"] +
            '/' + D_sample_details[wcs.sample]["genome"] + '.fa',
        parallel = lambda wcs, threads: threads - 1,
    log:
        "{output_path}/{sample}/logs/{sample}_mask_converted_bases.log",
    conda:
        "../envs/revelio.yaml"
    threads: lambda wildcards: get_cpus(1,64)
    resources:
        time_min=lambda wildcards, input, threads: get_time_min(wildcards, input, "mask_converted_bases", threads),
        cpus=lambda wildcards, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        "samtools calmd -b {input} {params.genome_fa} -@ {params.parallel} "
        "1> {output.calmd_bam} 2>> /dev/null \n\n"

        "samtools index -@ {params.parallel} {output.calmd_bam} &>> {log} \n\n"

        "python3 {workflow.basedir}/scripts/revelio.py -t {wildcards.output_path}/{wildcards.sample} "
        "-T {params.parallel} -Q {output.calmd_bam} {output.masked_bam} 2>> {log} \n\n"

        "samtools index -@ {params.parallel} {output.masked_bam} &>> {log}"


rule call_variants:
    input:
        masked_bam=lambda wcs: os.path.join(
            D_sample_details[wcs.sample]["output_path"],
            wcs.sample,
            "06_call_variants",
            wcs.sample + "_sd_calmd_masked.bam",
        ),
        index=lambda wcs: os.path.join(D_sample_details[wcs.sample]["genome_path"], D_sample_details[wcs.sample]["genome"] + ".fa.fai"),
    output:
        "{output_path}/{sample}/06_call_variants/{sample}_sd.vcf",
    params:
        ref=lambda wcs: os.path.join(D_sample_details[wcs.sample]["genome_path"],
            D_sample_details[wcs.sample]["genome"] + ".fa"),
        prefix=lambda wcs: os.path.join(wcs.output_path, wcs.sample,"06_call_variants/beds", D_sample_details[wcs.sample]["genome"]),
        beds=lambda wcs: calculate_vcf_files(wcs,"beds",""),
        vcfs=lambda wcs: calculate_vcf_files(wcs,"vcfs",".vcf"),        
        num_chroms=lambda wcs: len(D_sample_details[wcs.sample]["nchunks"]) - 1,
        chroms=lambda wcs: [chrom for chrom in D_sample_details[wcs.sample]["nchunks"]],
        chunks=lambda wcs: [D_sample_details[wcs.sample]["nchunks"][chrom] for chrom in D_sample_details[wcs.sample]["nchunks"]],
    log:
        "{output_path}/{sample}/logs/{sample}_call_variants.log",
    conda:
        "../envs/freebayes.yaml"
    threads: 53
    resources:
        time_min=lambda wildcards, input, threads: get_time_min(wildcards, input, "call_variants", threads),
        cpus=lambda wildcards, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        """
        chroms_list='{params.chroms}'
        declare -a chroms_array=($chroms_list)

        chunks_list='{params.chunks}'
        declare -a chunks_array=($chunks_list)

        for index in $(seq 0 {params.num_chroms}); do
            fasta_generate_regions.py --chromosomes ${{chroms_array[$index]}} --chunks --bed {params.prefix} {input.index} ${{chunks_array[$index]}}
        done
        
        beds_list='{params.beds}'
        declare -a beds_array=($beds_list)       
 
        for bed in "${{beds_array[@]}}"; do
            freebayes -f {params.ref} -t $bed.bed {input.masked_bam} > "${{bed//beds/vcfs}}.vcf" 2>> {log}
        done
        
        bcftools concat {params.vcfs} | vcfuniq > {output} 2>> {log}
        """
