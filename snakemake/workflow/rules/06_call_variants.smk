# Rule for masking converted bases in BAM files
rule mask_converted_bases:
    input:
        "{output_path}/{sample}/04_deduplicate_bam/{sample}_sd.bam",  # Input deduplicated BAM file
    output:
        # Output BAM files and their index files after calmd and masking
        calmd_bam="{output_path}/{sample}/06_call_variants/{sample}_sd_calmd.bam",
        calmd_bai="{output_path}/{sample}/06_call_variants/{sample}_sd_calmd.bam.bai",
        masked_bam="{output_path}/{sample}/06_call_variants/{sample}_sd_calmd_masked.bam",
        masked_bai="{output_path}/{sample}/06_call_variants/{sample}_sd_calmd_masked.bam.bai",
    params:
        # Location of the relevant genome fasta file
        genome_fa=lambda wcs: D_sample_details[wcs.sample]["genome_path"] +
            '/' + D_sample_details[wcs.sample]["genome"] + '.fa',
        # Parallel threads for running samtools, one less than the total available threads
        parallel = lambda wcs, threads: threads - 1,
    log:
        "{output_path}/{sample}/logs/{sample}_mask_converted_bases.log",
    conda:
        "../envs/revelio.yaml",
    threads: lambda wildcards: get_cpus(1,64),
    resources:
        # Resource requirements foddr running the rule
        time_min=lambda wildcards, input, threads: get_time_min(wildcards, input, "mask_converted_bases", threads),
        cpus=lambda wildcards, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        # Command to perform calmd operation using samtools, effectively recalibrating the MD and NM tags in BAM files
        "samtools calmd -b {input} {params.genome_fa} -@ {params.parallel} "
        "1> {output.calmd_bam} 2>> /dev/null \n\n"

        # Index the calmd BAM file
        "samtools index -@ {params.parallel} {output.calmd_bam} &>> {log} \n\n"

        # Python script to mask the converted bases in the calmd BAM file
        "python3 {workflow.basedir}/scripts/revelio.py -t {wildcards.output_path}/{wildcards.sample} "
        "-T {params.parallel} -Q {output.calmd_bam} {output.masked_bam} 2>> {log} \n\n"

        # Index the masked BAM file
        "samtools index -@ {params.parallel} {output.masked_bam} &>> {log}"


# Rule calls variants using the FreeBayes tool for variant calling.
rule call_variants:
    input:
        # Define the path to the masked BAM file.
        masked_bam=lambda wcs: os.path.join(
            D_sample_details[wcs.sample]["output_path"],
            wcs.sample,
            "06_call_variants",
            wcs.sample + "_sd_calmd_masked.bam",
        ),
        # Define the path to the genome index file.
        index=lambda wcs: os.path.join(D_sample_details[wcs.sample]["genome_path"], D_sample_details[wcs.sample]["genome"] + ".fa.fai"),
    output:
        # Define the path to the output VCF file.
        "{output_path}/{sample}/06_call_variants/{sample}_sd.vcf",
    params:
        # Define the path to the reference genome FASTA file.
        ref=lambda wcs: os.path.join(D_sample_details[wcs.sample]["genome_path"],
            D_sample_details[wcs.sample]["genome"] + ".fa"),
        # Define the path to the prefix for intermediate files.
        prefix=lambda wcs: os.path.join(wcs.output_path, wcs.sample,"06_call_variants/beds", D_sample_details[wcs.sample]["genome"]),
        # Calculate paths for various BED files and VCF files.
        beds=lambda wcs: calculate_vcf_files(wcs,"beds",""),
        vcfs=lambda wcs: calculate_vcf_files(wcs,"vcfs",".vcf"),
        # Calculate the number of chromosomes and their names.
        num_chroms=lambda wcs: len(D_sample_details[wcs.sample]["nchunks"]) - 1,
        chroms=lambda wcs: [chrom for chrom in D_sample_details[wcs.sample]["nchunks"]],
        chunks=lambda wcs: [D_sample_details[wcs.sample]["nchunks"][chrom] for chrom in D_sample_details[wcs.sample]["nchunks"]],
        # Determine if a region of interest (ROI) BED file exists.
        roi_bed=lambda wcs: D_sample_details[wcs.sample].get("roi_bed",""),
    log:
        "{output_path}/{sample}/logs/{sample}_call_variants.log",
    conda:
        "../envs/freebayes.yaml"
    threads: 53
    resources:
        # Specify resource requirements.
        time_min=lambda wildcards, input, threads: get_time_min(wildcards, input, "call_variants", threads),
        cpus=lambda wildcards, threads: threads,
        mem_mb=get_mem_mb,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition=""
    shell:
        """
        # If a region of interest (ROI) BED file is provided, copy it to the beds file.
        if [ -n '{params.roi_bed}' ]; then
            cp {params.roi_bed} {params.beds}.bed
        else
            # Extract chromosome and chunk information for generating BED files.
            chroms_list='{params.chroms}'
            declare -a chroms_array=($chroms_list)

            chunks_list='{params.chunks}'
            declare -a chunks_array=($chunks_list)

            # Loop through each chromosome and generate BED files using fasta_generate_regions.py.
            for index in $(seq 0 {params.num_chroms}); do
                fasta_generate_regions.py --chromosomes ${{chroms_array[$index]}} --chunks --bed {params.prefix} {input.index} ${{chunks_array[$index]}}
            done
        fi

        beds_list='{params.beds}'
        declare -a beds_array=($beds_list)

        # Loop through each BED file and call variants using FreeBayes.
        for bed in "${{beds_array[@]}}"; do
            # Variant calling with FreeBayes, outputting VCF files.
            freebayes -f {params.ref} -t $bed.bed {input.masked_bam} > "${{bed//beds/vcfs}}.vcf" 2>> {log}
        done

        # Concatenate and merge VCF files using bcftools and vcfuniq to produce the final VCF output.
        bcftools concat {params.vcfs} | vcfuniq > {output} 2>> {log}
        """

