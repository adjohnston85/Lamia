rule call_methylation:
    input:
        s_bam="{output_path}/{sample}/03_align_fastq/{sample}_merged_and_sorted.bam",
        sd_bam=lambda wildcards: (
            "{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe_duplex.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
            if "umi_len" in D_sample_details[wildcards.sample]
            else "{output_path}/{sample}/04_deduplicate_bam/{sample}_dedup_merged_and_sorted_se_pe.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
        ),
    output:
        # Outputs for sorted BAM
        "{output_path}/{sample}/05_call_methylation/{sample}_s_ROI_Conversion_CHH.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_s_ROI_Conversion_CHG.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_s_ROI_Conversion_CpG.bedGraph",
        # Outputs for deduplicated BAM
        "{output_path}/{sample}/05_call_methylation/{sample}_sd_ROI_Conversion_CHH.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_sd_ROI_Conversion_CHG.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_sd_ROI_Conversion_CpG.bedGraph",
        "{output_path}/{sample}/05_call_methylation/{sample}_sd_CpG.bedGraph",
        # M-bias plots
        CHH="{output_path}/{sample}/05_call_methylation/{sample}_CHH_mbias.txt",
        CHG="{output_path}/{sample}/05_call_methylation/{sample}_CHG_mbias.txt",
        CpG="{output_path}/{sample}/05_call_methylation/{sample}_CpG_mbias.txt",
    params:
        ROI_path=lambda wcs: D_sample_details[wcs.sample]["roi_bed"]
        if "roi_bed" in D_sample_details[wcs.sample]
        else D_sample_details[wcs.sample]["genome_path"] + "/CHH_ROI.bed",
        genome_fa=lambda wcs: D_sample_details[wcs.sample]["genome_path"] + '/' + D_sample_details[wcs.sample]["genome"] + ".fa",
        se_trim_options=lambda wcs: "--nOT {0},{0},0,0 --nOB {0},{0},0,0".format(
            D_sample_details[wcs.sample]["trim_profile"][0]),
        pe_trim_options=lambda wcs: "--nOT {0},0,{1},0 --nOB {0},0,{1},0".format(
            D_sample_details[wcs.sample]["trim_profile"][0],
            D_sample_details[wcs.sample]["trim_profile"][2]),
        sdd_bam=lambda wildcards: (
            "{output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe_duplex.bam".format(
                output_path=wildcards.output_path, sample=wildcards.sample
            )
            if "umi_len" in D_sample_details[wildcards.sample]
            else "None"
        ),
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

        # Process SE reads for sorted BAM with ROI
        MethylDackel extract {params.se_trim_options} --CHH --CHG --minOppositeDepth 5 --maxVariantFrac 0.2 \
        --opref {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_se_ROI_Conversion \
        -@ {threads} -l {params.ROI_path} {params.genome_fa} {input.s_bam} \
        --ignoreFlags 0x1 --ignoreFlags 2048 \
        &>> {log}

        # Process PE reads for sorted BAM with ROI
        MethylDackel extract {params.pe_trim_options} --CHH --CHG --minOppositeDepth 5 --maxVariantFrac 0.2 \
        --opref {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_pe_ROI_Conversion \
        -@ {threads} -l {params.ROI_path} {params.genome_fa} {input.s_bam} \
        --requireFlags 0x1 --ignoreFlags 2048 \
        &>> {log}

        # Merge SE and PE bedGraph files for sorted BAM with ROI
        cat {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_se_ROI_Conversion_CpG.bedGraph \
        {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_pe_ROI_Conversion_CpG.bedGraph \
        > {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_ROI_Conversion_CpG.bedGraph

        cat {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_se_ROI_Conversion_CHG.bedGraph \
        {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_pe_ROI_Conversion_CHG.bedGraph \
        > {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_ROI_Conversion_CHG.bedGraph

        cat {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_se_ROI_Conversion_CHH.bedGraph \
        {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_pe_ROI_Conversion_CHH.bedGraph \
        > {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_ROI_Conversion_CHH.bedGraph

        # Process SE reads for duduplicated BAM with ROI
        MethylDackel extract {params.se_trim_options} --CHH --CHG --minOppositeDepth 5 --maxVariantFrac 0.2 \
        --opref {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_se_ROI_Conversion \
        -@ {threads} -l {params.ROI_path} {params.genome_fa} {input.sd_bam} \
        --ignoreFlags 0x1 --ignoreFlags 2048 \
        &>> {log}

        # Process PE reads for duduplicated BAM with ROI
        MethylDackel extract {params.pe_trim_options} --CHH --CHG --minOppositeDepth 5 --maxVariantFrac 0.2 \
        --opref {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_pe_ROI_Conversion \
        -@ {threads} -l {params.ROI_path} {params.genome_fa} {input.sd_bam} \
        --requireFlags 0x1 --ignoreFlags 2048 \
        &>> {log}

        # Merge SE and PE bedGraph files for sorted BAM with ROI
        cat {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_se_ROI_Conversion_CpG.bedGraph \
        {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_pe_ROI_Conversion_CpG.bedGraph \
        > {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_ROI_Conversion_CpG.bedGraph

        cat {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_se_ROI_Conversion_CHG.bedGraph \
        {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_pe_ROI_Conversion_CHG.bedGraph \
        > {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_ROI_Conversion_CHG.bedGraph

        cat {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_se_ROI_Conversion_CHH.bedGraph \
        {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_pe_ROI_Conversion_CHH.bedGraph \
        > {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_ROI_Conversion_CHH.bedGraph

        # Process SE reads for deduplicated BAM with mergeContext
        MethylDackel extract {params.se_trim_options} --mergeContext \
        -@ {threads} {params.genome_fa} {input.sd_bam} \
        -o {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_se \
        --ignoreFlags 0x1 --ignoreFlags 2048 \
        &>> {log}

        # Process PE reads for deduplicated BAM with mergeContext
        MethylDackel extract {params.pe_trim_options} --mergeContext \
        -@ {threads} {params.genome_fa} {input.sd_bam} \
        -o {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_pe \
        --requireFlags 0x1 --ignoreFlags 2048 \
        &>> {log}

        # Merge SE and PE bedGraph files for deduplicated BAM
        cat {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_se_CpG.bedGraph \
        {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_pe_CpG.bedGraph \
        > {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_CpG.bedGraph

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

        # Conditionally process duplex BAM
        if [ "{params.sdd_bam}" != "None" ]; then
            # Process SE reads for duplex BAM with ROI
            MethylDackel extract {params.se_trim_options} --CHH --CHG --minOppositeDepth 5 --maxVariantFrac 0.2 \
            --opref {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_duplex_se_ROI_Conversion \
            -@ {threads} -l {params.ROI_path} {params.genome_fa} {params.sdd_bam} \
            --ignoreFlags 0x1 --ignoreFlags 2048 \
            &>> {log}

            # Process PE reads for duplex BAM with ROI
            MethylDackel extract {params.pe_trim_options} --CHH --CHG --minOppositeDepth 5 --maxVariantFrac 0.2 \
            --opref {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_duplex_pe_ROI_Conversion \
            -@ {threads} -l {params.ROI_path} {params.genome_fa} {params.sdd_bam} \
            --requireFlags 0x1 --ignoreFlags 2048 \
            &>> {log}

            # Merge SE and PE bedGraph files for duplex BAM with ROI
            cat {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_duplex_se_ROI_Conversion_CpG.bedGraph \
            {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_duplex_pe_ROI_Conversion_CpG.bedGraph \
            > {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_duplex_ROI_Conversion_CpG.bedGraph

            rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_duplex_se_ROI_Conversion_CpG.bedGraph \
            {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_duplex_pe_ROI_Conversion_CpG.bedGraph
        fi

        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_se_ROI_Conversion_CpG.bedGraph
        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_pe_ROI_Conversion_CpG.bedGraph
        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_se_ROI_Conversion_CHG.bedGraph
        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_pe_ROI_Conversion_CHG.bedGraph
        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_se_ROI_Conversion_CHH.bedGraph
        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_s_pe_ROI_Conversion_CHH.bedGraph
        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_se_ROI_Conversion_CpG.bedGraph
        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_pe_ROI_Conversion_CpG.bedGraph
        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_se_ROI_Conversion_CHG.bedGraph
        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_pe_ROI_Conversion_CHG.bedGraph
        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_se_ROI_Conversion_CHH.bedGraph
        # rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_pe_ROI_Conversion_CHH.bedGraph
        rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_se_CpG.bedGraph
        rm {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_pe_CpG.bedGraph     
        """