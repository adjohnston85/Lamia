rule cleanup:
    input:
        lambda wcs: D_sample_details[wcs.sample]["cleanup_inputs"]
    output:
        bedGraph="{output_path}/{sample}/{sample}_sd_CpG.bedGraph.gz",
        bam="{output_path}/{sample}/{sample}_sd.bam",
    resources:
        time_min=lambda wildcards, input, threads: get_time_min(wildcards, input, "cleanup", threads),
        mem_mb=get_mem_mb,
        cpus=1,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="--partition=io"
    run:
        print(D_sample_details[wildcards.sample]["cleanup"])
        if not os.path.exists(D_sample_details[wildcards.sample]["cleanup"]):
            shell(
                'if [ -e "{wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam/{wildcards.sample}_sd.bam" ]; then ' 
                    'mv {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam/{wildcards.sample}_sd.{{bam,bam.bai}} '
                        '{wildcards.output_path}/{wildcards.sample}/ ; '
                'fi \n'
                
                'rm {wildcards.output_path}/{wildcards.sample}/02_trim_fastq/*.fq.gz \n'
                'rm {wildcards.output_path}/{wildcards.sample}/{{03_align_fastq,04_deduplicate_bam}}/*.bam* \n'
                'rm -r {wildcards.output_path}/{wildcards.sample}/06_call_variants/{{beds,vcfs}} \n'
                
                'mkdir -p {wildcards.output_path}/{wildcards.sample}/QC/{{bam_deduplication,fastq_trimming,fastq_alignment}} \n'
                'mkdir -p {wildcards.output_path}/{wildcards.sample}/stats/methyldackel \n'
                'mkdir {wildcards.output_path}/{wildcards.sample}/MethylSeekR \n'
                'mkdir {wildcards.output_path}/{wildcards.sample}/browser_tracks \n'
    
                'mv {wildcards.output_path}/{wildcards.sample}/02_trim_fastq/ '
                    '{wildcards.output_path}/{wildcards.sample}/QC/fastq_trimming/ \n'
                
                'mv {wildcards.output_path}/{wildcards.sample}/03_align_fastq/ '
                    '{wildcards.output_path}/{wildcards.sample}/QC/fastq_alignment/ \n'
                
                'mv {wildcards.output_path}/{wildcards.sample}/04_deduplicate_bam/ '
                    '{wildcards.output_path}/{wildcards.sample}/QC/bam_deduplication/ \n'
                
                'mv {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_sd_CpG.bedGraph '
                    '{wildcards.output_path}/{wildcards.sample}/ \n'
                
                'mv {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}*.svg '
                    '{wildcards.output_path}/{wildcards.sample}/stats/methyldackel/ \n'
                
                'mv {wildcards.output_path}/{wildcards.sample}/05_call_methylation/{wildcards.sample}_ROI_Conversion*.bedGraph '
                    '{wildcards.output_path}/{wildcards.sample}/stats/methyldackel/ \n'
                
                'mv {wildcards.output_path}/{wildcards.sample}/06_call_variants/{wildcards.sample}_sd.vcf '
                    '{wildcards.output_path}/{wildcards.sample}/ \n'
                
                'mv {wildcards.output_path}/{wildcards.sample}/07_calculate_statistics/{wildcards.sample}*.txt '
                    '{wildcards.output_path}/{wildcards.sample}/stats/ \n'
                
                'mv {wildcards.output_path}/{wildcards.sample}/08_methylseekr_and_TDF/{wildcards.sample}_sd_CpG.tdf '
                    '{wildcards.output_path}/{wildcards.sample}/browser_tracks/ \n'
                
                'if grep -q "MethylSeekR not run, coverage < 10x" {wildcards.output_path}/{wildcards.sample}/08_methylseekr_and_TDF/{wildcards.sample}_PMD.bed; then '
                    'echo "MethylSeekR not run, inadequate coverage"; '
                    'touch {wildcards.output_path}/{wildcards.sample}/MethylSeekR/InadequateCoverage; '
                'else '
                    'mv {wildcards.output_path}/{wildcards.sample}/08_methylseekr_and_TDF/{{*UMR*,*PMD*,*CalculateFDR*,*Segmentation*,*AlphaDistribution*}} '
                        '{wildcards.output_path}/{wildcards.sample}/MethylSeekR/; '
                'fi \n'
                
                'find {wildcards.output_path}/{wildcards.sample} -type f \( -name "*.txt" -o -name "*.bed" -o -name "*.bedGraph" '
                    '-o -name "*.vcf" -o -name "*.tdf" -o -name "*.log" -o -name "*.std*" \) -exec gzip {{}} \;'
                )  
            
        for file_path in input:
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            if not os.path.exists(file_path):
                with open(file_path, 'w') as file:
                    pass

rule rsync:
    input:
        "{output_path}/{sample}/{sample}_sd_CpG.bedGraph.gz"
    output:
        "{output_path}/{sample}/rsync_complete.txt",
    params:
        rsync_path = lambda wcs: D_sample_details[wcs.sample]['rsync']
    resources:
        time_min=lambda wildcards, input, threads: get_time_min(wildcards, input, "rsync", threads),
        mem_mb=get_mem_mb,
        cpus=1,
        account=lambda wcs: D_sample_details[wcs.sample]['account'],
        email=lambda wcs: D_sample_details[wcs.sample]['email'],
        partition="--partition=io"
    shell:
        "mkdir -p {params.rsync_path}/{wildcards.sample}/logs\n"

        "rsync -avzh --whole-file --no-g --chmod=Dg=rwxs "
        "{wildcards.output_path}/{wildcards.sample} {params.rsync_path}/ "
        "&>> {params.rsync_path}/{wildcards.sample}/logs/{wildcards.sample}_rsync.log\n"

        "gzip {params.rsync_path}/{wildcards.sample}/logs/{wildcards.sample}_rsync.log\n"

        "touch {output}"
