def perform_cleanup(sample, output_path, umi_len=None):
    umi_len_check = "true" if umi_len else "false"
    shell(
        """
        # Move BAM files based on umi_len condition
        if [ "{umi_len_check}" = "false" ]; then
            [ -d "{output_path}/{sample}/03_align_fastq" ] && \
                mv {output_path}/{sample}/04_deduplicate_bam/{sample}_dedup_merged_and_sorted_se_pe.bam \
                {output_path}/{sample}/ ;
        else
            [ -d "{output_path}/{sample}/04_deduplicate_bam" ] && \
                mv {output_path}/{sample}/04_deduplicate_bam/{sample}_consensus_merged_and_sorted_se_pe.bam \
                {output_path}/{sample}/ ;
        fi

        # Check if directory exists before finding and removing files
        dir="{output_path}/{sample}/02_trim_fastq"; 
        [ -d "$dir" ] && find "$dir" -maxdepth 1 -name "*.fq.gz" -type f -exec rm {{}} +

        # Loop through subdirectories and check if they exist before removing files
        for subdir in 03_align_fastq 04_deduplicate_bam 07_calculate_statistics; do 
            dir="{output_path}/{sample}/$subdir"; 
            [ -d "$dir" ] && find "$dir" -maxdepth 1 \( -name "*.bam*" -o -name "*.bai*" -o -name "*.fastq*" \) \
                -exec rm {{}} +; 
        done

        # Check if directories exist before removing them
        for dir in {output_path}/{sample}/06_call_variants/{{beds,vcfs}}; do 
            [ -d "$dir" ] && rm -r "$dir"; 
        done

        # Create directories
        mkdir -p {output_path}/{sample}/methyldackel
        mkdir -p {output_path}/{sample}/stats
        mkdir -p {output_path}/{sample}/QC/{{bam_deduplication,fastq_trimming,fastq_alignment}} 
        mkdir -p {output_path}/{sample}/MethylSeekR 
        mkdir -p {output_path}/{sample}/browser_tracks 

        # Move files if source directories exist
        [ -d "{output_path}/{sample}/02_trim_fastq" ] && \
            rsync -av {output_path}/{sample}/02_trim_fastq/ \
                {output_path}/{sample}/QC/fastq_trimming 

        [ -d "{output_path}/{sample}/03_align_fastq" ] && \
            rsync -av {output_path}/{sample}/03_align_fastq/ \
                {output_path}/{sample}/QC/fastq_alignment 

        [ -d "{output_path}/{sample}/04_deduplicate_bam" ] && \
            rsync -av {output_path}/{sample}/04_deduplicate_bam/ \
                {output_path}/{sample}/QC/bam_deduplication 

        [ -d "{output_path}/{sample}/05_call_methylation" ] && \
            rsync -av {output_path}/{sample}/05_call_methylation/ \
                {output_path}/{sample}/methyldackel

        if [ -f {output_path}/{sample}/06_call_variants/{sample}_sd.vcf ]; then 
            mv {output_path}/{sample}/06_call_variants/{sample}_sd.vcf \
                {output_path}/{sample}/ ; 
        fi 

        [ -d "{output_path}/{sample}/07_calculate_statistics" ] && \
            rsync -av {output_path}/{sample}/07_calculate_statistics/ \
                {output_path}/{sample}/stats

        if [ -f {output_path}/{sample}/08_methylseekr_and_TDF/{sample}_sd_CpG.tdf ]; then 
            mv {output_path}/{sample}/08_methylseekr_and_TDF/{sample}_sd_CpG.tdf \
                {output_path}/{sample}/browser_tracks/ ; 
        fi 

        [ -d "{output_path}/{sample}/08_methylseekr_and_TDF" ] && \
            rsync -av {output_path}/{sample}/08_methylseekr_and_TDF/ \
                {output_path}/{sample}/MethylSeekR

        # Gzip specific file types
        find {output_path}/{sample} -type f \\( -name "*.txt" -o -name "*.bed" -o -name "*.bedGraph" \
            -o -name "*.vcf" -o -name "*.tdf" -o -name "*.log" -o -name "*.std*" \\) -print -exec gzip -f {{}} \\;

        rm -r {output_path}/{sample}/01_sequence_files
        rm -r {output_path}/{sample}/02_trim_fastq
        rm -r {output_path}/{sample}/03_align_fastq
        rm -r {output_path}/{sample}/04_deduplicate_bam
        rm -r {output_path}/{sample}/05_call_methylation
        rm -r {output_path}/{sample}/06_call_variants
        rm -r {output_path}/{sample}/07_calculate_statistics
        rm -r {output_path}/{sample}/08_methylseekr_and_TDF
        """
    )

def perform_rsync(sample, output_path, rsync_path):
    shell(
        """
        # Create necessary directories for rsync logs
        mkdir -p {rsync_path}/{sample}/logs

        # Perform rsync synchronization with options
        rsync -avzh --whole-file --no-g --chmod=Dg=rwxs \
        {output_path}/{sample} {rsync_path}/ \
        &>> {rsync_path}/{sample}/logs/{sample}_rsync.log

        # Gzip the rsync log file
        gzip -f {rsync_path}/{sample}/logs/{sample}_rsync.log

        # Indicate rsync completion
        touch {rsync_path}/{sample}/rsync_complete.txt
        """
    )

