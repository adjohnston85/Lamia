def perform_cleanup(sample, output_path):
    shell(
        """
        # Check if file exists before moving
        if [ -e "{output_path}/{sample}/03_align_fastq/{sample}_s.bam" ]; then 
            mv {output_path}/{sample}/03_align_fastq/{sample}_s.{{bam,bam.bai}} \
                {output_path}/{sample}/ ; 
        fi 

        # Check if directory exists before finding and removing files
        dir="{output_path}/{sample}/02_trim_fastq"; 
        [ -d "$dir" ] && find "$dir" -maxdepth 1 -name "*.fq.gz" -type f -exec rm {{}} +

        # Loop through subdirectories and check if they exist before removing files
        for subdir in 03_align_fastq 04_deduplicate_bam; do 
            dir="{output_path}/{sample}/$subdir"; 
            [ -d "$dir" ] && find "$dir" -maxdepth 1 -name "*.bam*" -type f -exec rm {{}} +; 
        done

        # Check if directories exist before removing them
        for dir in {output_path}/{sample}/06_call_variants/{{beds,vcfs}}; do 
            [ -d "$dir" ] && rm -r "$dir"; 
        done

        # Create directories
        mkdir -p {output_path}/{sample}/stats/methyldackel 
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

        # Move files if they exist
        if [ -f {output_path}/{sample}/05_call_methylation/{sample}_sd_CpG.bedGraph ]; then 
            mv {output_path}/{sample}/05_call_methylation/{sample}_sd_CpG.bedGraph \
                {output_path}/{sample}/ ; 
        fi 

        if ls {output_path}/{sample}/05_call_methylation/{sample}*.svg > /dev/null 2>&1; then 
            mv {output_path}/{sample}/05_call_methylation/{sample}*.svg \
                {output_path}/{sample}/stats/methyldackel/ ; 
        fi 

        if ls {output_path}/{sample}/05_call_methylation/{sample}_ROI_Conversion*.bedGraph > /dev/null 2>&1; then 
            mv {output_path}/{sample}/05_call_methylation/{sample}_ROI_Conversion*.bedGraph \
                {output_path}/{sample}/stats/methyldackel/ ; 
        fi 

        if [ -f {output_path}/{sample}/06_call_variants/{sample}_sd.vcf ]; then 
            mv {output_path}/{sample}/06_call_variants/{sample}_sd.vcf \
                {output_path}/{sample}/ ; 
        fi 

        if ls {output_path}/{sample}/07_calculate_statistics/*.txt > /dev/null 2>&1; then 
            mv {output_path}/{sample}/07_calculate_statistics/{sample}*.txt \
                {output_path}/{sample}/stats/ ; 
        fi 

        if [ -f {output_path}/{sample}/08_methylseekr_and_TDF/{sample}_sd_CpG.tdf ]; then 
            mv {output_path}/{sample}/08_methylseekr_and_TDF/{sample}_sd_CpG.tdf \
                {output_path}/{sample}/browser_tracks/ ; 
        fi 

        # Checks for inadequate coverage in MethylSeekR results; logs and flags if inadequate, otherwise moves relevant output files to MethylSeekR dir
        if grep -q "MethylSeekR not run, coverage < 10x" {output_path}/{sample}/08_methylseekr_and_TDF/{sample}_PMD.bed; then 
            echo "MethylSeekR not run, inadequate coverage"; 
            touch {output_path}/{sample}/MethylSeekR/InadequateCoverage; 
        else 
            mv {output_path}/{sample}/08_methylseekr_and_TDF/{{*UMR*,*PMD*,*CalculateFDR*,*Segmentation*,*AlphaDistribution*}} \
                {output_path}/{sample}/MethylSeekR/; 
        fi 

        # Gzip specific file types
        find {output_path}/{sample} -type f \\( -name "*.txt" -o -name "*.bed" -o -name "*.bedGraph" \
            -o -name "*.vcf" -o -name "*.tdf" -o -name "*.log" -o -name "*.std*" \\) -print -exec gzip -f {{}} \\;
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
        touch {output_path}/{sample}/rsync_complete.txt
        """
    )

