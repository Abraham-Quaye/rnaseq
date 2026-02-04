############### DOWNLOAD GRCH38 GENOME AND GTF FROM ENSEMBL ################
rule download_reference_files:
    input:
        "scripts/shell/get_reference_files.sh"
    output:
        genome = "raw_files/genome_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        zipped_genome = "raw_files/genome_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        gtf = "raw_files/annotations/Homo_sapiens.GRCh38.115.gtf.gz",
        transcriptome = "raw_files/genome_files/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    shell:
        """
        {input}

        echo "Unzipping Genome..."
        gunzip -k {output.genome}.gz
        """

################ GENERATE REQUISITE FILES FOR BUILDING SALMON GENOMIC INDEX ###############
rule make_salmon_genomic_index_files:
    input:
        genome = rules.download_reference_files.output.genome,
        z_genome = rules.download_reference_files.output.zipped_genome,
        trxptome = rules.download_reference_files.output.transcriptome
    output:
        decoy = "raw_files/genome_files/GRCh38_decoys.txt",
        gentrome = "raw_files/genome_files/GRCh38_gentrome.fa.gz"
    shell:
        """
        echo "Making Decoys file ...";
        grep "^>" {input.genome} | cut -d " " -f 1 | sed 's/>//' > {output.decoy}

        echo "Making Gentrome ..."
        # order is important: trxptome must be first in the concatenated file
        cat {input.trxptome} {input.z_genome} > {output.gentrome}
        """
        
################ BUILD GRCh38 GENOMIC INDEX FOR MAPPING WITH SALMON ###############
rule build_salmon_genomic_index:
    input:
        gentrome = rules.make_salmon_genomic_index_files.output.gentrome,
        decoy = rules.make_salmon_genomic_index_files.output.decoy
    output:
        directory("raw_files/salmon_GRCh38_genome_index")
    shell:
        """
        mkdir -p {output}
        
        echo "Bulding GRCh38 Genomic Index for Salmon ..."
        salmon index -t {input.gentrome} -d {input.decoy} -p 15 -i {output} -k 31
        """

################ BUILD GRCh38 GENOMIC INDEX FOR MAPPING WITH STAR ###############
rule build_star_genomic_index:
    input: 
        genome = rules.download_reference_files.output.genome,
        gtf = rules.download_reference_files.output.gtf
    output:
        star_genome_dir = directory("raw_files/star_GRCh38_genome_index")
    params:
        unzip_gtf = "raw_files/annotations/Homo_sapiens.GRCh38.115.gtf"
    shell:
        """
        mkdir -p {output.star_genome_dir}

        echo "Unzipping GTF ..."
        gunzip -k {input.gtf}

        echo "Building STAR Genomic Index Files ..."

        STAR --runThreadN 15 --runMode genomeGenerate --genomeDir {output.star_genome_dir} \
        --genomeFastaFiles {input.genome} --sjdbGTFfile {params.unzip_gtf} --sjdbOverhang 149
        """

########## CONCATENATE SAME SAME STRAND READS FROM DIFFERENT LANES #########
myocd_dir = "/Users/abrahamquaye/myocd_rnaseq"
rule concatenate_same_strand_reads:
    input:
        ctrl_1 = expand(f"{myocd_dir}/raw_files/siO_GFP_1_S13/siO_GFP_1_S13_L00{{lane}}_R{{strand}}_001.fastq.gz", \
        lane = [1, 2], strand = [1, 2]),
        ctrl_2 = expand(f"{myocd_dir}/raw_files/siO_GFP_3_S15/siO_GFP_3_S15_L00{{lane}}_R{{strand}}_001.fastq.gz", \
        lane = [1, 2], strand = [1, 2]),
        ctrl_3 = expand(f"{myocd_dir}/raw_files/siO_GFP_2_S14/siO_GFP_2_S14_L00{{lane}}_R{{strand}}_001.fastq.gz", \
        lane = [1, 2], strand = [1, 2]),
        treat_1 = expand(f"{myocd_dir}/raw_files/siO_MYOCD_1_S16/siO_MYOCD_1_S16_L00{{lane}}_R{{strand}}_001.fastq.gz", \
        lane = [1, 2], strand = [1, 2]),
        treat_2 = expand(f"{myocd_dir}/raw_files/siO_MYOCD_2_S17/siO_MYOCD_2_S17_L00{{lane}}_R{{strand}}_001.fastq.gz", \
        lane = [1, 2], strand = [1, 2]),
        treat_3 = expand(f"{myocd_dir}/raw_files/siO_MYOCD_3_S18/siO_MYOCD_3_S18_L00{{lane}}_R{{strand}}_001.fastq.gz", \
        lane = [1, 2], strand = [1, 2]),
        script = "scripts/shell/concatenate_reads.sh"
    output: 
        ctrl_1 = expand(f"{myocd_dir}/raw_files/merged_fastqs/siO_GFP_1_S13_R{{strand}}_merged.fastq.gz", \
        strand = [1, 2]),
        ctrl_2 = expand(f"{myocd_dir}/raw_files/merged_fastqs/siO_GFP_2_S14_R{{strand}}_merged.fastq.gz", \
        strand = [1, 2]),
        ctrl_3 = expand(f"{myocd_dir}/raw_files/merged_fastqs/siO_GFP_3_S15_R{{strand}}_merged.fastq.gz", \
        strand = [1, 2]),
        treat_1 = expand(f"{myocd_dir}/raw_files/merged_fastqs/siO_MYOCD_1_S16_R{{strand}}_merged.fastq.gz", \
        strand = [1, 2]),
        treat_2 = expand(f"{myocd_dir}/raw_files/merged_fastqs/siO_MYOCD_2_S17_R{{strand}}_merged.fastq.gz", \
        strand = [1, 2]),
        treat_3 = expand(f"{myocd_dir}/raw_files/merged_fastqs/siO_MYOCD_3_S18_R{{strand}}_merged.fastq.gz", \
        strand = [1, 2])
    shell:
        """
        {input.script}
        """

############### TRIM READS WITH TRIM GALORE ####################
rule trim_fastq_files:
    input:
        raw_fastqs = rules.concatenate_same_strand_reads.output,
        script = "scripts/shell/trim_reads.sh"
    output:
        t_fastqs = expand(f"{myocd_dir}/results/trimmedReads/{{sample}}_R{{strand}}_merged_val_{{strand}}.fq.gz", \
        sample = ["siO_GFP_1_S13", "siO_GFP_2_S14", "siO_GFP_3_S15", "siO_MYOCD_1_S16", \
        "siO_MYOCD_2_S17", "siO_MYOCD_3_S18"], strand = [1, 2])
    shell:
        "{input.script}"

############### QC TRIMMED READS WITH FASTQC ####################
rule FastQC_trimmed_reads:
    input:
        t_fastqs = rules.trim_fastq_files.output.t_fastqs
    output:
        t_fastqc = expand(f"{myocd_dir}/results/fastqc/{{sample}}_R{{strand}}_merged_val_{{strand}}_fastqc.html", \
        sample = ["siO_GFP_1_S13", "siO_GFP_2_S14", "siO_GFP_3_S15", "siO_MYOCD_1_S16", \
        "siO_MYOCD_2_S17", "siO_MYOCD_3_S18"], strand = [1, 2]),
        fastqc_dir = directory(f"{myocd_dir}/results/fastqc")
    shell:
        """
        mkdir -p {output.fastqc_dir}
        fastqc -t 16 --memory 1024 -o {output.fastqc_dir} {input.t_fastqs}
        """

#################### QC READS WITH MULTIQC ####################
mqc_dir = f"{myocd_dir}/results/multiqc"

rule MultiQC_all_fastqcs:
    input:
        rules.FastQC_trimmed_reads.output.fastqc_dir,
        rules.FastQC_trimmed_reads.output.t_fastqc
    output:
        f"{mqc_dir}/multiqc_report.html",
        directory(f"{mqc_dir}/multiqc_data")
    shell:
        """
        mkdir -p {mqc_dir}
        multiqc -f -o {mqc_dir} {input}
        """

####### MAP READS WITH STAR #############
rule STAR_map_quant_reads:
    input:
        t_fastqs = rules.trim_fastq_files.output.t_fastqs,
        star_index = rules.build_star_genomic_index.output.star_genome_dir,
        gtf = rules.build_star_genomic_index.params.unzip_gtf
    output:
        bams = expand(f"{myocd_dir}/results/star_mapped/{{sample}}/Aligned.sortedByCoord.out.bam", \
        sample = ["siO_GFP_1_S13", "siO_GFP_2_S14", "siO_GFP_3_S15", "siO_MYOCD_1_S16", \
        "siO_MYOCD_2_S17", "siO_MYOCD_3_S18"]),
        genecounts = expand(f"{myocd_dir}/results/star_mapped/{{sample}}/ReadsPerGene.out.tab", \
        sample = ["siO_GFP_1_S13", "siO_GFP_2_S14", "siO_GFP_3_S15", "siO_MYOCD_1_S16", \
        "siO_MYOCD_2_S17", "siO_MYOCD_3_S18"])
    params:
        trim_dir = f"{myocd_dir}/results/trimmedReads",
        star_dir = f"{myocd_dir}/results/star_mapped"
    shell:
        """
        mkdir -p {params.star_dir}

        forward_fastqs=( $(ls {input.t_fastqs} | grep '_val_1.fq.gz') )

        for i in ${{forward_fastqs[@]}}; do
            sample_id=$(basename ${{i}} | cut -d "_" -f 1-4);
            echo "Mapping reads for ${{sample_id}} ...";

            ulimit -n 100000

            STAR --runThreadN 14 --genomeDir {input.star_index} \\
            --readFilesIn ${{i}} {params.trim_dir}/${{sample_id}}_R2_merged_val_2.fq.gz \\
            --readFilesCommand gzcat \\
            --outFileNamePrefix {params.star_dir}/${{sample_id}}/ \\
            --outSAMtype BAM SortedByCoordinate \\
            --quantMode TranscriptomeSAM GeneCounts \\
            --sjdbGTFfile {input.gtf} \\
            --outSAMattributes NH HI AS nM XS \\
            --limitBAMsortRAM 60000000000;
        done
        """

#################### INDEX STAR ALIGNED BAM FILES ####################
rule index_star_bam_files:
    input:
        bams = rules.STAR_map_quant_reads.output.bams
    output:
        indexed_bams = expand(f"{myocd_dir}/results/star_mapped/{{sample}}/Aligned.sortedByCoord.out.bam.bai", \
        sample = ["siO_GFP_1_S13", "siO_GFP_2_S14", "siO_GFP_3_S15", "siO_MYOCD_1_S16", \
        "siO_MYOCD_2_S17", "siO_MYOCD_3_S18"])
    shell:
        """
        for bam in {input.bams}; do
            echo "Indexing ${{bam}} ..."
            samtools index -M -@ 10 ${{bam}}
        done
        """

#################### QUANTIFY READS WITH SALMON ####################
rule quantify_reads_salmon:
    input:
        t_fastqs = rules.trim_fastq_files.output.t_fastqs,
        salmon_index = rules.build_salmon_genomic_index.output
    output:
        expand(f"{myocd_dir}/results/salmon_quant/quant_siO_{{id}}/quant.sf", \
        id = ["GFP_1_S13", "GFP_2_S14", "GFP_3_S15", "MYOCD_1_S16", "MYOCD_2_S17", \
        "MYOCD_3_S18"])
    params:
        trim_dir = f"{myocd_dir}/results/trimmedReads",
        salmon_dir = f"{myocd_dir}/results/salmon_quant"
    shell:
        """
        mkdir -p {params.salmon_dir}

        forward_fastqs=( $(ls {input.t_fastqs} | grep '_val_1.fq.gz') )

        for i in ${{forward_fastqs[@]}}; do
            sample_id=$(basename ${{i}} | cut -d "_" -f 1-4);
            echo "Quantifying reads for ${{sample_id}} ...";

            salmon quant -i {input.salmon_index} -l A \\
            -1 ${{i}} \\
            -2 {params.trim_dir}/${{sample_id}}_R2_merged_val_2.fq.gz \\
            -p 15 --validateMappings -o {params.salmon_dir}/quant_${{sample_id}};
        done
        """
         
##################  DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS OF SALMON QUANT FILES #############
rule DESeq2_salmon_DE_analysis:
    input:
        salmon_quant = rules.quantify_reads_salmon.output,
        r_script = "scripts/r_code/deseq2_salmon_analysis.R",
        r_script2 = "scripts/r_code/DEG_plotting_functions.R",
    output:
        tables = expand(f"{myocd_dir}/results/r/tables/{{res_type}}_MYOCD_vs_GFP_DEGs.csv", \
        res_type = ["significant", "total"]),
        figs = expand(f"{myocd_dir}/results/r/figures/{{fig_type}}_MYOCD_vs_GFP.pdf", \
        fig_type = ["volcano", "heatmap", "pca", "dists"])
    shell:
        """
        {input.r_script}
        rm Rplots.pdf
        """

############# PLOT DEG BAR PLOTS #############
rule plot_deg_barplots:
    input:
        deg_files = rules.DESeq2_salmon_DE_analysis.output.tables,
        r_script = "scripts/r_code/plot_DEG_barplot.R"
    output:
        deg_barplots = f"{myocd_dir}/results/r/figures/DEG_levels_barplot.pdf"
    shell:
        "{input.r_script}"

############ FUNCTIONAL ENRICHMENT ANALYSES OF DEGs #############
# rule functional_enrichment_analysis:
#    input:
#        deg_files = rules.DESeq2_salmon_DE_analysis.output.sigs,
#        r_script = "scripts/r_code/enrichment_analysis.R",
#        rscript2 = "scripts/r_code/enrichment_analysis_functions.R" 
#    output:
#    shell:
#        """
#        {input.r_script}
#        """

############# RUN COMPLETE WORKFLOW #############
rule run_workflow:
    input:
        rules.MultiQC_all_fastqcs.output,
        rules.index_star_bam_files.output,
        rules.plot_deg_barplots.output,
        # rules.plot_deg_barplots.output,
        # rules.functional_enrichment_analysis.output
