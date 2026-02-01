
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

# #################### REMOVE OR TRIM LOW QUALITY READS ####################
prefix = "BB"
sample_nums = range(1, 34)
berges_dir = "/Users/abrahamquaye/berges_rnaseq"

rule trim_fastq_files:
    input:
        fastqs = expand(f"{berges_dir}/raw_files/raw_fastqs/LCS9697_{{id}}_Clean_Data{{strand}}.fq.gz", \
        id = [prefix + str(num) for num in sample_nums], strand = [1, 2]),
        script = "scripts/shell/trim_reads.sh"
    output:
        t_fastqs = expand(f"{berges_dir}/results/trimmedReads/LCS9697_{{id}}_Clean_Data{{strand}}_val_{{strand}}.fq.gz", \
        id = [prefix + str(num) for num in sample_nums], strand = [1, 2])
    shell:
        "{input.script}"

############### QC TRIMMED READS WITH FASTQC ####################
rule FastQC_trimmed_reads:
    input:
        t_fastqs = rules.trim_fastq_files.output.t_fastqs
    output:
        t_fastqc = expand(f"{berges_dir}/results/fastqc/LCS9697_{{id}}_Clean_Data{{strand}}_val_{{strand}}_fastqc.html", \
        id = [prefix + str(num) for num in sample_nums], strand = [1, 2])
    params:
        fastqc_dir = f"{berges_dir}/results/fastqc"
    shell:
        """
        mkdir -p {params.fastqc_dir}

        fastqc -t 16 --memory 1024 -o {params.fastqc_dir} {input.t_fastqs}
        """

#################### QC READS WITH MULTIQC ####################
mqc_dir = f"{berges_dir}/results/multiqc"
rule MultiQC_all_fastqcs:
    input:
        rules.FastQC_trimmed_reads.params.fastqc_dir,
        rules.FastQC_trimmed_reads.output.t_fastqc
    output:
        f"{mqc_dir}/multiqc_report.html",
        directory(f"{mqc_dir}/multiqc_data")
    shell:
        """
        mkdir -p {mqc_dir}
        multiqc -f -o {mqc_dir} {input}
        """

#################### QUANTIFY READS WITH SALMON ####################
rule quantify_reads_salmon:
    input:
        t_fastqs = rules.trim_fastq_files.output.t_fastqs,
        salmon_index = rules.build_salmon_genomic_index.output
    output:
        expand(f"{berges_dir}/results/salmon_quant/quant_{{id}}/quant.sf", \
        id = [prefix + str(num) for num in sample_nums])
    params:
        trim_dir = f"{berges_dir}/results/trimmedReads",
        salmon_dir = f"{berges_dir}/results/salmon_quant"
    shell:
        """
        mkdir -p {params.salmon_dir}

        forward_fastqs=( $(ls {input.t_fastqs} | grep '_val_1.fq.gz') )

        for i in ${{forward_fastqs[@]}}; do
            sample_id=$(basename ${{i}} | cut -d "_" -f2);
            echo "Quantifying reads for LCS9697_${{sample_id}} ...";

            salmon quant -i {input.salmon_index} -l A \\
            -1 {params.trim_dir}/LCS9697_${{sample_id}}_Clean_Data1_val_1.fq.gz \\
            -2 {params.trim_dir}/LCS9697_${{sample_id}}_Clean_Data2_val_2.fq.gz \\
            -p 15 --validateMappings -o {params.salmon_dir}/quant_${{sample_id}};
        done
        """
         
##################  DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS OF SALMON QUANT FILES #############
rule DESeq2_salmon_DE_analysis:
    input:
        salmon_quant = rules.quantify_reads_salmon.output,
        r_script = "scripts/r_code/deseq2_salmon_r77q_analysis.R",
        r_script2 = "scripts/r_code/DEG_plotting_functions.R",
    output:
        sigs = expand(f"{berges_dir}/results/r/tables/significant_{{treat}}{{tp}}hr_vs_mock72hr_DEGs.csv", \
        treat = ["wt", "r77q"], tp = [4, 8, 12, 24, 72]),
        total = expand(f"{berges_dir}/results/r/tables/total_{{treat}}{{tp}}hr_vs_mock72hr_DEGs.csv", \
        treat = ["wt", "r77q"], tp = [4, 8, 12, 24, 72]),
        r77qvswt_sigs = expand(f"{berges_dir}/results/r/tables/significant_r77q{{tp}}hr_vs_wt{{tp}}hr_DEGs.csv", \
        tp = [4, 8, 12, 24, 72]),
        r77qvswt_total = expand(f"{berges_dir}/results/r/tables/total_r77q{{tp}}hr_vs_wt{{tp}}hr_DEGs.csv", \
        tp = [4, 8, 12, 24, 72]),
        volcano_fig = expand(f"{berges_dir}/results/r/figures/volcano_{{treat}}{{tp}}hr_vs_mock72hr.pdf", \
        treat = ["wt", "r77q"], tp = [4, 8, 12, 24, 72]),
        r77qvswt_volcano_fig = expand(f"{berges_dir}/results/r/figures/volcano_r77q{{tp}}hr_vs_wt{{tp}}hr.pdf", \
        tp = [4, 8, 12, 24, 72]),
        heatmap_fig = expand(f"{berges_dir}/results/r/figures/heatmap_{{treat}}{{tp}}hr_vs_mock72hr.pdf", \
        treat = ["wt", "r77q"], tp = [4, 8, 12, 24, 72]),
        r77qvswt_heatmap_fig = expand(f"{berges_dir}/results/r/figures/heatmap_r77q{{tp}}hr_vs_wt{{tp}}hr.pdf", \
        tp = [8, 24, 72]),
        pca_fig_all = f"{berges_dir}/results/r/figures/PCA_all_vs_mock72hr.pdf",
        r77qvswt_pca_figs = expand(f"{berges_dir}/results/r/figures/PCA_r77q{{tp}}hr_vs_wt{{tp}}hr.pdf", \
        tp = [4, 8, 12, 24, 72]),
        corr_fig_all = f"{berges_dir}/results/r/figures/sampleCorrelation_all_vs_mock72hr.pdf",
        r77qvswt_corr_figs = expand(f"{berges_dir}/results/r/figures/sampleCorrelation_r77q{{tp}}hr_vs_wt{{tp}}hr.pdf", \
        tp = [4, 8, 12, 24, 72]) 
    shell:
        """
        {input.r_script}
        rm Rplots.pdf
        """

############# PLOT DEG BAR PLOTS #############
rule plot_deg_barplots:
    input:
        deg_files = rules.DESeq2_salmon_DE_analysis.output.sigs,
        r_script = "scripts/r_code/plot_DEG_barplot.R"
    output:
        deg_barplots = f"{berges_dir}/results/r/figures/DEG_levels_barplot.pdf"
    shell:
        "{input.r_script}"

############# FUNCTIONAL ENRICHMENT ANALYSES OF DEGs #############
rule functional_enrichment_analysis:
    input:
        deg_files = rules.DESeq2_salmon_DE_analysis.output.sigs,
        r_script = "scripts/r_code/enrichment_analysis.R",
        rscript2 = "scripts/r_code/enrichment_analysis_functions.R" 
    output:
        go_res_vs_mock = expand(f"{berges_dir}/results/r/tables/go_{{contr_name}}_totalDEG_sig.csv", \
        contr_name = ["wt4hr_vs_mock72hr", "wt8hr_vs_mock72hr", "wt12hr_vs_mock72hr", \
        "wt24hr_vs_mock72hr", "wt72hr_vs_mock72hr", \
        "r77q4hr_vs_mock72hr", "r77q8hr_vs_mock72hr", "r77q12hr_vs_mock72hr", \
        "r77q24hr_vs_mock72hr", "r77q72hr_vs_mock72hr"]),
        go_res_r77qvswt = expand(f"{berges_dir}/results/r/tables/go_r77q{{tp}}hr_vs_wt{{tp}}hr_totalDEG_sig.csv", \
        tp = [4, 8, 12, 24, 72]),
        go_res_up_vs_mock = expand(f"{berges_dir}/results/r/tables/go_{{contr_name}}_upDEG_sig.csv", \
        contr_name = ["wt4hr_vs_mock72hr", "wt8hr_vs_mock72hr", "wt12hr_vs_mock72hr", \
        "wt24hr_vs_mock72hr", "wt72hr_vs_mock72hr", \
        "r77q4hr_vs_mock72hr", "r77q8hr_vs_mock72hr", "r77q12hr_vs_mock72hr", \
        "r77q24hr_vs_mock72hr", "r77q72hr_vs_mock72hr"]),
        go_res_up_r77qvswt = expand(f"{berges_dir}/results/r/tables/go_r77q{{tp}}hr_vs_wt{{tp}}hr_upDEG_sig.csv", \
        tp = [4, 8, 12, 24, 72]),
        go_res_down_vs_mock = expand(f"{berges_dir}/results/r/tables/go_{{contr_name}}_downDEG_sig.csv", \
        contr_name = ["wt4hr_vs_mock72hr", "wt8hr_vs_mock72hr", "wt12hr_vs_mock72hr", \
        "wt24hr_vs_mock72hr", "wt72hr_vs_mock72hr", \
        "r77q4hr_vs_mock72hr", "r77q8hr_vs_mock72hr", "r77q12hr_vs_mock72hr", \
        "r77q24hr_vs_mock72hr", "r77q72hr_vs_mock72hr"]),
        go_res_down_r77qvswt = expand(f"{berges_dir}/results/r/tables/go_r77q{{tp}}hr_vs_wt{{tp}}hr_downDEG_sig.csv", \
        tp = [4, 8, 12, 24, 72]),
        kegg_res_vs_mock = expand(f"{berges_dir}/results/r/tables/kegg_{{contr_name}}_totalDEG_sigPathways.csv", \
        contr_name = ["wt4hr_vs_mock72hr", "wt8hr_vs_mock72hr", "wt12hr_vs_mock72hr", \
        "wt24hr_vs_mock72hr", "wt72hr_vs_mock72hr", \
        "r77q4hr_vs_mock72hr", "r77q8hr_vs_mock72hr", "r77q12hr_vs_mock72hr", \
        "r77q24hr_vs_mock72hr", "r77q72hr_vs_mock72hr"]),
        kegg_res_r77qvswt = expand(f"{berges_dir}/results/r/tables/kegg_r77q{{tp}}hr_vs_wt{{tp}}hr_totalDEG_sigPathways.csv", \
        tp = [4, 8, 12, 24, 72]),
        kegg_res_up_vs_mock = expand(f"{berges_dir}/results/r/tables/kegg_{{contr_name}}_upDEG_sigPathways.csv", \
        contr_name = ["wt4hr_vs_mock72hr", "wt8hr_vs_mock72hr", "wt12hr_vs_mock72hr", \
        "wt24hr_vs_mock72hr", "wt72hr_vs_mock72hr", \
        "r77q4hr_vs_mock72hr", "r77q8hr_vs_mock72hr", "r77q12hr_vs_mock72hr", \
        "r77q24hr_vs_mock72hr", "r77q72hr_vs_mock72hr"]),
        kegg_res_up_r77qvswt = expand(f"{berges_dir}/results/r/tables/kegg_r77q{{tp}}hr_vs_wt{{tp}}hr_upDEG_sigPathways.csv", \
        tp = [4, 8, 12, 24, 72]),
        kegg_res_down_vs_mock = expand(f"{berges_dir}/results/r/tables/kegg_{{contr_name}}_downDEG_sigPathways.csv", \
        contr_name = ["wt4hr_vs_mock72hr", "wt8hr_vs_mock72hr", "wt12hr_vs_mock72hr", \
        "wt24hr_vs_mock72hr", "wt72hr_vs_mock72hr", \
        "r77q4hr_vs_mock72hr", "r77q8hr_vs_mock72hr", "r77q12hr_vs_mock72hr", \
        "r77q24hr_vs_mock72hr", "r77q72hr_vs_mock72hr"]),
        kegg_res_down_r77qvswt = expand(f"{berges_dir}/results/r/tables/kegg_r77q{{tp}}hr_vs_wt{{tp}}hr_downDEG_sigPathways.csv", \
        tp = [4, 8, 12, 24, 72]),
        kegg_diagrams = directory(f"{berges_dir}/results/r/figures/r77q72hr_vs_wt72hr_kegg_pathway_diagrams")
    shell:
        """
        {input.r_script}
        """

############# RUN COMPLETE WORKFLOW #############
rule run_workflow:
    input:
        rules.MultiQC_all_fastqcs.output,
        rules.DESeq2_salmon_DE_analysis.output,
        rules.plot_deg_barplots.output,
        rules.functional_enrichment_analysis.output