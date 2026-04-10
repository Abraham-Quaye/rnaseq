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

proj_dir = "/Users/abrahamquaye/bm_fn_rnseq"

############### TRIM READS WITH TRIM GALORE ####################
rule trim_fastq_files:
    input:
        raw_fastqs = rules.concatenate_same_strand_reads.output,
        script = "scripts/shell/trim_reads.sh"
    output:
        t_fastqs = expand(f"{proj_dir}/results/trimmedReads/{{sample}}{{snum}}_R{{strand}}_001.fastq.gz", \
        sample = ["BM", "FN"], snum = range(1, 5), strand = [1, 2])
    shell:
        "{input.script}"

############### QC TRIMMED READS WITH FASTQC ####################
rule FastQC_trimmed_reads:
    input:
        t_fastqs = rules.trim_fastq_files.output.t_fastqs
    output:
        t_fastqc = expand(f"{proj_dir}/results/fastqc/{{sample}}{{snum}}_R{{strand}}_001.fastq.gz_fastqc.html", \
        sample = ["BM", "FN"], snum = range(1, 5), strand = [1, 2]),
        fastqc_dir = directory(f"{proj_dir}/results/fastqc")
    shell:
        """
        mkdir -p {output.fastqc_dir}
        fastqc -t 16 --memory 1024 -o {output.fastqc_dir} {input.t_fastqs}
        """

#################### QC READS WITH MULTIQC ####################
mqc_dir = f"{proj_dir}/results/multiqc"

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

#################### QUANTIFY READS WITH SALMON ####################
rule quantify_reads_salmon:
    input:
        t_fastqs = rules.trim_fastq_files.output.t_fastqs,
        salmon_index = rules.build_salmon_genomic_index.output
    output:
        expand(f"{proj_dir}/results/salmon_quant/quant_siO_{{id}}/quant.sf", \
        id = ["GFP_1_S13", "GFP_2_S14", "GFP_3_S15", "MYOCD_1_S16", "MYOCD_2_S17", \
        "MYOCD_3_S18"])
    params:
        trim_dir = f"{proj_dir}/results/trimmedReads",
        salmon_dir = f"{proj_dir}/results/salmon_quant"
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
        tables = expand(f"{proj_dir}/results/r/tables/{{res_type}}_MYOCD_vs_GFP_DEGs.csv", \
        res_type = ["significant", "total"]),
        figs = expand(f"{proj_dir}/results/r/figures/{{fig_type}}_MYOCD_vs_GFP.pdf", \
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
        deg_barplots = f"{proj_dir}/results/r/figures/DEG_levels_barplot.pdf"
    shell:
        "{input.r_script}"

############ FUNCTIONAL ENRICHMENT ANALYSES OF DEGs #############
rule functional_enrichment_analysis:
    input:
        deg_files = rules.DESeq2_salmon_DE_analysis.output.tables,
        r_script = "scripts/r_code/enrichment_analysis.R",
        rscript2 = "scripts/r_code/enrichment_analysis_functions.R" 
    output:
        kegg_res = expand(f"{proj_dir}/results/r/tables/kegg_MYOCD_vs_GFP_{{reg}}DEG_sigPathways.csv", \
        reg = ["up", "down", "total"]),
        go_res = expand(f"{proj_dir}/results/r/tables/go_MYOCD_vs_GFP_{{reg}}DEG_sig.csv", \
        reg = ["up", "down", "total"]),
        go_figs = expand(f"{proj_dir}/results/r/figures/go_MYOCD_vs_GFP_{{reg}}DEG_sig_dotplot.pdf", \
        reg = ["up", "down", "total"]),
        kegg_figs = expand(f"{proj_dir}/results/r/figures/kegg_MYOCD_vs_GFP_{{reg}}DEG_sigPathways_dotplot.pdf", \
        reg = ["up", "down", "total"]),
        kegg_diagrams = directory(f"{proj_dir}/results/r/figures/kegg_pathway_diagrams")
    shell:
        """
        {input.r_script}
        mv hsa*.pathview.png {output.kegg_diagrams}
        """

############# RUN COMPLETE WORKFLOW #############
rule run_workflow:
    input:
        rules.MultiQC_all_fastqcs.output,
        rules.index_star_bam_files.output,
        rules.plot_deg_barplots.output,
        rules.functional_enrichment_analysis.output
