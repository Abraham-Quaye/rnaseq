##########  GENERATE hg38 FASTA FROM R PACKAGE ####################
rule make_hg38_fasta:
    input:
        "scripts/r_code/make_hg38_fasta.R"
    output:
        "raw_files/genome_files/hg38_genome.fa"
    shell:
        "{input}"

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

# #################### BUILD hg38 GENOMIC INDEX FOR MAPPING WITH RSubread ######
rule build_rsubread_genomic_index:
    input:
        hg38 = rules.make_hg38_fasta.output,
        grch38 = rules.download_reference_files.output.genome,
        r_script = "scripts/r_code/build_genome_index.R"
    output:
        hg38 = directory("raw_files/rsubread_hg38_genome_index"),
        grch38 = directory("raw_files/rsubread_GRCh38_genome_index")
    shell:
        """
        mkdir -p {output.hg38} {output.grch38}
        {input.r_script}
        """

################ BUILD GENOMIC INDICES FOR MAPPING WITH HISAT2 ###############
rule build_hisat_genomic_index:
    input:
        grch38_genome = rules.download_reference_files.output.genome,
        hg38_genome = rules.make_hg38_fasta.output
    output:
        expand("raw_files/hisat_GRCh38_genome_index/grch38_tran.{n}.ht2", n = range(1,9)),
        expand("raw_files/hisat_hg38_genome_index/hg38_tran.{n}.ht2", n = range(1,9))
    shell:
        """
        mkdir -p raw_files/hisat_GRCh38_genome_index raw_files/hisat_hg38_genome_index

        echo "Building GRCh38 genomic index for Hisat2..."
        hisat2-build -p 15 {input.grch38_genome} raw_files/hisat_GRCh38_genome_index/grch38_tran

        echo "Building hg38 genomic index for Hisat2..."
        hisat2-build -p 15 {input.hg38_genome} raw_files/hisat_hg38_genome_index/hg38_tran
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

############# RUN COMPLETE WORKFLOW #############
rule run_workflow:
    input:
        rules.build_rsubread_genomic_index.output,
        rules.build_hisat_genomic_index.output,
        rules.build_star_genomic_index.output