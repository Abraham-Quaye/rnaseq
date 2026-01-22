#!/usr/bin/env bash
echo "Downloading GRCh38 genome fasta"
wget -O raw_files/genome_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz &&

echo "Downloading GRCh38 GTF annotation file"
wget -O raw_files/annotations/Homo_sapiens.GRCh38.115.gtf.gz https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

echo "Downloading GRCh38 transcriptome fasta file"
wget -O raw_files/genome_files/Homo_sapiens.GRCh38.cdna.all.fa.gz https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

if [ $? -eq 0 ]; then
    echo "Script executed successfully!"
else
    echo "Script failed."
fi