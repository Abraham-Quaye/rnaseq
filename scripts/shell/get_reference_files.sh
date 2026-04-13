#!/usr/bin/env bash
echo "Downloading GRCm39 genome fasta"
wget -O raw_files/genome_files/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz &&

echo "Downloading GRCm39 GTF annotation file"
wget -O raw_files/annotations/Mus_musculus.GRCm39.gtf.gz https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz

echo "Downloading GRCm39 transcriptome fasta file"
wget -O raw_files/genome_files/Mus_musculus.GRCm39.cdna.all.fa.gz https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

if [ $? -eq 0 ]; then
    echo "Script executed successfully!"
else
    echo "Script failed."
fi