#!/bin/bash

# Directories
raw_dir="../data/raw"
ref_genome_dir="../data/ref_genome"
results_dir="../results"
log_dir="../logs"

# Make sure output directories exist
mkdir -p $results_dir
mkdir -p $log_dir

# Running STARsolo for demultiplexing and quantification
echo "Running STARsolo for demultiplexing and quantification"
STAR \
--genomeDir $results_dir/index \
--readFilesIn $raw_dir/subset/L001_R2_001.fastq.gz,$raw_dir/subset/L002_R2_001.fastq.gz $raw_dir/subset/L001_R1_001.fastq.gz,$raw_dir/subset/L002_R1_001.fastq.gz \
--readFilesCommand zcat \
--soloType CB_UMI_Simple \
--runThreadN 4 \
--soloCBwhitelist $raw_dir/3M-february-2018.txt \
--soloUMIdedup 1MM_CR \
--soloCBmatchWLtype 1MM_multi \
--soloCBlen 16 \
--soloUMIlen 12 \
--soloStrand Forward \
--outSAMtype None \
--outFileNamePrefix STARsolo_ \
--sjdbGTFfile $ref_genome_dir/Homo_sapiens.GRCh38.112.gtf \
--soloUMIfiltering - \
--soloCellFilter None \
--soloFeatures Gene \
--soloOutFileNames ../results/STARsolo features.tsv barcodes.tsv matrix.mtx \
--runDirPerm All_RWX 


echo "STARsolo demultiplexing and quantification completed"
