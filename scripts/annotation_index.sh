#!/bin/bash

# Directories
raw_dir="../data/raw"
ref_genome_dir="../data/ref_genome"
results_dir="../results"
log_dir="../logs"

mkdir -p $results_dir/index
mkdir -p $log_dir

echo "indexing the reference genome"
STAR \
--runMode genomeGenerate \
--runThreadN 4 \
--genomeDir $results_dir/index \
--genomeFastaFiles $ref_genome_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile $ref_genome_dir/Homo_sapiens.GRCh38.112.gtf \
--genomeSAsparseD 3 \
--genomeSAindexNbases 12 > $log_dir/index.log 2>&1
echo "indexing completed"