#!/bin/bash

# Directories
results_dir="../results/Multi_QC"

mkdir -p $results_dir

echo "running Multi QC"
multiqc -n multiqc_report --module star --outdir $results_dir *.out
echo "quality control completed"
