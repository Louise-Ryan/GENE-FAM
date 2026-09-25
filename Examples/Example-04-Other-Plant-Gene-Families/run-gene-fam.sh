#!/bin/bash

set -e

echo "Running GENE-FAM for MADS-box genes..."

# Run GENE-FAM for MADS-box genes
perl GENE-FAM.pl \
    --prot-alignment PF00319_seed.txt \
    --nuc-alignment MADS_nhmmer_alignment.fa \
    --reference MADS_reference_file.fa

for output_dir in GCF*_outfiles; do
    mv "$output_dir" "MADS_${output_dir}"
done

echo "Running GENE-FAM for B3 genes..."

# Run GENE-FAM for B3 genes
perl GENE-FAM.pl \
    --prot-alignment B3-PF02362.26-Protein-alignment.seed \
    --nuc-alignment B3-Nucleotide-alignment.fa \
    --reference B3-Reference-File.fa \
    --automate-download no

for output_dir in GCF*_outfiles; do
    mv "$output_dir" "B3_${output_dir}"
done

echo "Summary of MADS-box and B3 genes across Arabidopsis and tomato:"
bash get-summary.sh