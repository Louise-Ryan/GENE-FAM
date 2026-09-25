#!/bin/bash

set -e


echo "Running GENE-FAM for MADS-box genes..."

# Run GENE-FAM for MADS-box genes
docker run --rm \
    -v "$PWD:/data" \
    -w /data \
    louiseryan314/gene-fam \
    --prot-alignment /data/PF00319_seed.txt \
    --nuc-alignment /data/MADS_nhmmer_alignment.fa \
    --reference /data/MADS_reference_file.fa \
    --automate-download yes

for output_dir in GCF*_outfiles; do
    mv "$output_dir" "MADS_${output_dir}"
done


echo "Running GENE-FAM for B3 genes..."

# Run GENE-FAM for B3 genes
docker run --rm \
    -v "$PWD:/data" \
    -w /data \
    louiseryan314/gene-fam \
    --prot-alignment /data/B3-PF02362.26-Protein-alignment.seed \
    --nuc-alignment /data/B3-Nucleotide-alignment.fa \
    --reference /data/B3-Reference-File.fa \
    --automate-download no

for output_dir in GCF*_outfiles; do
    mv "$output_dir" "B3_${output_dir}"
done


echo "Summary of MADS-box and B3 genes across Arabidopsis and tomato:"
bash get-summary.sh