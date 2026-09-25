#!/bin/bash

docker run --rm \
    -v "$PWD:/data" \
    -w /data \
    louiseryan314/gene-fam \
    --prot-alignment /data/TAAR-Protein-alignment.fa \
    --nuc-alignment /data/TAAR-Nucleotide-alignment.fa\
    --reference /data/TAAR-Reference-File.fa \
    --annotation-available yes \
    --default-phmmer-evalue no \
    --phmmer-evalue 1e-100 \
    --default-nhmmer-evalue no \
    --nhmmer-evalue 1e-100 \
    --augustus-species human

    
