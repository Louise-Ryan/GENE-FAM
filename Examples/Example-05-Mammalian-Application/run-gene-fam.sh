#!/bin/bash

perl GENE-FAM.pl \
    --prot-alignment TAAR-Protein-alignment.fa \
    --nuc-alignment TAAR-Nucleotide-alignment.fa \
    --reference TAAR-Reference-File.fa \
    --annotation-available yes \
    --default-phmmer-evalue no \
    --phmmer-evalue 1e-100 \
    --default-nhmmer-evalue no \
    --nhmmer-evalue 1e-100 \
    --augustus-species human
