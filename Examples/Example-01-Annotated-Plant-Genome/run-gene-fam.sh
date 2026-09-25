#!/bin/bash

perl GENE-FAM.pl \
    --prot-alignment PF00319_seed.txt \
    --nuc-alignment MADS_nhmmer_alignment.fa \
    --reference MADS_reference_file.fa
