#!/bin/bash

pperl GENE-FAM.pl \
    --prot-alignment TAAR-Protein-alignment.fa \
    --nuc-alignment TAAR-Nucleotide-alignment.fa \
    --reference TAAR-Reference-File.fa \
    --annotation-available yes \
    --default-phmmer-evalue no \
    --phmmer-evalue 1e-100 \
    --default-nhmmer-evalue no \
    --nhmmer-evalue 1e-100 \
    In this example, GENE-FAM will be run on the Puma concolor and Orcinus orca RefSeq genomes, with the aim of mining Trace Amine-Associated Receptor (TAAR) genes.

As these are RefSeq genomes with annotation files available, the genomes can be automatically downloaded using GENE-FAM.

The species of interest is specified in the `species.txt` file.

TAARs are members of the G protein-coupled receptor (GPCR) family and share sequence similarity with other GPCRs. 
To reduce the inclusion of non-TAAR GPCRs, a stringent E-value threshold is required for both phmmer and nhmmer searches. 
For this example, the default E-value settings need to be disabled and an E-value of 1e-100 specified for both searches.

Furthermore, as we are now investigating mammals, the AUGUSTUS training species needs to be changed from the default 'Arabidopsis' to 'human'.


To run GENE-FAM, you can either use the provided Bash script:

bash run-gene-fam.sh

or run the command directly:

perl GENE-FAM.pl \
    --prot-alignment TAAR-Protein-alignment.fa \
    --nuc-alignment TAAR-Nucleotide-alignment.fa \
    --reference TAAR-Reference-File.fa \
    --annotation-available yes \
    --default-phmmer-evalue no \
    --phmmer-evalue 1e-100 \
    --default-nhmmer-evalue no \
    --nhmmer-evalue 1e-100 \
    --augustus-species arabidopsis

To run GENE-FAM using Docker, you can either use the provided Bash script:

bash run-gene-fam-docker.sh

or run the Docker command directly:

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
