In this example, GENE-FAM will be run on four RefSeq plant genomes, with the aim of mining MADS-box genes and demonstrating the cross-species application of the pipeline.

As these are RefSeq genomes with annotation files available, the genomes can be automatically downloaded using GENE-FAM. 
As such, users only need to provide a list of species, making cross-species comparisons easy and efficient with minimal file preparation.

The species of interest are specified in the species.txt file:
Brachypodium distachyon
Selaginella moellendorffii
Amborella trichopoda
Malus domestica

To run GENE-FAM, you can either use the provided Bash script:

bash run-gene-fam.sh

or run the command directly:

perl GENE-FAM.pl \
    --prot-alignment PF00319_seed.txt \
    --nuc-alignment MADS_nhmmer_alignment.fa \
    --reference MADS_reference_file.fa
Running GENE-FAM using Docker

To run GENE-FAM using Docker, you can either use the provided Bash script:

bash run-gene-fam-docker.sh

or run the Docker command directly:

docker run --rm \
    -v "$PWD:/data" \
    -w /data \
    louiseryan314/gene-fam \
    --prot-alignment /data/PF00319_seed.txt \
    --nuc-alignment /data/MADS_nhmmer_alignment.fa \
    --reference /data/MADS_reference_file.fa

When complete, run get-summary.sh, to get a summary of MADS box genes mined across each species.
