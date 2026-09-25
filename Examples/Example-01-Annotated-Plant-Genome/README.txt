In this example, GENE-FAM will be run on the *Arabidopsis thaliana* RefSeq genome, with the aim of mining MADS-box genes.

As this is a RefSeq genome with annotation files available, the genome can be automatically downloaded using GENE-FAM.

The species of interest is specified in the `species.txt` file.

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