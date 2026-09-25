In this example, GENE-FAM will be run on an unannotated GenBank cannabis genome: Abacus strain (GCA_025232715.1). 

As this is not a RefSeq genome, you have to manually download and unzip the assembly. This can be done through the NCBI database or using the following commands:

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/232/715/GCA_025232715.1_Csat_AbacusV2/GCA_025232715.1_Csat_AbacusV2_genomic.fna.gz

gunzip GCA_025232715.1_Csat_AbacusV2_genomic.fna.gz

Alternatively 'bash download_genome.sh' will work for this example.

Once the genome is in your directory, run the following command:

bash run-gene-fam.sh

or run the command directly:

perl GENE-FAM.pl \
    --prot-alignment PF00319_seed.txt \
    --nuc-alignment MADS_nhmmer_alignment.fa \
    --reference MADS_reference_file.fa \
    --annotation-available no \
    --automate-download no


To run GENE-FAM using Docker, you can either use the provided Bash script:

bash run-gene-fam-docker.sh

or run the Docker command directly:

docker run --rm \
    -v "$PWD:/data" \
    -w /data \
    louiseryan314/gene-fam \
    --prot-alignment /data/PF00319_seed.txt \
    --nuc-alignment /data/MADS_nhmmer_alignment.fa \
    --reference /data/MADS_reference_file.fa \
    --annotation-available no \
    --automate-download no