In this example, TFAM will be run on an unannotated cannabis genome: Abacus strain (GCA_025232715.1). 

As this is not a RefSeq genome, you have to manually download and unzip the assembly. This can be done through the NCBI database or using the following commands:

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/232/715/GCA_025232715.1_Csat_AbacusV2/GCA_025232715.1_Csat_AbacusV2_genomic.fna.gz

gunzip GCA_025232715.1_Csat_AbacusV2_genomic.fna.gz

(alternatively bash download_genome.sh will work)

Once the genome is in your directory, run the following command:
perl TFAM.pl

You can ignore "smartmatch warnings" if they pop up.

Note that as the genome is unannotated, and hence no annotation files are available, the following variables were changed from default in the TFAM.pl script:
      - $annotation_available = "no"
      - $cds_available = "no"
      - $automate_download = "no"
      
      
