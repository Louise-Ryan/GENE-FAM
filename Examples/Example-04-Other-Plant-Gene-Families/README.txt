In this example, GENE-FAM will be run on two RefSeq plant genomes, Arabidopsis thaliana and Solanum lycopersicum, with the aim of mining both MADS-box and B3 genes.

As such, this example demonstrates the cross-species and cross-gene-family application of GENE-FAM.

As these are RefSeq genomes with annotation files available, the genomes can be automatically downloaded using GENE-FAM. 
The genomes are downloaded during the MADS-box analysis and then reused for the TCP analysis, avoiding the need to download the genomes twice.

The species of interest are specified in the `species.txt` file:
Arabidopsis thaliana
Solanum lycopersicum

The complete analysis can be run using the provided Bash script:
bash run-gene-fam.sh
This script runs GENE-FAM for both MADS-box and TCP genes and automatically generates a summary of the results.

The complete analysis can also be run using Docker using the provided Bash script:
bash run-gene-fam-docker.sh
The Docker script runs the same MADS-box and TCP analyses within the GENE-FAM Docker environment and automatically generates a summary of the results.
