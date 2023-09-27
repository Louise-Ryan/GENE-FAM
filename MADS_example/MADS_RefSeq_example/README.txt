In this example, GENE-FAM will be run on Apple and Arabidopsis to mine MADS box genes.

As these are all RefSeq genomes, with annotation files available, the genomes can be automatically downloaded using GENE-FAM.

The species of interest are specified in the species.txt file.

Run the following command:
perl GENE-FAM-RefSeq-example.pl

If everything is installed properly, this should work fine!

As annotation files are available for each species, the following parameters should be set as below:
   - $annotation_available = "yes"
   - $automate_download = "yes"
   - $species_list = "species.txt"
