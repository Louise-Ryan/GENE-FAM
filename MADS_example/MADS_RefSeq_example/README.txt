In this example, TFAM will be run on Cannabis, Apple and Arabidopsis looking for MADS box genes.

As these are all RefSeq genomes, with annotation files available, the genomes can be automatically downloaded with TFAM.

The species of interest are specified in the species.txt file.

Run the following command:
perl TFAM.pl

If everything is installed properly, this should work fine!

You can ignore "smartmatch warnings" if they pop up.

As annotation files are available for each species, the following parameters should be set as below:
   - $annotation_available = "yes"
   - $cds_available = "yes"
   - $automate_download = "yes"
   - $species_list = "species.txt"
