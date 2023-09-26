# TFAM
Transcription factor family annotation and mining tool

### <ins>Usage:</ins>
```
perl TFAM.pl
```
</p>
Note, the blat2hints.pl script must be in your working directory for augustus to work!

<br>


### <ins>Dependencies:</ins>
<ol type="1">
  
#### <li>HMMER:</li>
To install hmmer and easel miniapps, follow instructions from hmmer manual (pgs 17-18): <p>
http://eddylab.org/software/hmmer/Userguide.pdf

#### <li>BLAST:</li>
<b>Quick install with root privilages:</b>
```
sudo apt-get update
sudo apt-get -y install ncbi-blast+
```
<b>Install blast from source:</b> <p>
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

#### <li>AUGUSTUS:</li>
<b>Quick install with root privilages:</b>
```
sudo apt-get update
sudo apt-get install augustus
```
<b>Install AUGUSTUS from source </b><p>
Augustus dependencies: https://github.com/Gaius-Augustus/Augustus/blob/master/docs/INSTALL.md <p>
Build augustus: https://github.com/Gaius-Augustus/Augustus <p>

Make sure to set your AUGUSTUS_CONFIG_PATH variable by appending the following to your ~/.bashrc file: <p>
```
export AUGUSTUS_CONFIG_PATH=/my_path_to_AUGUSTUS/Augustus/config/    #where my_path_to_AUGUSTUS is dependent on where you cloned the augustus repo
```
#### <li>BLAT:</li>
Augustus requires blat to generate hints. Follow instructions here: <p>
https://bioinformaticsreview.com/20200822/installing-blat-a-pairwise-alignment-tool-on-ubuntu/ 
  </ol>
  
</p>
<br>

### <ins>Pipeline:</ins>
<p align="center">
<img src = "pipeline_image/TFAM_pipeline.png" width=2000>
</p>

<br>

### <ins> Input files, options and parameters:</ins>

<br>
<ol type="1">
<b>
  
<li> Input files and profile HMMs: </p>    
</b> 

Enter the name of your protein PFAM seed alignment here:
```
$pfam_seed = "protein_alignment.aln"; #Line 17
```
Enter the name of your nucleotide alignment here:
```
$nuc_alignmnent = "nucletide_alignment.aln"; #Line 18
```
If augustus prediction ($predict_new_hits) option is on, enter the name of your reference file here:
```
$reference_file = "reference_transcription_factors_mrna.fa"; #Line 21
```
Profile HMMs will be automatically built from your alignment files. To specify the names of these profile HMMs, the following variables need to be set:
```
my $phmm_profile = "my_name.hmm"; #Please note that this must end in ".hmm". #Line 24
my $nhmm_profile = "my_name.hmm"; #Please note that this must end in ".hmm". #Line 25
```


<b>
</ul>
<br>
  
<li> Options: </p>
</b>

If NCBI RefSeq annotations are available for your genome, set this variable to "yes". Set as "no" if no annotations are available, and you wish to mine the assembly only:
```
$annotation_available = "yes"; #Line 33
```
If yes, script will use a list of species to download assembly and annotation files:
```
$automate_download = "yes"; #Line 36
```
Enter the name of your species list file here:
```
$species_list = "species.txt"; #Line 37
```
If “yes”, the default e-value (1e-5) will be used for <b> protein hmmer </b>. For custom evalue, set this variable to "no":
```
$default_phmmer_evalue = "yes"; #Line 40
```

If $default_phmmer_evalue is "no", enter your custom evalue for <b> protein hmmer </b> here:
```
$phmmer_evalue = "1e-5"; #Line 41
```
If “yes”, the default e-value (1e-5) will be used for <b> nucleotide hmmer </b>. For custom evalue, set this variable to "no":
```
$default_nhmmer_evalue = "yes"; #Line 44
```

If $default_nhmmer_evalue is "no", enter your custom evalue for <b> nucleotide hmmer </b>  here:	
```
$nhmmer_evalue = "1e-5"; #Line 45
```
If you want to predict any new, unannotated, hits with augustus, set this to "yes". If augustus is not installed, set this as "no":
```
$predict_new_hits = "yes"; #Line 48
```
If using augustus, this is the species that augustus is trained on. Set your closely related species here:
```
$augustus_species = "arabidopsis"; #Line 49
```
If using augustus, this is the minimum identity required for alignment with a reference receptor to be used to generate prediction hints:
```
$minidentity = 60; #Line 50
```
If using augustus, this is the number of sequences from the reference file used to generate hints. This can either be set to a number greater than 0, or to "all" if you wish to use the entire reference file.
```
my $number_hints = "all";  #Line 51
```
If using augsutus, you can chose to append the mined NCBI sequences from each query species to the reference file to guide gene prediction. If you wish to do this, set the following variable to "yes". Otherwise this should be set to "no".
```
my $append_query = "no"; #Line 52
```
If new hits are predicted with augustus, they will be labelled with a prefix set using the following variable. For example, if this is set to "Hit", new opredictions will be labelled as "Hit1, Hit2 ...etc.". This can be adjusted to suit the use-case. 
```
my $hit_prefix = "Hit"; #Line 53
```
New augutsus predictions are checked to ensure that the hmmer identified region is retained in the prediction. The following parameter specifies the percentage of the hmmer identified region must be retained in the prediction to be considered valid. This number should be between 0 and 1.
```
my $domain_cover_threshold = 0.9; #Line 54
```
Augustus predictions are also scanned again using HMMER to ensure the domain is present. This HMM filter can use either the protein profile HMM or the nucleotide profile HMM. To select either option, specify the below variable as "protein" or "nucleotide".
```
my $hmm_filter_type = "protein"; #Line 55
```
For each novel hit identified with nhmmer on the genome assembly, the nucleotide region upstream and downstream of the hit are retrieved and fed into augustus for gene prediction. To specify the amount of nucleotides added to the 3' and 5' ends of the hit prior to prediction, please  specify the following variables: </p>

Add X nucleotides to end of sequence (3' end): 

```
$nhmmer_plus = 20000; #Line 58
```
Add X nucleotides to start of sequence (5' end):
```
$nhmmer_minus = 5000; #Line 59
```

Running nhmmer on the whole genome assembly can be an intensive task, requiring long run-times. To speed up the process for large genomes, an nhmmer database can be generated for the assembly. While, this dramatically speeds up run-times, it reduces sensitivity slightly - and hence we reccomend that this option is only switched on for large genomes. To specify this option, please set the following variable:
```
my $nhmmer_genome_database = "no"; #Line 62
```
If yes, all cds seqs with in-frame stop codons, or below threshold length, will be annotated as pseudogenes:
```
$pseudogene_check = "yes"; #Line 65
```
Coding sequences below this length are considered pseudogenes (nucleotide length):
```
$pseudogene_length = 300; #Line 66 
```

If you want to remove duplicates which may arise due to assembly error, set the following variable to "yes". This option will create a percent identity matrix to identify potential duplicates. The dupplicate on the largest contig is retained.
```
my $remove_duplicates = "no"; #Line 69
```
To specify the threshold percentage identity for which genes are considered duplicates, please set the following variable. Note that this is a percentage and should be a number between 0 and 1.
```
my $duplicate_threshold = 0.9; #Line 71
```
When identifying duplicates, the pipeline can make use of two distinct algoirithms - "pairwise" or "cluster". In the pairwise algorithm, duplicate pairs are identified as mutual best scoring hits in the percent identity matrix. Note that more than 2 members may exist in a pair, if each member shares the maximum identity score. Mutual best scores are only considered pairs if they exceed the $duplicate_threshold set above. The member in each pair which is located on the longest contig is retained. In the "cluster" algorithm, genes which share greater percent identity than the $duplicate_threshold are combined into clusters. The member in each cluster which is located on the longest contig is retained. To specify whether you want to use the "pairwise" or "cluster" algorithms, please set the below variable:
```
my $duplicate_type = "cluster"; #pairwise or cluster. 
```

<b>
</ul>
<br>






