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

### <ins>Options and parameters:</ins>
<ol type="1">
<b>
<li> Options: </p>
</b>

If NCBI RefSeq annotations are available for your genome, set this variable to "yes". Set as "no" if no annotations are available, and you wish to mine the assembly only:
```
$annotation_available = "yes";
```
Keep this as "yes" if cds files are available. Set as "no" if you only want to mine the mRNA files:
```
$cds_available = "yes"; 
```
If you want to predict any new, unannotated, hits with augustus, set this to "yes". If augustus is not installed, set this as "no":
```
$predict_new_hits = "yes";
```
If using augustus, this is the species that augustus is trained on. Set your closely related species here:
```
$augustus_species = "arabidopsis";
```
If using augustus, this is the minimum identity required for alignment with a reference receptor to be used to generate prediction hints:
```
$minidentity = 60;
```
If yes, all cds seqs with in-frame stop codons, or below threshold length, will be annotated as pseudogenes:
```
$pseudogene_check = "yes";
```
Coding sequences below this length are considered pseudogenes (nucleotide length):
```
$pseudogene_length = 300 ;
```

<b>
</ul>
<br>
<li> Input files: </p>    
</b> 

Enter the name of your protein PFAM seed alignment here:
```
$pfam_seed = "PF00319_seed.txt";
```
Enter the name of your nucleotide alignment here:
```
$nuc_alignmnent = "nucletide_alignment";
```
If augustus prediction ($predict_new_hits) option is on, enter the name of your reference file here:
```
$reference_file = "reference_transcription_factors_mrna.fa";
```
If yes, script will use a list of species to download assembly and annotation files:
```
$automate_download = "yes";
```
Enter the name of your species list file here:
```
$species_list = "species.txt";
```
<b>
</ul>
<br>
<li> Parameters: </p>    
</b>

If “yes”, the default e-value (1e-5) will be used for <b> protein hmmer </b>. For custom evalue, set this variable to "no":
```
$default_phmmer_evalue = "yes";
```

If $default_phmmer_evalue is "no", enter your custom evalue for <b> protein hmmer </b> here:
```
$phmmer_evalue = "1e-5";
```
If “yes”, the default e-value (1e-5) will be used for <b> nucleotide hmmer </b>. For custom evalue, set this variable to "no":
```
$default_nhmmer_evalue = "yes";
```

If $default_nhmmer_evalue is "no", enter your custom evalue for <b> nucleotide hmmer </b>  here:	
```
$nhmmer_evalue = "1e-5";
```

For each novel, unannotated hit, X nucleotides are taken from the 5' and 3' ends of the hmmer hit, and fed into augustus. </p>

Add X nucleotides to end of sequence (3' end) (cannot exceed 990,000 nt): 

```
$nhmmer_plus = 20000;
```
Add X nucleotides to start of sequence (5' end) (cannot exceed 990,000 nt):
```
$nhmmer_minus = 5000;
```
<b>
</ul>
<br>
<li> Profile hmm names: </p>    
</b> 

Set the name of your phmmer hmm profile. Recommend naming this the same as your protein PFAM seed alignment with the “.hmm” extension:
```
$phmm_profile = "PF00319.hmm";
```
Set the name of your nhmmer hmm profile. Recommend naming this the same as your nucleotide seed alignment with the “.hmm” extension:
```
$nhmm_profile = "nucleotide_alignment.hmm";
```

<br>
<b>
<li> Pseudogene vs Fuctional:</p>    
</b> 

If prediction is shorter than $pseudogene_length (set below), genes will receive the above annotation in their sequence header:
```
$annotation_short = "pseudogene_short";
```

Genes with cds lengths lower than this number be annotated as $pseudogene_short (above) regardless of conditions (1-6):
```
$pseudogene_length = “300”;
```






