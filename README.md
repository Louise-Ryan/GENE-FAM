# TFAM
Transcription factor family annotation and mining tool

### <ins>Usage:</ins>
<b>perl TFAM.pl</b>
</p>
Note, the blat2hints.pl script must be in your working directory for augustus to work!

<br>


### <ins>Dependencies:</ins>
<ol type="1">
  
#### <li>HMMER:</li>
To install hmmer and easel miniapps, follow instructions from hmmer manual (pgs 17-18): <p>
http://eddylab.org/software/hmmer/Userguide.pdf

#### <li>BLAST:</li>
Download latest version of BLAST here: <p>
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

#### <li>AUGUSTUS:</li>
Download AUGUSTUS here: <p>
https://github.com/Gaius-Augustus/Augustus 

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
<ul>
<li> <b> $annotation_available = "yes"</b> </p>If NCBI RefSeq annotations are available for your genome, set this variable to "yes". Set as "no" if no annotations are available, and you wish to mine the assembly only. </li>
</p>
<li> <b> $cds_available = "yes" = "yes"</b> </p> Keep this as "yes" if cds files are available. Set as "no" if you only want to mine the mRNA files.
</p>
<li> <b> $predict_new_hits = "yes"</b> </p> If you want to predict any new, unannotated, hits with augustus, set this to "yes". If augustus is not installed, set this as "no". </li>
</p>
<li> <b> $augustus_species = "arabidopsis"</b> </p> If using augustus, set your closely related species here. This is the species that augustus is trained on. </li>
</p>
<li> <b>$minidentity = 60 </b> </p> If using augustus, this is the minimum identity required for alignment with a reference receptor to be used to generate prediction hints. </li>
</p>
<li> <b> $pseudogene_check = "yes" </b> <p> If yes, all cds seqs with in-frame stop codons, or below threshold length, will be annotated as pseudogenes. </li>
</p>
<li> <b> $pseudogene_length = 300 </b></p> Coding sequences below this length are considered pseudogenes (nucleotide length). </li>
</p>
<b>
</ul>
<br>
<li> Input files: </p>    
</b> 
<ul> <b> $pfam_seed = "PF00319_seed.txt" </b> <p> Enter the name of your protein PFAM seed alignment here. </li>
</p>
<li> <b> $nuc_alignmnent = "nucletide_alignment" </b> <p> Enter the name of your nucleotide alignment here.
</p>
<li> <b> $reference_file = "reference_transcription_factors_mrna.fa" </b> <p>  If augustus prediction ($predict_new_hits) option is on, enter the name of your reference file here.
</p>
<li> <b> $automate_download = "yes" </b> <p> If yes, script will use a list of species to download assembly and annotation files.
</p>
<li> <b> $species_list = "species.txt" </b> <p> Enter the name of your species list file here.
</p>
<b>
</ul>
<br>
<li> Parameters: </p>    
</b> 
<ul> 
<li> <b> $default_phmmer_evalue = "yes" </b> <p> If “yes”, the default e-value will be used (1e-5). If “no”, you can enter a custom evalue in the 
</p>
<li> <b> $phmmer_evalue = "1e-5" </b> <p> If $default_phmmer_evalue is "no", enter your custom evalue here.
</p>
<li> <b> $default_nhmmer_evalue = "yes" </b> <p> If “yes”, the default e-value will be used (1e-5). If “no”, you can enter a custom evalue in the $nhmmer_evalue variable below.
</p>  
<li> <b> $nhmmer_evalue = "1e-5" </b> <p> If $default_nhmmer_evalue is "no", enter your custom evalue here.	
</p>
<li> <b> $nhmmer_plus = 20000 </b> <p> Add X nucleotides to end of sequence (3' end) (cannot exceed 990,000 nt) 
</p>
<li> <b> $nhmmer_minus = 5000 </b> <p> Add X nucleotides to start of sequence (5' end) (cannot exceed 990,000 nt)
</p>
<b>
</ul>
<br>
<li> Profile hmm names: </p>    
</b> 
<ul> 
<li> <b> $phmm_profile = "PF00319.hmm" </b> This is the name of your phmmer hmm profile. Recommend naming this the same as your protein PFAM seed alignment with the “.hmm” extension.
</p>
<li> <b> $nhmm_profile = "nucleotide_alignment.hmm" </b> This is the name of your nhmmer hmm profile. Recommend naming this the same as your nucleotide seed alignment with the “.hmm” extension.
</p>
<b>
</ul>
<br>
<li> Pseudogene vs Fuctional:</p>    
</b> 
<ul>
<li> <b> $annotation_short = "pseudogene_short" </b> </p> If prediction is shorter than $pseudogene_length (set below), genes will receive the above annotation in their sequence header. 
</p>
<li> <b> $pseudogene_length = “300”. </b> </p> Genes with cds lengths lower than this number be annotated as $pseudogene_short (above) regardless of conditions (1-6).
</p>




