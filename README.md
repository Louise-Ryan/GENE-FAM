# TFAM
Transcription factor family annotation and mining tool


### <ins>Overview:</ins>

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

### <ins>Options and parameters:</ins>
<ol type="1">
<b>
<li> Options: </p>
</b>
<ul>
<li> <b> $annotation_available = "yes"</b> </p>If NCBI RefSeq annotations are available for your genome, set this variable to "yes". Set as "no" if no annotations are available, and you wish to mine the assembly only.
</p>
<li> <b> $cds_available = "yes" = "yes"</b> </p> Keep this as "yes" if cds files are available. Set as "no" if you only want to mine the mRNA files.
</p>
<li> <b> $predict_new_hits = "yes"</b> </p> If you want to predict any new, unannotated, hits with augustus, set this to "yes". If augustus is not installed, set this as "no".
</p>
<li> <b> $augustus_species = "arabidopsis"</b> </p> If using augustus, set your closely related species here. This is the species that augustus is trained on.
</p>
<li> <b>$minidentity = 60 </b> </p> If using augustus, this is the minimum identity required for alignment with a reference receptor to be used to generate prediction hints.
</p>
<li> <b> $pseudogene_check = "yes" </b></p> If yes, all cds seqs with in-frame stop codons, or below threshold length, will be annotated as pseudogenes.
</p>
<li> <b> $pseudogene_length = 300 </b></p> Coding sequences below this length are considered pseudogenes (nucleotide length).
</p>
<b>
<li> Input files: </p>
</b>
<ul>

