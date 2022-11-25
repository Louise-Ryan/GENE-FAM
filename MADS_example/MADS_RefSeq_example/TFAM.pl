#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

###########################################################################
#USER PARAMETERS:                                                         #
###########################################################################
#Annotation files available?
my $annotation_available = "yes"; #If NCBI annotations are available for your genome set below variable to "yes". Set as "no" if no annotations are available, and you wish to mine the assembly only.
my $cds_available = "yes"; #keep this as yes if cds files are available. Set as "no" if you only want to mine the mrna files.
my $predict_new_hits = "yes"; #If you want to predict any new, unannotated, hits with augutus, set this to "yes". If augustus is not installed, keep this as "no".
my $augustus_species = "arabidopsis"; #If using augustus, set your closely related species here. This is the species that augustus is trained on.
my $minidentity = 60; #If using augutsus, this is the minimum identity required for alignment with a reference receptor to be used to generate prediction hints.
my $pseudogene_check = "yes"; #If yes, all cds seqs with in frame stop codons, or below threshold length, will be annotated as pseudogenes.
my $pseudogene_length = 300; #coding sequences below this length are considered pesuogenes (nucloetide length).

############################################################################
#Input Files
#phmmer and nhmmer alignments for transcription factor domain
my $pfam_seed = "PF00319_seed.txt"; #protein PFAM seed alignment
my $nuc_alignment = "MADS_nhmmer_alignment.fa"; #Nucleotide alignment
#Genome, Nucleotide and Protein file extension names:
my $genome_suffix = "genomic.fna"; #genome file (default for ncbi) #automatically downloaded if $automate_download is "yes".
my $nt_transcript_suffix = "rna.fna"; #nucleotide file (default for ncbi) #automatically downloaded if $automate_download is "yes".
my $prot_transcript_suffix = "protein.faa"; #protein file (default for ncbi) #automatically downloaded if $automate_download is "yes".
my $cds_suffix = "cds_from_genomic.fna"; #cds file (default for ncbi) #automatically downloaded if $automate_download is "yes".
#Augustus reference file:
my $reference_file = "MADS_reference_file.fa"; #if augustus option is on, enter reference file name here.
#Automate ncbi download?
#FOR REFSEQ GENOMES ONLY!#
#If yes, script will use a list of species to download assembly and annotation files
my $automate_download = "yes";
my $species_list = "species.txt"; #list species in species.txt file to download files for each species

###########################################################################
#User parameters and options:
#phmmer evalue threshold:
my $default_phmmer_evalue = "yes"; #yes: default; no: Use custom evalue (set $phmmer_evalue variable below)
my $phmmer_evalue = "1e-5"; #If $default_phmmer_evalue is "no", use this custom evalue
#nhmmer evalue threshold:
my $default_nhmmer_evalue = "yes"; #yes: default #no: user specified (use $nhmmer_evalue variable below)
my $nhmmer_evalue = "1e-5"; #If $default_nhmmer_evalue is "no", use this custom evalue
#Range for unannotated hits - plus or minus X nucleotides
my $nhmmer_plus = 20000; #Add X nucleotides to end of sequence (3' end) #cant exceed 990,000
my $nhmmer_minus = 5000; #Add X nucleotides to start of sequence (5' end) #cant exceed 990,000

###########################################################################
#output files
#hmm profile names
my $phmm_profile = "MADSp.hmm"; #nhmmer profile: use hmmer to build hmm profile from $pfam_seed
my $nhmm_profile = "MADSn.hmm"; #phmmer profile: use hmmer to build hmm profile from $nuc_alignment
#Output sequence file names
my $final_cds_nuc = "_cds_nuc.fa"; #final cds seq file - includes augustus predictions if turned on
my $final_cds_prot = "_cds_prot.fa"; #final cds seq file - includes augustus predictions if turned on
my $nucleotide_longest_transcripts="_longest_transcripts_rna.fa"; #Longest rna transcripts 
my $cds_nucleotide_seqfile="_longest_transcripts_cds_nucleotide.fa"; #Longest cds transcripts (nucleotide)
my $cds_protein_seqfile = "_longest_transcripts_cds_protein.fa"; #Longest cds transcripts (protein)
my $nhmmer_unnanotated_seqfile = "_unannotated_newhits_from_assembly.fa"; #Unannotated hits from nhmmer on assembly

###########################################################################
#Annotations: Functional vs Pseudogene
my $annotation_short = "pseudogene_short"; #If prediction is shorter than $pseudogene_length, gene will be annotated as pseudogene regardless of conditions (1-6).
my $annotation_1 = "functional"; #START codon && no in frame stop codons.........: ATG -----------
my $annotation_2 = "functional"; #no START codon && no stop codons in any frame..: ---------------
my $annotation_3 = "functional"; #no START codon && no in frame stop codons......: ---------------
my $annotation_4 = "pseudogene"; #START codon && in frame stop codon.............: ATG------TGA---
my $annotation_5 = "pseudogene"; #no START codon && stop codons in all frames....: ---------TGA---
my $annotation_6 = "pseudogene"; #no START codon && in frame stop codon..........: ---------TGA---

#############################################################################
#MAIN CODE                                                                  #
#############################################################################

#if automate download is yes, download files for species:
if ($automate_download eq "Yes" || $automate_download eq "yes"){
    print "\nAutomatically downloading genome files from ncbi for your query species ... \n\n";
    open(SPECIESFILE, $species_list);
    my @species = <SPECIESFILE>;
    close SPECIESFILE;
    chomp(@species);
    `wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt`; #RefSeq summary file
    foreach my $target(@species){
	print $target."\n";
	open(SUMMARY, "assembly_summary_refseq.txt");
	while(<SUMMARY>) { #loop through summary file
	    my $line=$_;
	    if ($line =~m/genome[\t]+[0-9]+[\t]+[0-9]+[\t]+$target[\t]/i){
		print $line."\n";
		if($line=~m/(ftp[\S]+)/){ #Store FTP link (FTP://...link) 
		    my $link=$1; 
		    if($link=~m/(GCF_+[0-9]+\.[0-9])/){ #Match and store genome acession number for summary files
			my $accession =$1." "; 
			if($link=~m/\/(GCF\_[\S]+)/){
			    my $accession = $1;
			    my $file_cds =$accession."_".$cds_suffix.".gz"; #Append file extension
			    my $file_genome = $accession."_".$genome_suffix.".gz";
			    my $file_rna = $accession."_".$nt_transcript_suffix.".gz";
			    my $file_prot = $accession."_".$prot_transcript_suffix.".gz";
			    my $final_link_cds = $link."/".$file_cds; #Append file extension and link to ftp link
			    `wget $final_link_cds`; #wget link to download the genome
			    sleep(4);
			    my $final_link_genome = $link."/".$file_genome;
			    `wget $final_link_genome`;
			    sleep(4);
			    my $final_link_rna = $link."/".$file_rna;
			    `wget $final_link_rna`;
			    sleep(4);
			    my $final_link_prot = $link."/".$file_prot;
			    `wget $final_link_prot`;
			    sleep(4);
			    `gunzip *gz`; #unzip file
			}
		    }
		}
	    }
	}
    }
}


#Read in genome, nucleotide and protein files: 
my @genomes=(<*$genome_suffix>); #read in multiple genome names
my @nucleotide_transcripts =(<*$nt_transcript_suffix>);
my @protein_transcripts = (<*$prot_transcript_suffix>);
my @cds_transcripts = ();
if ($cds_available eq "yes" || $cds_available eq "Yes"){
    @cds_transcripts = (<*$cds_suffix>);
}
#hmmer outfiles
my $phmmer_out = "_phmmer.out"; #genome name will be appedned to this file within code so file will look like: genome_phmmer.out
my $nhmmer_out ="_nhmmer.out"; #nhmmer outfile
#Build hmms:
if ($annotation_available eq "yes" || $annotation_available eq "Yes"){
    `hmmbuild $phmm_profile $pfam_seed`;
}
`hmmbuild $nhmm_profile $nuc_alignment`; 
#Run on each genome in directory
foreach my $genome(@genomes){
    my @genome_files = ();
    my @genome_IDs = ();
    unless ($genome =~ m/$cds_suffix/i){
	print $genome."\n";
	if ($genome =~ m/([\S]+).*\_$genome_suffix/){
	    my $genome_ID = $1;
	    my $wd = getcwd;
	    my $outdir = $wd."/".$genome_ID."_outfiles";
	    `mkdir $outdir`;
	    my $subdir = $outdir."/hmmer_files";
	    my $subdir2 = $outdir."/isoform_files";
	    my $nhmmer_nucleotide_sequences = $genome_ID.$nhmmer_unnanotated_seqfile;
	    my $nhmmer_file = $genome_ID."_assembly".$nhmmer_out; #nhmmer out file
	    my $cds_nuc = $genome_ID.$cds_nucleotide_seqfile;
	    my $cds_prot = $genome_ID.$cds_protein_seqfile;
	    my $cds_final_nuc = $genome_ID.$final_cds_nuc;
	    my $cds_final_prot = $genome_ID.$final_cds_prot;
	    my $phmmer_prot_isoforms = $genome_ID."_transcripts_all_isoforms_protein";
	    my $unique_longest_transcripts_out = $genome_ID.$nucleotide_longest_transcripts;
	    my $transcript_nucleotide_isoforms = $genome_ID."_transcripts_all_isoforms_nucleotide";
	    my $transcript_db = $genome_ID."_transcript_database";
	    my $transcript_blastout = $genome_ID."_tblastn.out";
	    my $nblast_db = $genome_ID."_nblast_db";
	    my $nblastout = $genome_ID."_nblast.out";
	    my $phmmer_file = $genome_ID.$phmmer_out; #phmmer out file
	    my $nhmmer_transcript_file = $genome_ID."_transcripts".$nhmmer_out;
	    my $nblast_cds_db = $genome_ID."_nblast_cds_db";
	    my $nblast_cds = $genome_ID."_nblast_cds.out";
	    my $nhmmer_assembly_log = $genome_ID."_nhmmer_assembly_domain_log.txt";
	    my $nhmmer_mads_domain_seqs = $genome_ID."_Domain_sequences_nuc.fa";
	    my @hit_log = ();
	    if ($annotation_available eq "yes" || $annotation_available eq "Yes"){
		my $nucleotide_ID = "";
		my $protein_ID = "";
		my $cds = "";
		foreach my $nucleotide(@nucleotide_transcripts){
		    if ($nucleotide =~ m/([\S]+).*\_$nt_transcript_suffix/){
			$nucleotide_ID = $1;
			if ($genome_ID eq $nucleotide_ID){
			    #print $nucleotide."\n";
			    push (@genome_files, $nucleotide);
			}
		    }
		}
		foreach my $protein(@protein_transcripts){
		    if ($protein =~ m/([\S]+).*\_$prot_transcript_suffix/){
			$protein_ID = $1;
			if ($genome_ID eq $protein_ID){ 
			   # print $protein."\n";
			    push(@genome_files, $protein);
			}
		    }
		}
		if ($cds_available eq "yes" || $cds_available eq "Yes"){
		    foreach my $cds(@cds_transcripts){
			if ($cds =~ m/([\S]+).*\_$cds_suffix/){
			    my $cds_ID = $1;
			    if ($genome_ID eq $cds_ID){
				#print $cds."\n";
				push (@genome_files, $cds);
			    }
			}
		    }
		}
		#phmmer:
		my $nucleotide = $genome_files[0];
		my $protein = $genome_files[1];
		print $nucleotide."\n";
		print $protein."\n";
		if ($cds_available eq "yes" || $cds_available eq "Yes"){
		    $cds = $genome_files[2];
		    print $cds."\n";   
		}
		if($default_phmmer_evalue eq "yes"){
		    `hmmsearch $phmm_profile $protein >> $phmmer_file`;
		}else{
		    `hmmsearch --incE $phmmer_evalue $phmm_profile $protein >> $phmmer_file`;
		}
		print "\nrunning phmmer on protein annotations...\n\n";
		my $phmmer_results;
		open(PHMMER, $phmmer_file);
		{
		    local $/; #set delimiter to nothing, enables phmmer file to be read in as one chunk
		    $phmmer_results = <PHMMER>; #store phmmer results
		}
		close PHMMER;
		`esl-sfetch --index $protein`; #index protein transcripts
		my $prot_ID = "";
		my @phmm_array = split(/\>\>/, $phmmer_results);
		my $phmmer_hit_chunk = $phmm_array[0];
		my @phmmer_array2 = split("Description", $phmmer_hit_chunk);
		my $phmmer_hit_chunk2 =  $phmmer_array2[1];
		my @phmmer_hits = ();
		if ($phmmer_hit_chunk2 =~ m/.*inclusion[\s]threshold.*/){
		    my @phmmer_array3 = split("------ inclusion threshold ------", $phmmer_hit_chunk2);
		    my $sig_phmmer_hits = $phmmer_array3[0];
		    @phmmer_hits = split("\n", $sig_phmmer_hits);
		    shift(@phmmer_hits); #remove rubbish element 1
		    shift(@phmmer_hits); #remove rubbish element 2
		    pop(@phmmer_hits); #remove empty line at end
		}
		else{
		    my @phmmer_array3 = split("\n\nDomain", $phmmer_hit_chunk2);
		    my $sig_phmmer_hits = $phmmer_array3[0];
		    @phmmer_hits = split("\n", $sig_phmmer_hits);
		    shift(@phmmer_hits); #remove rubbish element 1
		    shift(@phmmer_hits); #remove rubbish element 2
		}
		foreach my $phit(@phmmer_hits){
		    $phit =~ s/[\s]+/\|/g;
		    my @phmmdetails = split(/\|/, $phit);
		    shift @phmmdetails;
		    my $prot_ID = $phmmdetails[8];
		    #print $prot_ID."\n";
		    my $cmd = "esl-sfetch $protein $prot_ID >> $phmmer_prot_isoforms"; #fetch protein sequence
		    #print $cmd."\n";
		    `esl-sfetch $protein "$prot_ID" >> $phmmer_prot_isoforms`; #execute command
		}				
		my $tblastn_results;
		`makeblastdb -in $nucleotide -dbtype="nucl" -out $transcript_db`;
		`tblastn -db $transcript_db -query $phmmer_prot_isoforms -out $transcript_blastout`;
		open(TBLASTN, $transcript_blastout);
		{
		    local $/; #  change line delimter to nothing - read in file as one chunk
		    $tblastn_results = (<TBLASTN>);
		}
		close TBLASTN;
		`esl-sfetch --index $nucleotide`; #index nucleotide transcripts
		my @tblastn_hits = split(/Query\=/, $tblastn_results);
		shift @tblastn_hits;
		my @transcripts = (); #new
		my @locs = ();
		foreach my $query(@tblastn_hits){
		    my @queryhits = split(/\>/, $query);
		    my $continue = 1;
		    for(my $i=0;$i<scalar(@queryhits);$i++){
			my $tophit =">".$queryhits[$i]; #not element zero, this is the string before first hit
			if ($tophit =~ m/Identities[\s]+\=[\s]+([0-9]+)\/([0-9]+)[\s]+\(([0-9]+).*\)/){
			    my $identity = $3;
			    if ($identity eq 100){
				if ($tophit =~ m/(\>[\S]+)?\s.*/){
				    my $trans_id = $1;
				    $trans_id =~ s/\>//g;
				    my $loc = "";
				    if ($tophit =~ m/\(([\S]+)\)/){
					$loc = $1;
				    }
				    if ($trans_id ~~ @transcripts){
				    }
				    else{	
					#print $trans_id."\n";
					push @transcripts, $trans_id; #new
					unless($loc ~~ @locs){
					    push @locs, $loc;
					}
					my $cmd = "esl-sfetch $nucleotide $trans_id >> $transcript_nucleotide_isoforms";
					#print $cmd."\n";
					`esl-sfetch $nucleotide $trans_id >> $transcript_nucleotide_isoforms`;
				    }
				}
			    }
			}
		    }
		}
		## mine nhmmer on transcripts here, then do longest transcript per loc
		#print "\n nblast on transcripts ...\n";
		if($default_nhmmer_evalue eq "yes"){
		    `nhmmer $nhmm_profile $nucleotide >> $nhmmer_transcript_file`;
		}else{
		    `nhmmer --incE $nhmmer_evalue $nhmm_profile $nucleotide >> $nhmmer_transcript_file`;
		}
		print "running nhmmer on mRNA transcripts ...\n\n";
		my $nhmmer_t_results = "";
		open(NHMMER_T, $nhmmer_transcript_file);
		{
		    local$/;
		    $nhmmer_t_results = <NHMMER_T>;
		}
		close NHMMER_T;
		my @nhmm_t_array = split(/\>\>/, $nhmmer_t_results);
		my $nhmmer_t_hit_chunk = $nhmm_t_array[0];
		my @nhmm_t_array2 = split("Description", $nhmmer_t_hit_chunk);
		my $nhmmer_t_hit_chunk2 = $nhmm_t_array2[1];
		my @nhmmer_t_hits = ();
		if ($nhmmer_t_hit_chunk2 =~ m/.*inclusion[\s]threshold.*/){
		    my @nhmmer_t_array3 = split("------ inclusion threshold ------", $nhmmer_t_hit_chunk2);
		    my $sig_nhmmer_t_hits = $nhmmer_t_array3[0];
		    @nhmmer_t_hits = split("\n", $sig_nhmmer_t_hits);
		    shift(@nhmmer_t_hits);
		    shift(@nhmmer_t_hits);
		    pop(@nhmmer_t_hits);
		}
		else{
		    my @nhmmer_t_array3 = split("\n\n", $nhmmer_t_hit_chunk2);
		    my $sig_nhmmer_t_hits = $nhmmer_t_array3[0];
		    @nhmmer_t_hits = split("\n", $sig_nhmmer_t_hits);
		    shift(@nhmmer_t_hits);
		    shift(@nhmmer_t_hits);
		}
		foreach my $trans_hit(@nhmmer_t_hits){
		    $trans_hit =~ s/[\s]+/\|/g;
		    my @nhmmtdetails = split(/\|/, $trans_hit);
		    shift(@nhmmtdetails);
		    my $transcript_ID = $nhmmtdetails[3];
		    unless($transcript_ID ~~ @transcripts){
			#print "new transcript: $transcript_ID \n";
			my $cmd = "esl-sfetch $nucleotide $transcript_ID \>\> $transcript_nucleotide_isoforms";
			#print $cmd."\n";
			`esl-sfetch $nucleotide $transcript_ID >> $transcript_nucleotide_isoforms`;
		   }
		}
		#################################################################################################
		#pull only unique locs, longest transcript
		my $transcript_file;
		open(TRANSCRIPTS, $transcript_nucleotide_isoforms);
		{
		    local $/;
		    $transcript_file = <TRANSCRIPTS>;
		}
		close TRANSCRIPTS;
		my @transcript_seqs = split(/\>/, $transcript_file);
		foreach my $LOC(@locs){
		    my @uniquelocs = ();
		    #print $LOC."-----------------------\n";
		    my %transcript_lengths;
		    foreach my $tseq(@transcript_seqs){
			if($tseq =~ m/(.*)\n([\S\n]+)/){
			    my $theader = ">".$1;
			    my $tseq = $2;
			    my $tID ="";
			    my $locID ="";
			    my $tlength = length($tseq);
			    if ($theader =~ m/\>([\S]+)?\s.*/){
				$tID = $1;
			    }
			    if ($theader =~ m/\(([\S]+)\)/){
				$locID = $1;
			    }
			    if($locID eq $LOC){
				$transcript_lengths{$tID} = $tlength; #store header as key and evalue as value in hash as pair
			    }	
			}
		    }
		    foreach my $traID(sort { $transcript_lengths{$a} <=> $transcript_lengths{$b} or $a cmp $b } keys %transcript_lengths){
			#print $transcript_lengths{$traID}."\n";
			push(@uniquelocs, $traID);
		    }
		    my $longest_transcript = pop(@uniquelocs);
		    `esl-sfetch $nucleotide $longest_transcript >> $unique_longest_transcripts_out`;
		}
		#################################################################################################
		#blast nucleotide transcripts to assembly
		`makeblastdb -in $genome -dbtype="nucl" -out $nblast_db`; #make blast database using query genome
		`blastn -db $nblast_db -query $unique_longest_transcripts_out -out $nblastout`; #blast assembly with nhmmer results
		#get coordinates
		my $nblast_results;
		my @block_coordinates = ();
		my $contig;
		my $start;
		my $end;
		my $coordinates;
		open(BLASTN, $nblastout);
		{
		    local $/;
		    $nblast_results = <BLASTN>;
		}
		close BLASTN;
		my @blastn_hits = split(/Query\=/, $nblast_results);
		shift @blastn_hits;
		foreach my $query(@blastn_hits){
		    my @start_coords = ();
		    my @end_coords = ();
		    my $tmp;
		    my @queryhits = split(/\>/, $query);
		    my $tophit =">".$queryhits[1]; #not element zero, this is the string before first
		    my @exon_hits = split(/Score/, $tophit);
		    my $contig_details = $exon_hits[0];
		    if ($contig_details =~ m/\>([\S]+)?\s.*/){
			$contig = $1;
		    }
		    shift(@exon_hits);
		    foreach my $exon(@exon_hits){
			if ($exon =~ m/Identities[\s]+\=[\s]+([0-9]+)\/([0-9]+)[\s]+\(([0-9]+).*\)/){
			    #print "this is an exon:\n$exon\n"; 
			    if($exon=~m/Sbjct[\s]+([0-9]+)/){
				$start = $1;
			    }
			    while($exon=~s/Sbjct[\s]+[0-9]+[\s]+[\S]+[\s]+([0-9]+)//){
				$end = $1;
			    }
			    if($start > $end){
				$tmp = $start;
				$start = $end;
				$end = $tmp;
			    }
			    push(@start_coords, $start);
			    push(@end_coords, $end);
			}#we want gene start and end, not exon start and end - so we block out intronic regions too
		    }
		    my @sorted_start = sort { $a <=> $b } @start_coords;
		    my @sorted_end = sort { $a <=> $b } @end_coords;
		    my $true_start = shift(@sorted_start); #gene start 
		    my $true_end = pop(@sorted_end); #gene end
		    $coordinates = $contig."|".$true_start."|".$true_end;
		    #print $coordinates."\n";
		    push @block_coordinates, $coordinates;
		}
		#####################################################################################################
		#mine cds and output
		#mine corresponding proteins (longest proteins) and output
		if ($cds_available eq "yes" || $cds_available eq "Yes"){
		    `makeblastdb -in $cds -dbtype="nucl" -out $nblast_cds_db`;
		    `blastn -db $nblast_cds_db -query $unique_longest_transcripts_out -out $nblast_cds`;
		    `esl-sfetch --index $cds`;
		    my $cds_blasthits = "";
		    my @identifiers_cds = ();
		    open(CDSBLAST, $nblast_cds);
		    {
			local $/;
			$cds_blasthits = <CDSBLAST>;
		    }
		    close CDSBLAST;
		    my @cds_hits = split(/Query\=/, $cds_blasthits);
		    shift @cds_hits;
		    foreach my $cds_hit(@cds_hits){
			#my @protein_cds_info = ();
			my @cdshits = split(/\>/, $cds_hit);
			my $cds_1 = $cdshits[1];
			my $top_cds = ">".$cds_1;
			if ($top_cds =~ m/Identities[\s]+\=[\s]+([0-9]+)\/([0-9]+)[\s]+\(([0-9]+).*\)/){
			    my $identity = $3;
			    if ($identity eq 100){
				if ($top_cds =~ m/(\>[\S]+)?\s.*/){
				    my $identifier = $1;
				    $identifier =~ s/\>//g;
				    unless($identifier ~~ @identifiers_cds){
					#print $identifier."\n";
					`esl-sfetch $cds \"$identifier\" >> $cds_nuc`;
					push(@identifiers_cds, $identifier);
				    }
				}
			    }
			}
		    }
		    my @identifiers_protein_cds = ();
		    open(CDS,  $cds_nuc);
		    {
			local $/ = ">"; #read in by seq
			while(<CDS>){
			    my $seq = $_;
			    if ($seq =~ m/.*\[protein\_id\=([\S]+)?\].*/){
				my $protein_cds = $1;
				unless($protein_cds ~~ @identifiers_protein_cds){
				    #print "get $protein_cds !\n";
				    `esl-sfetch $protein $protein_cds >> $cds_prot`;
				    push(@identifiers_protein_cds, $protein_cds);
				}
			    }
			}
		    }
		}
		close CDS;
		#####################################################################################################
		#nhmmer on assembly to discover new hits
		#nhmmer then discount anything which overlaps with existing coordinates
		print "running nhmmer on whole assembly to pull new hits ...\n\n";
		`esl-sfetch --index $genome`; #index genomefile
		my @nhmmer_coordinates = ();
		if($default_nhmmer_evalue eq "yes"){
		    `nhmmer $nhmm_profile $genome >> $nhmmer_file`;
		}else{
		    `nhmmer --incE $nhmmer_evalue $nhmm_profile $genome >> $nhmmer_file`;
		}
		my $nhmmer_results;
		open(NHMMER, $nhmmer_file);
		{
		    local $/;
		    $nhmmer_results = <NHMMER>;
		}
		close NHMMER;
		my @nhmm_array = split(/\>\>/, $nhmmer_results);
		my $nhmmer_hit_chunk = $nhmm_array[0];
		my @nhmmer_array2 = split("Description", $nhmmer_hit_chunk);
		my $nhmmer_hit_chunk2 =  $nhmmer_array2[1];
		my @nhmmer_hits = ();
		if ($nhmmer_hit_chunk2 =~ m/.*inclusion[\s]threshold.*/){
		    my @nhmmer_array3 = split("------ inclusion threshold ------", $nhmmer_hit_chunk2);
		    my $sig_nhmmer_hits = $nhmmer_array3[0];
		    @nhmmer_hits = split("\n", $sig_nhmmer_hits);
		    shift(@nhmmer_hits); #remove rubbish element
		    shift(@nhmmer_hits); #remove rubbish element
		    pop(@nhmmer_hits); #remove empty line at end
		}else{
		    my @nhmmer_array3 = split("\n\n", $nhmmer_hit_chunk2);
		    my $sig_nhmmer_hits = $nhmmer_array3[0];
		    #print $sig_nhmmer_hits."\n";
		    @nhmmer_hits = split("\n", $sig_nhmmer_hits);
		    shift(@nhmmer_hits); #remove rubbish element
		    shift(@nhmmer_hits); #remove rubbish element
		}
		foreach my $nhit(@nhmmer_hits){
		    my $ntmp;
		    $nhit =~ s/[\s]+/\|/g;
		    my @nhmmdetails = split(/\|/, $nhit);
		    shift @nhmmdetails;
		    my $ncontig = $nhmmdetails[3]; #print $ncontig."\n";
		    my $nstart =  $nhmmdetails[4]; #print $nstart."\n";
		    my $nend =  $nhmmdetails[5]; #print $nend."\n";
		    my $nhmm_hit = $ncontig."|".$nstart."|".$nend;
		    push @nhmmer_coordinates, $nhmm_hit;
		}
		open(OUT, ">>$nhmmer_assembly_log");
		foreach my $nhmm_locus(@nhmmer_coordinates){
		    my $overlap = "";
		    my $ntmp = "";
		    my $strand = "";
		    my @details = split(/\|/, $nhmm_locus);
		    my $contig_n = $details[0]; #need to add Xnts 
		    my $start_n = $details[1]; #need to add Xnts
		    my $end_n = $details[2];
		    if($start_n > $end_n){
			#print "reversed\n";
			$ntmp = $start_n;
			$start_n = $end_n;
			$end_n = $ntmp;
			$strand = "rev";
		    }
		    else{
			$strand = "pos";
		    }
		    my $nhit_length = $end_n - $start_n;
		    foreach my $stored_hits(@block_coordinates){
			my $min_coord = "";
			my $max_coord = "";
			my $max_span = "";
			my $total_length = "";
			my @stored_details = split(/\|/, $stored_hits);
			my $contig_stored = $stored_details[0];
			my $start_stored = $stored_details[1];
			my $end_stored = $stored_details[2];
			my $stored_hit_length = $end_stored - $start_stored;
			if ($contig_n eq $contig_stored){
			    if ($start_n < $start_stored){
				$min_coord = $start_n;
			    }
			    else{
				$min_coord = $start_stored;
			    }
			    if ($end_n > $end_stored){
				$max_coord = $end_n;
			    }
			    else{
				$max_coord = $end_stored;
			    }
			    $max_span = $max_coord - $min_coord;
			    $total_length = $nhit_length + $stored_hit_length;
			    if($max_span<$total_length){
				$overlap = "Y";
			    }else{ #push coordinates to stored and use esl-sfetch to pull hit +- X nts
			    }
			}
		    }
		    if($overlap ne "Y"){
			my $contig_out = "contig.fa";
			my $contig_seq = "";
			my $contig_length = "";
			my $cmd_1 = "esl-sfetch $genome $contig_n >> $contig_out";
			#print $cmd_1."\n";
			`$cmd_1`;
			open(IN, $contig_out);
			{
			    local $/;
			    $contig_seq = <IN>;
			}
			close IN;
			`rm $contig_out`;
			if($contig_seq =~ m/(.*)\n([\S\n]+)/){
			    my $chead = $1;
			    my $cseq = $2;
			    #print $chead."\n";
			    $cseq =~ s/\n//g;
			    $contig_length = length($cseq);
			    #print $contig_length."\n";
			}
			my $nhmm_dets = $nhmm_locus."|".$contig_length;
			print OUT $nhmm_dets."\n";
			#print $nhmm_dets."\n";
			#print "new hit: ";
			my $dom_range = $start_n."\.\.".$end_n;
			if ($strand eq "rev"){
			    $start_n -= $nhmmer_plus; #3' end is $start (hence - plus)
			    if($start_n < 1){
				$start_n = 1;
			    }
			    $end_n += $nhmmer_minus; #5' end is $end (hence + minus)
			    if($end_n > $contig_length){
				$end_n = $contig_length;
			    }
			    # print "reverse strand\n";
			    my $range = $start_n."\.\.".$end_n;
			    my $cmd = "esl-sfetch -c $range -r $genome $contig_n >> $nhmmer_nucleotide_sequences";
			    #print $cmd."\n";
			    `$cmd`;
			    my $cmd2 = "esl-sfetch -c $dom_range -r $genome $contig_n >> $nhmmer_mads_domain_seqs ";
			    `$cmd2`;
			}
			elsif($strand eq "pos"){
			    $start_n -= $nhmmer_minus;
			    if($start_n < 1){
				$start_n = 1;
			    }
			    $end_n += $nhmmer_plus;
			    if($end_n > $contig_length){
				$end_n = $contig_length;
			    }
			    #print "positive strand\n";
			    my $range = $start_n."\.\.".$end_n;
			    my $cmd = "esl-sfetch -c $range $genome $contig_n >> $nhmmer_nucleotide_sequences";
			    #print $cmd."\n";
			    `$cmd`;
			    my $cmd2 = "esl-sfetch -c $dom_range $genome $contig_n >> $nhmmer_mads_domain_seqs";
			    `$cmd2`;
			}
			#print $contig_n."...". $start_n."...".$end_n."\n";
		    }
		}
		`mkdir $subdir`;
		`mkdir $subdir2`;
		`mv $nhmmer_file $phmmer_file $nhmmer_transcript_file $subdir`; #$nhmm_profile #phmm_profile add this to remove, I just don't have alignment 
		`mv $phmmer_prot_isoforms $transcript_nucleotide_isoforms $subdir2`;
		if ($cds_available eq "yes" || $cds_available eq "Yes"){
		    `cp $cds_nuc $cds_final_nuc`;
		    `cp $cds_prot $cds_final_prot`;
		}
		`rm $transcript_db* $nblast_db* $nblast_cds_db* *nblast.out *tblastn.out $nblast_cds`;
		`rm *ssi`;
	    }
	    else{ #only run nhmmer on assembly
		print "running nhmmer on whole assembly to pull hits ...\n\n";
		`esl-sfetch --index $genome`;
		my @nhmmer_coordinates = ();
		if($default_nhmmer_evalue eq "yes"){
		    `nhmmer $nhmm_profile $genome >> $nhmmer_file`;
		}else{
		    `nhmmer --incE $nhmmer_evalue $nhmm_profile $genome >> $nhmmer_file`;
		}
		my $nhmmer_results;
		open(NHMMER, $nhmmer_file);
		{
		    local $/;
		    $nhmmer_results = <NHMMER>;
		}
		close NHMMER;
		my @nhmm_array = split(/\>\>/, $nhmmer_results);
		my $nhmmer_hit_chunk = $nhmm_array[0];
		my @nhmmer_array2 = split("Description", $nhmmer_hit_chunk);
		my $nhmmer_hit_chunk2 =  $nhmmer_array2[1];
		my @nhmmer_hits = ();
		if ($nhmmer_hit_chunk2 =~ m/.*inclusion[\s]threshold.*/){
		    my @nhmmer_array3 = split("------ inclusion threshold ------", $nhmmer_hit_chunk2);
		    my $sig_nhmmer_hits = $nhmmer_array3[0];
		    @nhmmer_hits = split("\n", $sig_nhmmer_hits);
		    shift(@nhmmer_hits); #remove rubbish element
		    shift(@nhmmer_hits); #remove rubbish element
		    pop(@nhmmer_hits); #remove empty line at end
		}else{
		    my @nhmmer_array3 = split("\n\n", $nhmmer_hit_chunk2);
		    my $sig_nhmmer_hits = $nhmmer_array3[0];
		    @nhmmer_hits = split("\n", $sig_nhmmer_hits);
		    shift(@nhmmer_hits); #remove rubbish element
		    shift(@nhmmer_hits); #remove rubbish element
		}
		foreach my $nhit(@nhmmer_hits){
		    my $ntmp;
		    $nhit =~ s/[\s]+/\|/g;
		    my @nhmmdetails = split(/\|/, $nhit);
		    shift @nhmmdetails;
		    my $ncontig = $nhmmdetails[3]; #print $ncontig."\n";
		    my $nstart =  $nhmmdetails[4]; #print $nstart."\n";
		    my $nend =  $nhmmdetails[5]; #print $nend."\n";
		    my $nhmm_hit = $ncontig."|".$nstart."|".$nend;
		    push @nhmmer_coordinates, $nhmm_hit;
		}
		open(OUT, ">>$nhmmer_assembly_log");
		foreach my $nhmm_locus(@nhmmer_coordinates){
		    #print OUT $nhmm_locus."\n";
		    my $ntmp = "";
		    my $strand = "";
		    my @details = split(/\|/, $nhmm_locus);
		    my $contig_n = $details[0];
		    my $start_n = $details[1];
		    my $end_n = $details[2];
		    if($start_n > $end_n){
			$ntmp = $start_n;
			$start_n = $end_n;
			$end_n = $ntmp;
			$strand = "rev";
		    }
		    else{
			$strand = "pos";
		    }
		    #ensure that limit is not surpassed:
		    my $contig_out = "contig.fa";
		    my $contig_seq = "";
		    my $contig_length = "";
		    my $cmd_1 = "esl-sfetch $genome $contig_n >> $contig_out";
		    #print $cmd_1."\n";
		    `$cmd_1`;
		    open(IN, $contig_out);
		    {
			local $/;
			$contig_seq = <IN>;
		    }
		    close IN;
		    `rm $contig_out`;
		    if($contig_seq =~ m/(.*)\n([\S\n]+)/){
			my $chead = $1;
			my $cseq = $2;
			#print $chead."\n";
			$cseq =~ s/\n//g;
			$contig_length = length($cseq);
			#print $contig_length."\n";
		    }
		    my $nhmm_dets = $nhmm_locus."|".$contig_length;
		    print OUT $nhmm_dets."\n";
		    if($strand eq "rev"){
			my $range_domain = $start_n."\.\.".$end_n;
			my $cmd1 = "esl-sfetch -c $range_domain -r $genome $contig_n >> $nhmmer_mads_domain_seqs";
			`$cmd1`;
			$start_n -= $nhmmer_plus; #3' end is $start (hence - plus)
			$end_n += $nhmmer_minus; #5' end is $end (hence + minus)
			if ($start_n < 1){
			    $start_n = 1;
			}
			if ($end_n > $contig_length){
			    #print "exceeds contig range, reverting to max contig bounds!\n";
			    $end_n = $contig_length;
			}
			#print "reverse strand\n";
			my $range = $start_n."\.\.".$end_n;
			my $cmd = "esl-sfetch -c $range -r $genome $contig_n >> $nhmmer_nucleotide_sequences";
			#print $cmd."\n";
			`$cmd`;
		    }
		    elsif($strand eq "pos"){
			my $range_domain = $start_n."\.\.".$end_n;
			my $cmd1 = "esl-sfetch -c $range_domain $genome $contig_n >> $nhmmer_mads_domain_seqs";
			`$cmd1`;
			$start_n -= $nhmmer_minus;
			$end_n += $nhmmer_plus;
			if ($start_n < 1){
			    $start_n = 1;
			}
			if ($end_n > $contig_length){
			   # print "exceeds contig range, reverting to max contig bounds!\n";
			    $end_n = $contig_length;
			}
			#print "positive strand\n";
			my $range = $start_n."\.\.".$end_n;
			my $cmd = "esl-sfetch -c $range $genome $contig_n >> $nhmmer_nucleotide_sequences";
			#print $cmd."\n";
			`$cmd`;
		    }
		    #print $contig_n."...". $start_n."...".$end_n."\n";
		}
		close OUT;
		my $subdir = $outdir."/hmmer_files";
		`mkdir $subdir`;
		`mv $nhmmer_file $subdir`;
		`rm *ssi`;
	    }
	    if ($predict_new_hits eq "Yes" || $predict_new_hits eq "yes"){
		print "running augustus to predict new unannotated hits ...\n\n";
		`esl-sfetch --index $reference_file`;
		my @seqs = ();
		my $sequences = "";
		my $seq_db_pre = $genome_ID."_";
		my $hit_no = 0;
		open(SEQS, $nhmmer_nucleotide_sequences);
		{
		    local $/;
		    $sequences = <SEQS>;
		}
		close SEQS;
		@seqs = split (/\>/, $sequences);
		open(IN,  $nhmmer_assembly_log);
		my $domain_det_list = "";
		{
		    local $/;
		    $domain_det_list = <IN>;
		}
		close IN;
		my @domain_details = split(/\n/, $domain_det_list);
		#foreach my $detail(@domain_details){
		#    $detail =~ s/\n//;
		#    print $detail."\n";
		#}
		foreach my $seq(@seqs){
		    if($seq =~ m/(.*)\n([\S\n]+)/){
			my $seq_header = $1;
			my $seq_nt = $2;
			my $domain_info = $domain_details[0];
			my @info = split(/\|/, $domain_info);
			my $contig_name = $info[0];
			my $domain_info_start = $info[1];
			my $domain_info_end = $info[2];
			my $contig_max = $info[3];
			if($domain_info_start > $domain_info_end){
			    my $bu = $domain_info_start;
			    $domain_info_start = $domain_info_end;
			    $domain_info_end = $bu;
			}
			#print $domain_info."\n";
			#print "This is domain start $domain_info_start \n";
			#print "This is domain end $domain_info_end \n";
			#print "This is max contig $contig_max \n";
			shift(@domain_details);
			$hit_no +=1;
			my $hit_annotation = "Hit".$hit_no;
			push(@hit_log, $hit_annotation);
			my $contig = "";
			my $start = "";
			my $end = "";
			my $domain_scoord = "";
			my $domain_ecoord = "";
			my $domain_length = "";
			my $domain_start = "";
			my $domain_end = "";
			if ($seq_header =~ m/^([\S]+)?\/([0-9]+)\-([0-9]+)\s.*/){
			    $contig = $1; $start = $2; $end = $3;
			}
			my $a = $domain_info_start - $nhmmer_minus; #position a 
			if($a < 1){
			    my $x = 1 - $a;
			    my $y = $nhmmer_minus - $x;
			    $domain_scoord = 1 + $y;
			}
			elsif($a >= 1){
			    $domain_scoord = 1 + $nhmmer_minus;
			}
			$domain_length = ($domain_info_end - $domain_info_start) +1;
			$domain_ecoord = $domain_scoord + $domain_length;
			$domain_start = $domain_scoord;
			$domain_end = $domain_ecoord;
			$seq_header =~ s/\//\_/g;  #NC_044377.1/77156829-77158035 Cannabis sativa chromosome 6,
			$seq_header =~ s/\-/\_/g;
			$seq_header = ">Hit".$hit_no."_".$seq_header;
			my $tmp_out = "tmp.fa";
			open(OUT, ">$tmp_out");
			print OUT $seq_header."\n".$seq_nt;
			close OUT;
			my $reference = "";
			my $psl = "Hit".$hit_no."_ref.psl";
			my $hints = "Hit".$hit_no."_hints.gff";
			my $reference_outseq = "Hit".$hit_no."_ref.fa";
			my $prediction_gff = "Hit".$hit_no."prediction_out.gff";
			`blat -minIdentity=$minidentity $tmp_out $reference_file $psl`;
			`perl blat2hints.pl --in=$psl --out=$hints`;
			system("augustus --species=$augustus_species --strand=forward --codingseq=on --softmasking=0 --hintsfile=$hints --extrinsicCfgFile=extrinsic.ME.cfg $tmp_out > $prediction_gff");
			my $prediction_in = "";
			open(AUG, $prediction_gff);
			{
			    local $/;
			    $prediction_in = <AUG>;
			}
			close AUG;
			my @intron_coords = ();
			my @exon_coords = ();
			my @predictions = split(/#[\s]start[\s]gene/, $prediction_in);
			shift @predictions;
			foreach my $pred(@predictions){
			    #`echo \"-----------------------------\n\" >> predictions_log_ALL.gff`;
			    #`echo \"$hit_details\" >> predictions_log_ALL.gff`;
			    #`echo \"$pred\" >> predictions_log_ALL.gff`;
			    #print "prediction: $pred \n";
			    my @pred_split = split (/#[\s]protein[\s]sequence/, $pred);
			    my $gene_cds_details = $pred_split[0];
			    my $protein_details = $pred_split[1];
			    my @split_again = split(/#[\s]coding[\s]sequence/, $gene_cds_details);
			    my $gene_details = $split_again[0];
			    my $cds_details = $split_again[1];
			    my @gene_dets = split(/\n/, $gene_details);
			    shift(@gene_dets);
			    my $pred_start = "";
			    my $pred_end = "";
			    my $maximum_end = "";
			    my $minimum_start = "";
			    my $total_length3 = "";
			    my $max_span3 = "";
			    my $prediction_length = "";
			    my $prediction_found = 0;
			    my @prediction_coords = ();
			    foreach my $gd(@gene_dets){
				$gd =~ s/[\s]+/\|/g;
				if ($gd =~ m/\|gene\|([0-9]+)\|([0-9]+)/){
				    $pred_start = $1; $pred_end = $2;
				    push(@prediction_coords, $pred_start);
				    push(@prediction_coords, $pred_end);
				}
			    }
			    $pred_start = $prediction_coords[0];
			    $pred_end = $prediction_coords[1];
			    #print "Prediction start: $pred_start \n";
			    #print "Prediction end: $pred_end \n";
			    if($pred_end > $domain_end){
				$maximum_end = $pred_end;
			#	print "Maximum end is pred end: $pred_end \n";
			    }else{
				$maximum_end = $domain_end;
			#	print "Maximum end is domain end: $domain_end \n";
			    }
			    if($pred_start < $domain_start){
				$minimum_start = $pred_start;
			#	print "Minimum start is pred start: $pred_start \n";
			    }else{
				$minimum_start = $domain_start;
			#	print "Minimum start is domain start: $domain_start \n";
			    }
			    $prediction_length = ($pred_end - $pred_start) +1;
			    $total_length3 = $prediction_length + $domain_length;
			    $max_span3 = ($maximum_end - $minimum_start) +1;
			    if ($max_span3 < $total_length3){
			#	print "we want this prediction!\n";
				#print $gd."\n";
				$prediction_found = 1;
				my $hit_details = "Hit".$hit_no." prediction:\n";
				`echo \"-----------------------------\n\" >> predictions_log.gff`;
				`echo \"$hit_details\" >> predictions_log.gff`;
				`echo \"$pred\" >> predictions_log.gff`;
			    }else{
				$prediction_found = 0;
			    }
			    $protein_details =~ s/\n//g;
			    $cds_details =~ s/\n//g;
			    if ($prediction_found  == 1){
				my $cds_header = ">Hit".$hit_no."_new_augustus_prediction\n";
				my $hit_annotation = "Hit".$hit_no;
				open(CDS_NT, ">>$cds_final_nuc");
				if ($cds_details =~ m/\[([^\]]+)\]/){
				    my $cds_seq = $1;
				    $cds_seq =~ s/\#//g;
				    $cds_seq =~ s/\s//g;
				    $cds_seq =~ s/.{80}\K/\n/g;
				    $cds_seq = uc($cds_seq);
				    #print $cds_header;
				    #print "cds:\n$cds_seq\n";
				    print CDS_NT $cds_header.$cds_seq."\n";
				}
				close CDS_NT;
				open(CDS_PROT, ">>$cds_final_prot");
				if ($protein_details =~ m/\[([^\]]+)\]/){
				    my $protein_seq = $1;
				    $protein_seq =~ s/\#//g;
				    $protein_seq =~ s/\s//g;
				    $protein_seq =~ s/.{80}\K/\n/g;
				   # print "protein:\n$protein_seq\n";
				   print CDS_PROT $cds_header.$protein_seq."\n";
				}
				close CDS_PROT;
			    }
			}
		    }
		}
		`rm Hit* *ssi tmp.fa`;
		my @domain_seqs = "";
		if(-e $nhmmer_mads_domain_seqs){
		    open(DOMAINS,$nhmmer_mads_domain_seqs);
		    {
			local $/ = ">";
			while(<DOMAINS>){
			    my $dom = $_;
			    if($dom =~ m/\>/){
				$dom =~ s/\>//g;
				$dom = ">".$dom;
			    }
			    else{
				$dom = ">".$dom;
			    }
			    push (@domain_seqs, $dom);
			    #my $dom = $_;
			    #push(@domain_seqs, $dom); 
			}
		    }
		    close DOMAINS;
		#foreach my $n(@hit_log){
		#    print "n: $n \n";
		#}
		    my $tmp_doms = "tmp_doms.fa";
		    open(DOMS, ">>$tmp_doms");
		    foreach my $domain(@domain_seqs){
			if($domain=~m/(\>[^\n]+)\n([\S\n]+)/){
			    my $dom_header = $1;
			    my $dom_seq = $2;
			    my $hit_label = $hit_log[0];
			    shift(@hit_log);
			    $dom_header =~ s/\>//g;
			    $dom_header = ">".$hit_label."_".$dom_header."\n";
			    print DOMS $dom_header.$dom_seq;
			    #print $domain."\n";
			}
		    }
		    `mv $tmp_doms $nhmmer_mads_domain_seqs`;
		}
	    }
	    if ($pseudogene_check eq "yes" || $pseudogene_check eq "Yes"){
		print "Annotating hits with functional or pseudogene status ...\n\n";
		open(CDS, $cds_final_nuc);
		my $out = "tmp_cds.fa";
		my $final_seqs = "";
		{
		    local $/;
		    $final_seqs = <CDS>;
		}
		close CDS;
		my @finalseqs = split (/\>/, $final_seqs);
		foreach my $cds(@finalseqs){
		    my $annotation = "";
		    my @status = ();
		    if($cds =~ m/(.*)\n([\S\n]+)/){
			my $cds_header = $1;
			my $cds_seq = $2;
			$cds_seq =~ s/\n//g; #remove new lines
			#print $cds_seq."\n";
			my $cds_length = length($cds_seq);
			#$cds_seq =~ s/.{80}\K/\n/g;
			@status=&checkframe($cds_seq);
			my $stat = $status[0];
			if($stat eq "1"){ 
			    $annotation=$annotation_1; #START codon && no in frame stop codons
			}
			if($stat eq "2"){
			    $annotation=$annotation_2; #no START codon && no stop codons in any frame   
			}
			if($stat eq "3"){
			    $annotation=$annotation_3; #no START codon && no in frame stop codons   
			}
			if($stat eq "4"){
			    $annotation=$annotation_4; #START codon && in frame stop codon  
			}
			if($stat eq "5"){
			    $annotation=$annotation_5; #no START codon && in frame stop codon 
			}
			if($stat eq "6"){
			    $annotation=$annotation_6; #no STARt codon and && stop codons in all frames
			}
			if ($cds_length < $pseudogene_length){
			    $annotation = $annotation_short;
			}
			$cds_header = ">".$cds_header."_".$annotation."\n";
			$cds_seq =~ s/.{80}\K/\n/g;
			open(OUT_FINAL, ">>$out");
			print OUT_FINAL $cds_header.$cds_seq."\n";
			close OUT_FINAL;
		    }
		}
		`mv $out $cds_final_nuc`;
	    }
	    if(-e $nhmmer_nucleotide_sequences){
		my $subdir4 = $outdir."/unnanotated_hits";
		`mkdir $subdir4`;
		`mv $nhmmer_nucleotide_sequences $subdir4`;
	    }
	    `mv $cds_final_nuc $cds_final_prot $outdir`;
	    if ($annotation_available eq "Yes" || $annotation_available eq "yes"){
		my $subdir3 = $outdir."/mined_annotation_files";
		`mkdir $subdir3`;
		`mv $cds_nuc $cds_prot $subdir3`;
		`mv $unique_longest_transcripts_out $subdir3`;
		`mv $subdir2 $subdir3`;
	    }
	    if($predict_new_hits eq "Yes" || $predict_new_hits eq "yes"){
		`mv predictions_log.gff $outdir`;
	    }
	    if(-e $nhmmer_nucleotide_sequences){
		`mv $nhmmer_nucleotide_sequences $outdir`;
	    }
	    if( -e $nhmmer_assembly_log){
		`rm $nhmmer_assembly_log `;
	    }
	    if(-e $nhmmer_mads_domain_seqs){
		`mv $nhmmer_mads_domain_seqs $outdir`;
	    }
	}
    }
}
if ($annotation_available eq "yes" || $annotation_available eq "Yes"){
    `for dir in *outfiles; do cp $phmm_profile \$dir/hmmer_files;done`;
    `for dir in *outfiles; do cp $nhmm_profile \$dir/hmmer_files;done`;
    `rm $phmm_profile`;
    `rm $nhmm_profile`;
}
else{
    `for dir in *outfiles; do cp $nhmm_profile \$dir/hmmer_files;done`;
    `rm $nhmm_profile`;
}
print "Complete!\n";

sub checkframe{
    my $status=0; #holds annotation status 
    my $dnaseqf1=uc($_[0]); #takes prediction as input (frame 1)
    my @framedata=();
    my $ignore_seq=0;
    if($dnaseqf1=~m/TGA$/){
	$dnaseqf1=~s/TGA$//; #remove stop codon
    }
    elsif($dnaseqf1=~m/TAA$/){
	$dnaseqf1=~s/TAA$//; #remove stop codon
    }
    elsif($dnaseqf1=~m/TAG$/){
	$dnaseqf1=~s/TAG$//; #remove stop codon
    }
    else{ #if no stop codon do nothing
    }
    my $dnaseqf2=substr($dnaseqf1,1,(length($dnaseqf1)-1)); #frame2
    my $dnaseqf3=substr($dnaseqf1,2,(length($dnaseqf1)-2)); #frame3
    my $f1stop=0; my $f2stop=0; my $f3stop=0; #declare variables to track stop codons in each frame
    my @framecheck=();
    $framecheck[0]=$dnaseqf1;$framecheck[1]=$dnaseqf2;$framecheck[2]=$dnaseqf3; #store counts in @framecheck array
    for(my $k=0;$k<scalar(@framecheck);$k++){
	for(my $n=0;$n<length($framecheck[$k])-2;$n+=3){ #read in codons for each frame
	    my $codon=substr($framecheck[$k],$n,3);
	    if($codon eq "TGA" || $codon eq "TAA" || $codon eq "TAG"){
		if($k==0){
		    $f1stop++; #stop codon in frame 1
		}
		elsif($k==1){
		    $f2stop++; #stop codon in frame 2
		}
		elsif($k==2){
		    $f3stop++; #stop codon in frame 3
		}
	    }
	}
    }
    if($dnaseqf1=~m/^ATG/){ #if start codon
	if($f1stop==0){
	    $status=1; ##1 means start codon, with no in frame stop codons
	}
	else{
	    $status=4; #4 means start codon, but contains at least 1 in frame stop codon
	}
    }
    else{ #no start codon
	if($f1stop==0 && $f2stop==0 && $f3stop==0){
	    $status=2; #2 means no start but no stop codons in any frame
	}
	elsif($f1stop>0 && $f2stop>0 && $f3stop>0){
	    $status=5; #5 means no start and stop codons in in all frames 
	}
	elsif($f1stop>0){ #at least frame1 has stop codon
	    $status =6; #6 means no START codon and in frame stop codon
	}
	else{ #only frames 2 and 3 have stop codon
	    $status=3; #3 means no START codon and no in frame stop codon (stop codon in frame 2 and 3) 
	}
    }
    push @framedata,$status; #push status to @framedata array
    return @framedata; #subroutine returns this array with annotation status
}
