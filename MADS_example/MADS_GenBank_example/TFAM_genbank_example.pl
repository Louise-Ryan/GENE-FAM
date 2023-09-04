#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Scalar::Util qw(looks_like_number);

###########################################################################
#USER PARAMETERS:                                                         #
###########################################################################
#Annotation files available?
my $annotation_available = "no"; #If NCBI annotations are available for your genome set below variable to "yes". Set as "no" if no annotations are available, and you wish to mine the assembly only.
my $cds_available = "no"; #keep this as yes if cds files are available. Set as "no" if you only want to mine the mrna files.
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
my $automate_download = "no";
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

#######################################
# 1. Check input files and parameters:
#######################################

######################################
# 1.1. Check variables are specified:
######################################

#1.1.1. yes/no options
my @yes_no_scalar = ();
unless($annotation_available =~ m/^yes$/i || $annotation_available =~ m/^no$/i){
    print "The \$annotation_available parameter is not set correctly. Please specify as \"Yes\" or \"No\" and retry. Aborting job!\n";
    push(@yes_no_scalar, $annotation_available);
}
unless($cds_available =~ m/^yes$/i || $cds_available =~ m/^no$/i){
    print "The \$cds_available parameter is not set correctly. Please specify as \"yes\" or \"no\" and retry. Aborting job!\n";
    push(@yes_no_scalar, $cds_available);
}
unless($predict_new_hits =~ m/^yes$/i || $predict_new_hits =~ m/^no$/i){
    print "The \$predict_new_hits parameter is not set correctly. Please specify as \"yes\" or \"no\" and retry. Aborting job!\n";
    push(@yes_no_scalar, $predict_new_hits);
}
unless($pseudogene_check =~ m/^yes$/i || $pseudogene_check =~ m/^no$/i){
    print "The \$pseudogene_check parameter is not set correctly. Please specify as \"yes\" or \"no\" and retry. Aborting job!\n";
    push(@yes_no_scalar, $pseudogene_check);
}
unless($default_nhmmer_evalue =~ m/^yes$/i || $default_nhmmer_evalue =~ m/^no$/i){
    print "The \$default_nhmmer_evalue parameter is not set correctly. Please specify as \"yes\" or \"no\" and retry. Aborting job!\n";
    push(@yes_no_scalar, $default_nhmmer_evalue);
}
unless($default_phmmer_evalue =~ m/^yes$/i || $default_phmmer_evalue =~ m/^no$/i){
    print "The \$default_phmmer_evalue parameter is not set correctly. Please specify as \"yes\" or \"no\" and retry. Aborting job!\n";
    push(@yes_no_scalar, $default_phmmer_evalue);
}
if(scalar(@yes_no_scalar) > 0){
    die;
}


#1.1.2. Output files and augustus variable
my @check_output_vars = ();
unless($phmm_profile){
    print "Did you forget to specify the \$phmm_profile variable? Aborting job!\n";
    push(@check_output_vars, $phmm_profile);
}
unless($nhmm_profile){
      print "Did you forget to specify the \$nhmm_profile variable? Aborting job!\n";
      push(@check_output_vars, $nhmm_profile);
}
unless($final_cds_nuc){
    print "Did you forget to specify the \$final_cds_nuc variable? Aborting job!\n";
    push(@check_output_vars, $final_cds_nuc);
}
unless($final_cds_prot){
    print "Did you forget to specify the \$final_cds_prot variable? Aborting job!\n";
    push(@check_output_vars, $final_cds_prot);
}
unless($nucleotide_longest_transcripts){
    print "Did you forget to specify the \$phmm_profile variable? Aborting job!\n";
    push(@check_output_vars, $phmm_profile);
}
unless($cds_nucleotide_seqfile){
    print "Did you forget to specify the \$nucleotide_longest_transcripts variable? Aborting job!\n";
    push(@check_output_vars, $nucleotide_longest_transcripts);
}
unless($cds_protein_seqfile){
    print "Did you forget to specify the \$cds_protein_seqfile variable? Aborting job!\n";
    push(@check_output_vars, $cds_protein_seqfile);
}
unless($nhmmer_unnanotated_seqfile){
    print "Did you forget to specify the \$nhmmer_unnanotated_seqfile variable? Aborting job!\n";
    push(@check_output_vars, $nhmmer_unnanotated_seqfile);
}
unless($augustus_species){
    print "Did you forget to specify the \$augustus_species variable? Aborting job!\n";
    push(@check_output_vars, $augustus_species);
}
if(scalar(@check_output_vars) > 0){
    die;
}
if($phmm_profile eq $nhmm_profile){
    print "\$phmm_profile cannot be assigned the same name as \$nhmm_profile. Please rename variables and try again. Abort job!\n";
    die;
}

#1.1.3. Numeric values
my @numeric_check = ();
unless(looks_like_number($minidentity)){
    print "The \$minidentity variable is not numeric. Please ensure a numeric value is set for this variable and try again. Abort job!\n";
    push(@numeric_check, $minidentity);
}
unless(looks_like_number($pseudogene_length)){
    print "The \$pseudogene_length variable is not numeric. Please ensure a numeric value is set for this variable and try again. Abort job!\n";
    push(@numeric_check, $pseudogene_length);
}
unless(looks_like_number($nhmmer_plus)){
    print "The \$nhmmer_plus variable is not numeric. Please ensure a numeric value is set for this variable and try again. Abort job!\n";
    push(@numeric_check, $nhmmer_plus);
}
unless(looks_like_number($nhmmer_minus)){
    print "The \$nhmmer_minus variable is not numeric. Please ensure a numeric value is set for this variable and try again. Abort job!\n";
    push(@numeric_check, $nhmmer_minus);
}
if($default_phmmer_evalue =~ m/^no$/i){
    unless(looks_like_number($phmmer_evalue)){
	print "The \$phmmer_evalue variable is not numeric. Please ensure a numeric value is set for this variable and try again. Abort job!\n";
	push(@numeric_check, $phmmer_evalue);
    }
}
if($default_nhmmer_evalue =~ m/^no$/i){
    unless(looks_like_number($nhmmer_evalue)){
	print "The \$nhmmer_evalue variable is not numeric. Please ensure a numeric value is set for this variable and try again. Abort job!\n";
	push(@numeric_check, $nhmmer_evalue);
    }
}

if(scalar(@numeric_check > 0)){
    die;
}

#1.1.4. Pseudogene names
if($pseudogene_check =~ m/^yes$/i){
    unless($annotation_short && $annotation_1 && $annotation_2 && $annotation_3 && $annotation_4 && $annotation_5 && $annotation_6){
	print "The annotation variables are not set correctly. Please ensure all \$annotation variables are specified and try again. Abort job!\n";
	die;
    }
}


################################
# 1.2. Check input files exist:
################################

unless(-e $pfam_seed){
    print "$pfam_seed does not exist! Did you specify the \$pfam_seed variable? Aborting job!\n";
    die;
}
unless( -e $nuc_alignment){
    print "$nuc_alignment does not exist! Did you specify the \$nuc_alignment variable? Aborting job!\n";
    die;
}
unless(-e $reference_file){
    print "$reference_file does not exist! Did you specify the \$reference_file variable? Aborting job!\n";
    die;
}
if($automate_download eq "Yes" || $automate_download eq "yes"){
    unless(-e $species_list){
	print "$species_list does not exist! Did you specify the \$species_list variable? Aborting job!\n";
	die;
    }
}
else{
    my @check_suffix = ();
    my @suffix_scalar = ();
    push(@check_suffix, $genome_suffix);
    unless($annotation_available =~ m/^no$/i){
	push(@check_suffix, $nt_transcript_suffix);
	push(@check_suffix, $prot_transcript_suffix);
	push(@check_suffix, $cds_suffix);
    }
    foreach my $suffix_value(@check_suffix){
	my @matching_files = glob("*$suffix_value");
	unless(@matching_files) {
	    print "No files ending in $suffix_value exist in working directory. Please check file names and ensure suffix values are correctly assigned in the script. Aborting job!\n";
	    push (@suffix_scalar, $suffix_value);
	}
    }
    if(scalar(@suffix_scalar) > 0){
	die;
    }
}


#######################################
# 2. Automatically download files:
#######################################

#If automate download is set to yes, automatically download genome files
if ($automate_download eq "Yes" || $automate_download eq "yes"){
    my @status =&downloadGenomes($species_list);
    print "\n";
    if(scalar(@status) > 0){
	print "ERROR: Failed to download files automatically. Please manually download the files listed below and set the \$automate_download variable to \"no\". Aborting job!\n";
	foreach my $f(@status){
	    print "Failed to download $f\n";
	}
	die;
    }
    else{
	print "Files downloaded successfully! \n";
    }
}


##############################
# 3.0 Prepare files
#############################

#Read in genome, nucleotide and protein files: 
my @genomes=(<*$genome_suffix>); #read in multiple genome names
my @nucleotide_transcripts =(<*$nt_transcript_suffix>);
my @protein_transcripts = (<*$prot_transcript_suffix>);
my @cds_transcripts = ();
if ($cds_available eq "yes" || $cds_available eq "Yes"){
    @cds_transcripts = (<*$cds_suffix>);
}

#HMMER files
my $phmmer_out = "_phmmer.out"; #genome name will be appedned to this file within code so file will look like: genome_phmmer.out
my $nhmmer_out ="_nhmmer.out"; #nhmmer outfile
#Build hmms:
if ($annotation_available eq "yes" || $annotation_available eq "Yes"){
    `hmmbuild $phmm_profile $pfam_seed`;
}
`hmmbuild $nhmm_profile $nuc_alignment`; 


####################################
# 4. Run TFAM:
####################################

foreach my $genome(@genomes){
    my @genome_files = ();
    my @genome_IDs = ();
    
    ####################################################
    # 4.1. Prepare output directories and  file names:
    ####################################################

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
	    if ($annotation_available =~ m/^yes$/i){
		my $nucleotide_ID = "";
		my $protein_ID = "";
		my $cds = "";

		# RNA file
		foreach my $nucleotide(@nucleotide_transcripts){
		    if ($nucleotide =~ m/([\S]+).*\_$nt_transcript_suffix/){
			$nucleotide_ID = $1;
			if ($genome_ID eq $nucleotide_ID){
			    #print $nucleotide."\n";
			    push (@genome_files, $nucleotide);
			}
		    }
		}

		# Protein file
		foreach my $protein(@protein_transcripts){
		    if ($protein =~ m/([\S]+).*\_$prot_transcript_suffix/){
			$protein_ID = $1;
			if ($genome_ID eq $protein_ID){ 
			   # print $protein."\n";
			    push(@genome_files, $protein);
			}
		    }
		}

		# CDS file
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
		
		############################################
		# 4.2. PHMMER on protein annotation files:
		############################################
		
		my $nucleotide = $genome_files[0];
		my $protein = $genome_files[1];
		print $nucleotide."\n";
		print $protein."\n";
		if ($cds_available =~ m/^yes$/i){
		    $cds = $genome_files[2];
		    print $cds."\n";   
		}
		
		##############
		# Run phmmer #
		##############
		
		print "\nrunning phmmer on protein annotations...\n\n";
		if($default_phmmer_evalue =~ m/^yes$/i){
		    `hmmsearch $phmm_profile $protein >> $phmmer_file`;
		}
		else{
		    `hmmsearch --incE $phmmer_evalue $phmm_profile $protein >> $phmmer_file`;
		}

		################################################
		# Index the protein transcripts for esl-sfetch #
		################################################
		
		`esl-sfetch --index $protein`;
		
		########################
		# Parse phmmer results #
		########################
		
		my @phmmer_hits =&parse_hmmer($phmmer_file);

		
		##################################################
		# Write out hits to file (protein: all isoforms) #
		##################################################
		
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

		
		##############################################################
		# 4.3. TBLASTN: Protein --> RNA
		# Pull corresponding nucleotide rna sequences
		# Database: rna annotations
		# Query: target protein annotations identified with PHMMER
		##############################################################

		###############
		# run tblastn #
		###############
		
		`makeblastdb -in $nucleotide -dbtype="nucl" -out $transcript_db`;
		`tblastn -db $transcript_db -query $phmmer_prot_isoforms -out $transcript_blastout`;

		#############################
		# Index the rna annotations #
		#############################
		
		`esl-sfetch --index $nucleotide`; #index nucleotide transcripts

		########################
		# Parse tblastn output #
		########################
	
		my ($ref_array1, $ref_array2) =&parse_tblastn($transcript_blastout);
		
		# Dereference the arrays
		my @transcripts = @$ref_array1;
		my @locs = @$ref_array2;
	    
		
		#########################################################
		# Write out rna annotations to file (rna: all isoforms) #
		#########################################################

		foreach my $trans_id(@transcripts){
		    my $cmd = "esl-sfetch $nucleotide $trans_id >> $transcript_nucleotide_isoforms";
		    #print $cmd."\n";
		    `esl-sfetch $nucleotide $trans_id >> $transcript_nucleotide_isoforms`;
		}

		##########################################################################################
		# 4.4. NHMMER on rna transcripts
		# Safety net, pick up any additional hits which were missed using phmmer or using tblastn 
		##########################################################################################

		##############
		# Run NHMMER #
		##############

		print "running nhmmer on mRNA transcripts ...\n\n";
		if($default_nhmmer_evalue =~ m/^yes$/i){
		    `nhmmer $nhmm_profile $nucleotide >> $nhmmer_transcript_file`;
		}
		else{
		    `nhmmer --incE $nhmmer_evalue $nhmm_profile $nucleotide >> $nhmmer_transcript_file`;
		}

		########################
		# Parse NHMMER results #
		########################
	        
		my @nhmmer_t_hits =&parse_hmmer($nhmmer_transcript_file);
	

		#################################################################
		# Write out any new rna transcripts to file (rna, all isoforms) #
		#################################################################
		
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

		#############################################
		# 4.5. Pull longest isoform per locus       #
		# RNA --> CDS nucleotide --> CDS protein    #
		#############################################

		###################################################
		# Parse fasta (convert to parsefasta subroutine) #
		###################################################

		#pull only unique locs, longest transcript
	
		my @transcript_seqs =&parse_fasta($transcript_nucleotide_isoforms);

		
		foreach my $LOC(@locs){
		    my @uniquelocs = ();
		    #print $LOC."-----------------------\n";
		    my %transcript_lengths;
		    foreach my $tseq(@transcript_seqs){
			if($tseq =~ m/\>(.*)\n([\S\n]+)/){
			    my $theader = ">".$1;
			    my $tseq = $2;
			    my $tID ="";
			    my $locID ="";
			    
			    ##########################
			    #  Get sequence length   #
			    ##########################
			    
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

		    ########################################
		    #Sort hash to order isoforms by length #
		    ########################################
		    
		    foreach my $traID(sort { $transcript_lengths{$a} <=> $transcript_lengths{$b} or $a cmp $b } keys %transcript_lengths){
			#print $transcript_lengths{$traID}."\n";
			push(@uniquelocs, $traID);
		    }
		    
		    #######################################################################
		    # Pull the longest isoform and write to file (rna, longest isoforms)  #
		    #######################################################################
		    
		    my $longest_transcript = pop(@uniquelocs);
		    `esl-sfetch $nucleotide $longest_transcript >> $unique_longest_transcripts_out`;
		}


		################################################################
		# Blast longest rna transcripts to assembly to get coordinates #
		################################################################

		##########################################################
		# Blast longest rna transcripts against genome assembly  #
		##########################################################
		
		`makeblastdb -in $genome -dbtype="nucl" -out $nblast_db`; #make blast database using query genome
		`blastn -db $nblast_db -query $unique_longest_transcripts_out -out $nblastout`; #blast assembly with nhmmer results

		###########################################
		# Parse blastn results to get coordinates #
		###########################################

		my @block_coordinates =&get_blastn_coords($nblastout);
	       

		#######################
		# mine cds and output #
		#######################
		
		###############################################################
		# Blast longest rna transcripts against cds_from_genomic file #
		###############################################################
		
		if ($cds_available =~ m/^yes$/i){
		    `makeblastdb -in $cds -dbtype="nucl" -out $nblast_cds_db`;
		    `blastn -db $nblast_cds_db -query $unique_longest_transcripts_out -out $nblast_cds`;
		    `esl-sfetch --index $cds`;

		    ########################
		    # Parse blastn results #
		    ########################

		    my @identifiers_cds =&parse_blast($nblast_cds);
		    
		    
		    #########################################################
		    # Write CDS annotations to file (CDS, longest isoforms) #
		    #########################################################

		    foreach my $identifier_name(@identifiers_cds){
			`esl-sfetch $cds \"$identifier_name\" >> $cds_nuc`;
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

				    #################################################################
				    # Use CDS identifiers to pull corresponding proteins            #
				    # Write protein annotations to file (protein, longest isoforms) #
				    #################################################################
				    
				    `esl-sfetch $protein $protein_cds >> $cds_prot`;
				    push(@identifiers_protein_cds, $protein_cds);
				}
			    }
			}
		    }
		}
		close CDS;

		####################################################################
		# 4.6. NHMMER on genome assembly to predict novel unannotated hits #
		####################################################################
		
		#nhmmer on assembly to discover new hits
		#nhmmer then discount anything which overlaps with existing coordinates


		##########################
		# Run nhmmer on assembly #
		##########################
		
		print "running nhmmer on whole assembly to pull new hits ...\n\n";
		`esl-sfetch --index $genome`; #index genomefile
		my @nhmmer_coordinates = ();
		if($default_nhmmer_evalue eq "yes"){
		    `nhmmer $nhmm_profile $genome >> $nhmmer_file`;
		}else{
		    `nhmmer --incE $nhmmer_evalue $nhmm_profile $genome >> $nhmmer_file`;
		}

		########################
		# Parse nhmmer results #
		########################
		
		my @nhmmer_hits=&parse_hmmer($nhmmer_file);

		
		####################
		# Get coordinates  #
		####################
		
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

		#####################################################
		# Infer forward/reverse strand based on coordinates #
		#####################################################
		
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

		    ###################################################
		    # Check if hit overlaps with existing annotations #
		    # Only report as 'novel' if no overlap            #
		    ###################################################
		    
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
		    
		    ####################################################################
		    # If novel hit (not overlapping), print nucleotide range to file   #
		    # Also print domain range to file                                  #
		    ####################################################################
		    
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

			########################################################
			# Reverse strand: Esl-sfetch range and output to files #
			########################################################
			
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

			#########################################################
			# Positive strand: Esl-sfetch range and output to files #
			#########################################################
			
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

		###################################################
		# Make outout directories and tidy existing files #
		###################################################
		
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

	    ############################################################
	    # If no annotations available, run nhmmer on assembly      #
	    # No need to check for overlap as no existing annotations  #
	    ############################################################
	    
	    else{
		
		##########################
		# Run NHMMER on assembly #
		##########################
		
		print "running nhmmer on whole assembly to pull hits ...\n\n";
		`esl-sfetch --index $genome`;
		my @nhmmer_coordinates = ();
		if($default_nhmmer_evalue =~ m/^yes$/i){
		    `nhmmer $nhmm_profile $genome >> $nhmmer_file`;
		}
		else{
		    `nhmmer --incE $nhmmer_evalue $nhmm_profile $genome >> $nhmmer_file`;
		}

		########################
		# Parse NHMMER results #
		########################
		
		my @nhmmer_hits =&parse_hmmer($nhmmer_file);

		
		####################
		# Get coordinates  #
		####################
		
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

		
		#####################################################
		# Infer forward/reverse strand based on coordinates #
		#####################################################
		
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

		    ####################################################################
		    # Print nucleotide range to file                                   #
		    # Also print domain range to file                                  #
		    ####################################################################
		    
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

		    ########################################################
		    # Reverse strand: Esl-sfetch range and output to files #
		    ########################################################
		    
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

		    ########################################################
		    # Forward strand: Esl-sfetch range and output to files #
		    ########################################################
		    
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

		###################################################
		# Make outout directories and tidy existing files #
		###################################################
		
		my $subdir = $outdir."/hmmer_files";
		`mkdir $subdir`;
		`mv $nhmmer_file $subdir`;
		`rm *ssi`;
	    }

	    #######################################
	    # 4.7. Predict new hits with AUGUSTUS #
	    #######################################
	    
	    if ($predict_new_hits eq "Yes" || $predict_new_hits eq "yes"){
		print "running augustus to predict new unannotated hits ...\n\n";

		###################################################################
		# Index the reference file, this is used to guide gene prediction #
		###################################################################
		
		`esl-sfetch --index $reference_file`;

		##################################################
		# Parse fasta: nucleotide +/- range for augustus #
		##################################################
		
		my $sequences = "";
		my $seq_db_pre = $genome_ID."_";
		my $hit_no = 0;
	
		my @seqs =&parse_fasta($nhmmer_nucleotide_sequences);
		
		##################################
		# Get domain ranges from outfile #
		##################################
		
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

		######################################################################
		# Iterate over each novel hit and feed into augustus for predictions #
		######################################################################
		
		foreach my $seq(@seqs){
		    if($seq =~ m/\>(.*)\n([\S\n]+)/){
			my $seq_header = $1;
			my $seq_nt = $2;

			###########################
			# Get domain coordinates  #
			###########################
			
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

			########################### 
			# Prepare sequence header #
			###########################
			
			$seq_header =~ s/\//\_/g;  #NC_044377.1/77156829-77158035 Cannabis sativa chromosome 6,
			$seq_header =~ s/\-/\_/g;
			$seq_header = ">Hit".$hit_no."_".$seq_header;
			my $tmp_out = "tmp.fa";


			#######################################################
			# Print hit +/- range to file and feed into augustus  #
			#######################################################
			
			open(OUT, ">$tmp_out");
			print OUT $seq_header."\n".$seq_nt;
			close OUT;

			############################
			# Prepare augustus files   #
			############################
			my $reference = "";
			my $psl = "Hit".$hit_no."_ref.psl";
			my $hints = "Hit".$hit_no."_hints.gff";
			my $reference_outseq = "Hit".$hit_no."_ref.fa";
			my $prediction_gff = "Hit".$hit_no."prediction_out.gff";

			##########################################################################
			# Use blat and blat2hints to generate hints for augustus gene prediction #
			##########################################################################
			`blat -minIdentity=$minidentity $tmp_out $reference_file $psl`;
			`perl blat2hints.pl --in=$psl --out=$hints`;

			#################
			# Run AUGUSTUS  #
			#################
			
			system("augustus --species=$augustus_species --strand=forward --codingseq=on --softmasking=0 --hintsfile=$hints --extrinsicCfgFile=extrinsic.ME.cfg $tmp_out > $prediction_gff");

			##############################
			# Parse AUGUSTUS predictions #
			##############################
			
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

			    ##########################
			    # Get prediction details #
			    ##########################
			    
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

			    #####################################################################
			    # Check that prediction overlaps with target DOMAIN (NHMMER COORDS) #
			    #####################################################################
			    
			    if($pred_end > $domain_end){
				$maximum_end = $pred_end;
				# print "Maximum end is pred end: $pred_end \n";
			    }else{
				$maximum_end = $domain_end;
				# print "Maximum end is domain end: $domain_end \n";
			    }
			    if($pred_start < $domain_start){
				$minimum_start = $pred_start;
				# print "Minimum start is pred start: $pred_start \n";
			    }
			    else{
				$minimum_start = $domain_start;
				# print "Minimum start is domain start: $domain_start \n";
			    }
			    $prediction_length = ($pred_end - $pred_start) +1;
			    $total_length3 = $prediction_length + $domain_length;
			    $max_span3 = ($maximum_end - $minimum_start) +1;
			    if ($max_span3 < $total_length3){
				# print "we want this prediction!\n";
				# print $gd."\n";
				$prediction_found = 1;
				my $hit_details = "Hit".$hit_no." prediction:\n";
				`echo \"-----------------------------\n\" >> predictions_log.gff`;
				`echo \"$hit_details\" >> predictions_log.gff`;
				`echo \"$pred\" >> predictions_log.gff`;
			    }
			    else{
				$prediction_found = 0;
			    }
			    $protein_details =~ s/\n//g;
			    $cds_details =~ s/\n//g;

			    #################################################################
			    # If prediction overlaps with NHMMER coordinates, print to file #
			    #################################################################

			    if ($prediction_found  == 1){
				my $cds_header = ">Hit".$hit_no."_new_augustus_prediction\n";
				my $hit_annotation = "Hit".$hit_no;

				#######################
				# CDS nucleotide file #
				#######################
				
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

				#######################
				# CDS protein file #
				#######################
				
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
		
		##############################
		# Remove the temporary files #
		##############################
		
		`rm Hit* *ssi tmp.fa`;
		my @domain_seqs = "";

		######################################
		# Print out domain sequences to file #
		######################################
		
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

	    ###############################################################
	    # 4.8. Assign pseudogene/functional status to all predictions #
	    ###############################################################
	    
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

	    ################################
	    # 4.9. Sort output directories #
	    ################################
	    
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

##################################
# 5.0 CLEAN up working directory #
##################################

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


################
# SUBROUTINES: #
################

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


sub downloadGenomes {
    my $splist = $_[0];
    print "\nAutomatically downloading genome files from ncbi for your query species ... \n\n";
    open(SPECIESFILE, $splist);
    my @species = <SPECIESFILE>;
    close SPECIESFILE;
    chomp(@species);
    my $assembly_summary = "assembly_summary_refseq.txt";
    if(-e $assembly_summary){
	`rm $assembly_summary`;
    }
    `wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt`; #RefSeq summary file
    my @status = ();
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
			    $file_cds =~ s/\.gz//g;
			    $file_genome =~ s/\.gz//g;
			    $file_rna =~ s/\.gz//g;
			    $file_prot =~ s/\.gz//g;
			    if(-e $file_cds){
			    }
			    else{
				push(@status, $file_cds);
			    }
			    if(-e $file_genome){
			    }
			    else{
				push(@status, $file_genome);
			    }
			    if(-e $file_rna){
			    }
			    else{
				push(@status, $file_rna);
			    }
			    if(-e $file_prot){
			    }
			    else{
				push(@status, $file_prot);
			    }
			}
		    }
		}
	    }
	}
    }
    return @status;
}


sub parse_hmmer{
    my $hmmer_file = $_[0];
    print "This is hmmer_file: $hmmer_file \n";
    my $hmmer_results;
    open(HMMER, $hmmer_file);
    {
	local $/; #set delimiter to nothing, enables phmmer file to be read in as one chunk
	$hmmer_results = <HMMER>; #store phmmer results
    }
    close HMMER;
    my @hmm_array = split(/\>\>/, $hmmer_results);
    my $hmmer_hit_chunk = $hmm_array[0];
    my @hmmer_array2 = split("Description\n", $hmmer_hit_chunk);
    my $hmmer_hit_chunk2 =  $hmmer_array2[1];
    my @hmmer_hits = ();
    if ($hmmer_hit_chunk2 =~ m/.*inclusion[\s]threshold.*/){
	my @hmmer_array3 = split("------ inclusion threshold ------", $hmmer_hit_chunk2);
	my $sig_hmmer_hits = $hmmer_array3[0];
	@hmmer_hits = split("\n", $sig_hmmer_hits);
	shift(@hmmer_hits); #remove rubbish element 
	pop(@hmmer_hits); #remove empty line at end
    }
    else{
	my @hmmer_array3 = split("\n\nDomain", $hmmer_hit_chunk2);
	my $sig_hmmer_hits = $hmmer_array3[0];
	@hmmer_hits = split("\n", $sig_hmmer_hits);
	shift(@hmmer_hits); #remove rubbish element 1
    }
    return @hmmer_hits;
}

sub parse_fasta{ #returns sequences stored in an array
    my $seqfile = $_[0];
    my @seqs = ();
    open(SEQS, "$seqfile");
    {
	local $/ = ">";
	while(<SEQS>){
	    my $seq = $_;
	    if($seq =~ m/\>/){
		$seq =~ s/\>//g;
		$seq = ">".$seq;
	    }
	    else{
		$seq = ">".$seq;
	    }
	    push (@seqs, $seq);
	}
    }
    close SEQS;
    shift @seqs;
    return @seqs;
}

sub parse_blast{
    my $blast_out_file = $_[0];
    my $blasthits = "";
    my @identifiers = ();
    open(BLASTFILE, $blast_out_file);
    {
	local $/;
	$blasthits = <BLASTFILE>;
    }
    close BLASTFILE;
    my @blast_hits = split(/Query\=/, $blasthits);
    shift @blast_hits;
    foreach my $hit(@blast_hits){
	my @hits = split(/\>/, $hit);
	my $val_1 = $hits[1];
	my $top_hit = ">".$val_1;
	if ($top_hit =~ m/Identities[\s]+\=[\s]+([0-9]+)\/([0-9]+)[\s]+\(([0-9]+).*\)/){
	    my $identity = $3;
	    if ($identity eq 100){
		if ($top_hit =~ m/(\>[\S]+)?\s.*/){
		    my $identifier = $1;
		    $identifier =~ s/\>//g;
		    unless($identifier ~~ @identifiers){
			push(@identifiers, $identifier);
			
		    }
		}
	    }
	}
    }
    return @identifiers;
}

sub parse_tblastn{
    my $tblastn_outfile = $_[0];
    my $tblastn_results;
    my @transcripts_array = (); #new
    my @locs_array = ();
    open(TBLASTN, $tblastn_outfile);
    {
        local $/; #  change line delimter to nothing - read in file as one chunk
        $tblastn_results = (<TBLASTN>);
    }
    close TBLASTN;
    my @tblastn_hits = split(/Query\=/, $tblastn_results);
    shift @tblastn_hits;
    foreach my $query(@tblastn_hits){
        my @queryhits = split(/\>/, $query);
        for(my $i=0;$i<scalar(@queryhits);$i++){
	    my $tophit =">".$queryhits[$i]; #not element zero, this is the string before first hit
	    if ($tophit =~ m/Identities[\s]+\=[\s]+([0-9]+)\/([0-9]+)[\s]+\(([0-9]+).*\)/){
		my $identity = $3;
		if ($identity eq 100){
		    if ($tophit =~ m/(\>[\S]+)?\s.*/){
			my $trans_id_name = $1;
			$trans_id_name =~ s/\>//g;
			my $loc_name = "";
			if ($tophit =~ m/\(([\S]+)\)/){
			    $loc_name = $1;
			}
			if ($trans_id_name ~~ @transcripts_array){
			}
			else{	
			    #print $trans_id."\n";
			    push @transcripts_array, $trans_id_name; #new
			    unless($loc_name ~~ @locs_array){
				push @locs_array, $loc_name;
			    }
			}
		    }
		}
	    }
	}
    }
    return (\@transcripts_array, \@locs_array);
}


sub get_blastn_coords{
    my $nblast_outfile = $_[0];
    my $nblast_results;
    my @block_coords = ();
    my $contig;
    my $start;
    my $end;
    my $coordinates;
    open(BLASTN, $nblast_outfile);
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
	push @block_coords, $coordinates;
    }
    return @block_coords;
}
