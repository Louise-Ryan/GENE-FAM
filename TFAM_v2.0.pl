#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Scalar::Util qw(looks_like_number);

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
my $gff_suffix = "genomic.gff";
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
# Threshold for domain cover in augustus predictions
my $domain_cover_threshold = 0.9; #Augustus predictions which fail to cover this percentage of the hmmer identified region will be discarded

###########################################################################
#output files
#hmm profile names
my $phmm_profile = "MADSp.hmm"; #nhmmer profile: use hmmer to build hmm profile from $pfam_seed
my $nhmm_profile = "MADSn.hmm"; #phmmer profile: use hmmer to build hmm profile from $nuc_alignment
#Output sequence file names
my $final_cds_nuc = "_cds_nuc.fa"; #final cds seq file - includes augustus predictions if turned on
my $final_cds_prot = "_cds_prot.fa"; #final cds seq file - includes augustus predictions if turned on
my $nucleotide_longest_transcripts="_longest_isoforms_rna.fa"; #Longest rna transcripts 
my $cds_nucleotide_seqfile="_longest_isoforms_cds_nucleotide.fa"; #Longest cds transcripts (nucleotide)
my $cds_protein_seqfile = "_longest_isoforms_cds_protein.fa"; #Longest cds transcripts (protein)
my $nhmmer_unnanotated_seqfile = "_unannotated_newhits_from_assembly.fa"; #Unannotated hits from nhmmer on assembly
my $hit_prefix = "Ath";

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
my $blat2hints_file = "blat2hints.pl";
if($predict_new_hits =~ m/^yes$/i){
    unless(-e $blat2hints_file){
	print "blat2hints.pl is not in current working directory. Please copy this file to your working directory and try again. Aborting job!\n";
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
my @gff_files = (<*$gff_suffix>);

#HMMER files
my $phmmer_out = "_phmmer.out"; #genome name will be appedned to this file within code so file will look like: genome_phmmer.out
my $nhmmer_out ="_nhmmer.out"; #nhmmer outfile

#Build hmms:
`hmmbuild $phmm_profile $pfam_seed`;
`hmmbuild $nhmm_profile $nuc_alignment`; 


####################################
# 4. Run TFAM:
####################################

foreach my $genome(@genomes){
    my @genome_files = ();
    my @genome_IDs = ();

    # Declare variables
    my @domain_details = ();
    my @nhmmer_evalues = ();
    my %cds_functional_hash;
    
    ####################################################
    # 4.1. Prepare output directories and  file names:
    ####################################################

    unless ($genome =~ m/$cds_suffix/i){
	print $genome."\n";
	if ($genome =~ m/([\S]+).*\_$genome_suffix/){
	    my $genome_ID = $1;

	    # Directories
	    my $wd = getcwd;
	    my $outdir = $wd."/".$genome_ID."_outfiles";
	    `mkdir $outdir`;
	    my $subdir = $outdir."/hmmer_files";
	    my $subdir2 = $outdir."/isoform_files";

	    # All isoforms files
	    my $transcript_nucleotide_isoforms = $genome_ID."_all_isoforms_mRNA.fa"; #mRNA
	    my $cds_all_isoforms = $genome_ID."_all_isoforms_cds_nucleotide.fa"; #nucleotide CDS
	    my $phmmer_prot_isoforms = $genome_ID."_all_isoforms_cds_protein.fa"; #protein CDS


	    # Longest isoforms files
	    my $unique_longest_transcripts_out = $genome_ID.$nucleotide_longest_transcripts; #mRNA 
	    my $cds_nuc = $genome_ID.$cds_nucleotide_seqfile; #CDS (no augustus)
	    my $cds_prot = $genome_ID.$cds_protein_seqfile; # CDS protein (no augustus)
	    my $cds_final_nuc = $genome_ID.$final_cds_nuc; # CDS (with augustus)
	    my $cds_final_prot = $genome_ID.$final_cds_prot; # CDS (with augustus)

	    
	    # Hmmer out files
	    my $phmmer_file = $genome_ID.$phmmer_out; #phmmer out file
	    my $nhmmer_transcript_file = $genome_ID."_transcripts".$nhmmer_out;
	    my $nhmmer_nucleotide_sequences = $genome_ID.$nhmmer_unnanotated_seqfile;
	    my $nhmmer_file = $genome_ID."_assembly".$nhmmer_out; #nhmmer out file
	    my $nhmmer_cds_file = $genome_ID."_cds_nucleotide".$nhmmer_out;
	   
	    # LOG file
	    my @hit_log = ();

	    #######################
	    # Prepare summary TSV #
	    ######################

	    my %details_hash;
	    my $tsv_summary = $genome_ID."_summary.tsv";
	    my $tsv_details = "";

	    $tsv_details.="Name\t";
	    $tsv_details.="NCBI Annotation (Yes/No)\t";
	    $tsv_details.="Augustus Prediction (Yes/No)\t";
	    $tsv_details.="HMM Filter (Pass/Fail)\t";
	    $tsv_details.="Remove as duplicate (Yes/No)\t";
	    $tsv_details.="Contig\t";
	    $tsv_details.="Contig Length\t";
	    $tsv_details.="Strand\t";
	    $tsv_details.="Locus Name\t";
	    $tsv_details.="mRNA Identifier\t";
	    $tsv_details.="mRNA Start Coordinate\t";
	    $tsv_details.="mRNA End Coordinate\t";
	    $tsv_details.="CDS Identifier\t";
	    $tsv_details.="CDS Start Coordinate\t";
	    $tsv_details.="CDS End coordinate\t";
	    $tsv_details.="CDS Length\t";
	    $tsv_details.="Status (Functional/Pseudogene)\t";
	    $tsv_details.="Hmmer Full Sequence Evalue (Protein)\t";
	    $tsv_details.="Hmmer Domain Evalue (Protein)\t";
	    $tsv_details.="nhmmer Evalue (mRNA)\t";
	    $tsv_details.="nhmmer Evalue (CDS)\t";
	    $tsv_details.="nhmmer Evalue (Genomic)\t";
	    $tsv_details.="nhmmer Start (mRNA)\t";
	    $tsv_details.="nhmmer End (mRNA)\t";
	    $tsv_details.="nhmmer Start (CDS)\t";
	    $tsv_details.="nhmmer End (CDS)\t";
	    $tsv_details.="nhmmer Start (Genomic)\t";
	    $tsv_details.="nhmmer End (Genomic)\n";

	    open(TSV_FILE, ">$tsv_summary");
	    print TSV_FILE $tsv_details;
	    close TSV_FILE;
	    
	    
	    ###################################
	    # Prepare annotation file names  #
	    ###################################
	    
	    if ($annotation_available =~ m/^yes$/i){
		my $nucleotide_ID = "";
		my $protein_ID = "";
		my $cds = "";

		# RNA file
		foreach my $nucleotide(@nucleotide_transcripts){
		    if ($nucleotide =~ m/([\S]+).*\_$nt_transcript_suffix/){
			$nucleotide_ID = $1;
			if ($genome_ID eq $nucleotide_ID){
			    push (@genome_files, $nucleotide);
			}
		    }
		}

		# Protein file
		foreach my $protein(@protein_transcripts){
		    if ($protein =~ m/([\S]+).*\_$prot_transcript_suffix/){
			$protein_ID = $1;
			if ($genome_ID eq $protein_ID){
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
				push (@genome_files, $cds);
			    }
			}
		    }
		}

		# GFF file
		foreach my $gff(@gff_files){
		    if ($gff =~ m/([\S]+).*\_$gff_suffix/){
			my $gff_ID = $1;
			if ($genome_ID eq $gff_ID){
			    push (@genome_files, $gff);
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
		my $gff = $genome_files[3];
		print $gff."\n";
	
		
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
		
		
		########################
		# Parse phmmer results #
		########################

		my %protein_hmm_key;

		my @phmmer_hits =&parse_hmmer($phmmer_file);

		my @protein_key_seen = ();
		my @protein_identifiers = ();
		foreach my $phit(@phmmer_hits){
		    #print "This is a phit: $phit \n";
		    $phit =~ s/[\s]+/\|/g;
		    my @phmmdetails = split(/\|/, $phit);
		    shift @phmmdetails;

		    #PHMMER details
		    my @protein_phmm_details = ();
		    my $full_seq_evalue = $phmmdetails[0];
		    my $domain_evalue = $phmmdetails[3];
		    push(@protein_phmm_details, $full_seq_evalue);
		    push(@protein_phmm_details, $domain_evalue);

		    #PHMMER ID
		    my $prot_ID = $phmmdetails[8];

		    #Populate protein phmmer hash
		    unless($prot_ID ~~ @protein_key_seen){
			#push(@protein_key_seen, $prot_ID);
			$protein_hmm_key{$prot_ID} = \@protein_phmm_details;
		    }

		    #Pull protein IDs
		    push(@protein_identifiers, $prot_ID);
		}
		
		
		
		####################################################
		# Convert Protein hits --> mRNA + CDS nucleotide   #
		####################################################
		
		################################
		# Parse GFF and pull gene info #
		################################
		
		my $id_type = "protein_id";
		my $details =&parse_gff($gff, $id_type, \@protein_identifiers);

		my @locus_IDs = ();
		my @mRNA_IDs = ();
		my @CDS_IDs = ();
		
		my %returned_hash = %$details;
		my %info_hash = map { $_ => $returned_hash{$_} } keys %returned_hash;
		
		foreach my $key(keys %info_hash){
		    #print "Values for $key: \n";
		    my @info = split(/\|/, $info_hash{$key});
		    push(@locus_IDs, $info[6]);
		    push(@mRNA_IDs, $info[7]);
		    push(@CDS_IDs, $info[8]);
		    #print join("\n", @info), "\n";
		    #print "=" x 50, "\n";
		}

		
		###################################
		# 4.3. NHMMER on mRNA transcripts #
		###################################

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

		my %mrna_hmm_key;
		
		my @nhmmer_t_hits =&parse_hmmer($nhmmer_transcript_file);

		my @mrna_key_seen = ();
		my @new_mRNAs = ();
		my @mRNA_gff_IDs = ();
		foreach my $trans_hit(@nhmmer_t_hits){

		    #print "This is trans hit: $trans_hit \n";
		    
		    $trans_hit =~ s/[\s]+/\|/g;
		    my @nhmmtdetails = split(/\|/, $trans_hit);
		    shift(@nhmmtdetails);

		    #mRNA nhmmer details
		    my @mrna_nhmm_details = ();
		    my $nmrna_evalue = $nhmmtdetails[0];
		    my $nmrna_start = $nhmmtdetails[4];
		    my $nmrna_stop = $nhmmtdetails[5];

		    push(@mrna_nhmm_details, $nmrna_evalue);
		    if($nmrna_stop < $nmrna_start){
			my $tmp = $nmrna_start;
			$nmrna_start = $nmrna_stop;
			$nmrna_stop = $tmp;
		    }
		    push(@mrna_nhmm_details, $nmrna_start);
		    push(@mrna_nhmm_details, $nmrna_stop);

		    #mRNA ID
		    my $transcript_ID = $nhmmtdetails[3];

		    #Populate mrna nhmmer hash
		    unless($transcript_ID ~~ @mrna_key_seen){
			$mrna_hmm_key{$transcript_ID} = \@mrna_nhmm_details;
		    }

		    #Pull mRNA IDs
		    unless($transcript_ID ~~ @mRNA_IDs){
			#print "new transcript: $transcript_ID \n";
			push(@new_mRNAs, $transcript_ID);
			my $gff_ID = "rna-".$transcript_ID;
			push(@mRNA_gff_IDs, $gff_ID);
		   }
		}

		
		###################################################
		# new mRNA hits --> CDS protein + CDS nucleotide  #
		###################################################

		################################
		# Parse GFF and pull gene info #
		################################

		if(@mRNA_gff_IDs){
		    my $id_type = "Parent";
		    my $mRNA_details =&parse_gff($gff, $id_type, \@mRNA_gff_IDs);
		    
		    my %mRNA_info_hash = %$mRNA_details;
		    foreach my $key(keys %mRNA_info_hash){
			$info_hash{$key} = $mRNA_info_hash{$key}; #merge hashes
			#print "Values for $key: \n";
			my @mRNA_info = split(/\|/, $mRNA_info_hash{$key});
			push(@locus_IDs, $mRNA_info[6]);
			push(@mRNA_IDs, $mRNA_info[7]);
			push(@CDS_IDs, $mRNA_info[8]);
			#print join("\n", @mRNA_info), "\n";
			#print "=" x 50, "\n";
		    }
		}
		
		
		#####################################
		# 4.4. Run NHMMER on CDS nucleotide #
		#####################################

		print "running nhmmer on CDS nucleotide ...\n\n";
		if($default_nhmmer_evalue =~ m/^yes$/i){
		    `nhmmer $nhmm_profile $cds >> $nhmmer_cds_file`;
		}
		else{
		    `nhmmer --incE $nhmmer_evalue $nhmm_profile $cds >> $nhmmer_cds_file`;
		}

		
		############################
		# Parse NHMMER CDS results #
		############################

		my %cds_hmm_key;
 
		my @nhmmer_cds_hits =&parse_hmmer($nhmmer_cds_file);

		my @cds_key_seen = ();
		my @new_CDS_IDs = ();
		foreach my $nCDS_hit(@nhmmer_cds_hits){
		    #print "This is a cds hit: \n$nCDS_hit\n";
		    $nCDS_hit =~ s/[\s]+/\|/g;
		    my @nhmm_cds_details = split(/\|/, $nCDS_hit);
		    shift(@nhmm_cds_details);

		    #CDS nhmmer details
		    my @cds_nhmm_details = ();
		    my $ncds_evalue = $nhmm_cds_details[0];
		    my $ncds_start = $nhmm_cds_details[5];
		    my $ncds_stop = $nhmm_cds_details[6];
		    
		    if($ncds_stop < $ncds_start){
			my $tmp = $ncds_start;
			$ncds_start = $ncds_stop;
			$ncds_stop = $tmp;
		    }
		    
		    push(@cds_nhmm_details, $ncds_evalue);
		    push(@cds_nhmm_details, $ncds_start);
		    push(@cds_nhmm_details, $ncds_stop);

		    #CDS Identifier
		    my $new_CDS_ID = $nhmm_cds_details[4];
		    if($new_CDS_ID =~ m/cds\_([\S]+)/){
			$new_CDS_ID = $1;
			if($new_CDS_ID =~ m/(\_[0-9]+$)/){
			    $new_CDS_ID =~ s/$1//g;
			}
		    }

		    #Populate nucleotide CDS nhmmer hash
		    unless($new_CDS_ID ~~ @cds_key_seen){
			$cds_hmm_key{$new_CDS_ID} = \@cds_nhmm_details;
		    }

		    #Pull CDS Identifiers
		    unless($new_CDS_ID ~~ @CDS_IDs){
			#print "new CDS: $new_CDS_ID \n";
			push(@new_CDS_IDs, $new_CDS_ID);
		    }
		}

		
		################################
		# Parse GFF and pull gene info #
		################################

		if(@new_CDS_IDs){
		    my $id_type = "protein_id"; #CDS Identifiers are same for protein and nucleotide
		    my $nCDS_details =&parse_gff($gff, $id_type, \@new_CDS_IDs);
		    
		    my %nCDS_info_hash = %$nCDS_details;
		    foreach my $key(keys %nCDS_info_hash){
			$info_hash{$key} = $nCDS_info_hash{$key}; #merge hashes
			#print "Values for $key: \n";
			my @nCDS_info = split(/\|/, $nCDS_info_hash{$key});
			push(@locus_IDs, $nCDS_info[6]);
			push(@mRNA_IDs, $nCDS_info[7]);
			push(@CDS_IDs, $nCDS_info[8]);
			#print join("\n", @nCDS_info), "\n";
			#print "=" x 50, "\n";
		    }
		}

		
		################################################################
		# Write out to files:                                          #
		# 1) Protein, all isoforms.                                    #
		# 2) mRNA, all isoforms,                                       #
		# 3) CDS, nuceleotide all isoforms                             #
		################################################################

		my $parse_type = "";

		# Parse mRNA
		$parse_type = "mRNA";
		my @mRNA_all_headers =&mine_seqs($nucleotide, $transcript_nucleotide_isoforms, \@mRNA_IDs, $parse_type);
	

		# Parse protein
		$parse_type = "prot_CDS";
		my @CDS_protein_all_headers =&mine_seqs($protein, $phmmer_prot_isoforms, \@CDS_IDs, $parse_type);
		
		
		# Parse nucleotide
		$parse_type = "nuc_CDS";
		my @CDS_nuc_all_headers =&mine_seqs($cds, $cds_all_isoforms, \@CDS_IDs, $parse_type);


		
		#############################################
		# 4.5. Pull longest isoform per locus       #
		# CDS nucleotide, CDS protein               #
		#############################################

		#################################
		#   Get CDS lengths from file   #
		#   Check CDS functional status #
		#################################
		
		my %cds_length_hash;
		#my %cds_functional_hash;
		my %ncbi_cds_coords;
		
		my %cds_seqs =&parse_fasta_hash($cds_all_isoforms);

		foreach my $key(keys %cds_seqs){
		    my $cds_seq = $cds_seqs{$key};
		    my $cds_header = $key;
		    $cds_seq =~ s/\n//g;
		    $cds_header =~ s/\n//g;

		    ##################
		    # Get CDS length #
		    ##################

		    my $cds_length = length($cds_seq);
		    my $cdsID = "";
		    if($cds_header =~ m/\[protein\_id\=([^\]]+)/){
			$cdsID = $1;
		    }
		    else{
			#print $cds_header." ...no id\n";
		    }
		    $cds_length_hash{$cdsID} = $cds_length;

		    ###################
		    # Get CDS Coords  #
		    ###################

		    # Pull coords from header
		    my @cds_ncbi_coordinates = $cds_header =~ /(\d+)\.\.(\d+)/g;

		    # Sort array
		    @cds_ncbi_coordinates = sort { $a <=> $b } @cds_ncbi_coordinates;

		    #Get the minimum and maximum values
		    my $cds_genomic_start = shift(@cds_ncbi_coordinates);
		    my $cds_genomic_end = pop(@cds_ncbi_coordinates);
		    my $ncbi_genomic_locus = $cds_genomic_start."|".$cds_genomic_end;
		    #print "This is locus: $ncbi_genomic_locus \n";
		    
		    # Store coords in array:
		    $ncbi_cds_coords{$cdsID} = $ncbi_genomic_locus;

		    
		    #Get functional status
		    if($pseudogene_check =~ m/^yes$/i){
			my @status=&checkframe($cds_seq);
			my $stat = $status[0];
			my $annotation = "";
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
			$cds_functional_hash{$cdsID} = $annotation;
		    }
		}

		# Print out lengths and status
		#foreach my $x(keys %cds_length_hash){
		#    print $x.":".$cds_length_hash{$x}."\n";
		#    print $x.":".$cds_functional_hash{$x}."\n";
		#}

		
		###################################
		# Get longest isoforms per locus  #
		###################################

		my @locus_seen = ();
		my @longest_isoforms = ();
		my @longest_mrnas = ();
		foreach my $rnaid(keys %info_hash){
		    my @information = split(/\|/, $info_hash{$rnaid});
		    my $locus = $information[6];
		    my %locus_lengths;
		    unless($locus ~~ @locus_seen){
			push(@locus_seen, $locus);
			#print "========================================\n";
			#print "Looking for longest isoform for $locus ..\n";
			foreach my $mrnaid(keys %info_hash){
			    my @information2 = split(/\|/, $info_hash{$mrnaid});
			    my $locus2 = $information2[6];
			    my $cds_name = $information2[8];
			    if($cds_name){ 
				if($locus eq $locus2){
				    #print $mrnaid."\n";
				    #print $cds_name."\n";
				    #print "CDS name: $cds_name\n";
				    my $cds_length = $cds_length_hash{$cds_name};
				    #print "CDS length: $cds_length\n";
				    $locus_lengths{$mrnaid} = $cds_length;
				    #print $mrnaid.":".$cds_length."\n";
				}
			    }
			}
			my @sorted_keys = sort { $locus_lengths{$a} <=> $locus_lengths{$b} or $a cmp $b } keys %locus_lengths;
			my $longest_id = pop(@sorted_keys);
			my $longest_length = $locus_lengths{$longest_id};
			my @longest_details = split(/\|/, $info_hash{$longest_id});
			my $longest_cds = $longest_details[8];
			push(@longest_mrnas, $longest_id);
			push(@longest_isoforms, $longest_cds);
			#print "Longest id is $longest_id with length of $longest_length\n";
		    }
		}

		#########################################
		# Print longest isoforms to files:      #
		# 1) mRNA, longest isoforms             #
		# 2) CDS nucleotide, longest isoforms   #
		# 3) CDS protein, longest isoforms      #
		#########################################
	 
		#parse mRNA
		$parse_type = "mRNA";
		my @mRNA_longest_headers =&mine_seqs($transcript_nucleotide_isoforms, $unique_longest_transcripts_out, \@longest_mrnas, $parse_type);
	
		# Parse protein
		$parse_type = "prot_CDS";
		my @CDS_protein_longest_headers =&mine_seqs($phmmer_prot_isoforms, $cds_prot, \@longest_isoforms, $parse_type);
		
		
		# Parse nucleotide
		$parse_type = "nuc_CDS";
		my @CDS_nuc_longest_headers =&mine_seqs($cds_all_isoforms, $cds_nuc, \@longest_isoforms, $parse_type);

		
		############################################################
		# Get mRNA Coordinates                                     #
		# Allows check for overlap when searching for novel hits   #
		############################################################

		my @contig_IDs = ();
		my @block_coordinates = ();
		foreach my $mrnaid(@longest_mrnas){
		    my @details = split(/\|/, $info_hash{$mrnaid});
		    my $contig = $details[0];
		    unless($contig ~~ @contig_IDs){
			push(@contig_IDs, $contig);
		    }
		    my $mrna_start = $details[2];
		    my $mrna_end = $details[3];
		    my $coordinates = $contig."|".$mrna_start."|".$mrna_end;
		    #print $coordinates."\n";
		    push (@block_coordinates, $coordinates);
		}
		

		######################
		# Get contig lengths #
		######################

		my $contig_hash_ref =&get_contig_lengths($genome);
		my %contig_lengths = %$contig_hash_ref;
		
		#foreach my $contig_key(keys %contig_lengths){
		#    print $contig_key.":".$contig_lengths{$contig_key}."\n";
		#}


		###################################
		# PRINT ANNOTATION SUMMARY TO TSV #
		###################################

		open(TSV, ">>$tsv_summary");

		# Sort the hash for consistency
		my @sorted_ncbi_keys = sort keys %info_hash;
		
		# Iterate through the sorted keys
		foreach my $ncbi_hit (@sorted_ncbi_keys) {
		
	
		#foreach my $ncbi_hit(keys %info_hash){
		    my $tsv_entry;
		    my @ncbi_hit_info = split(/\|/, $info_hash{$ncbi_hit});

		    # info hash
		    my $contig_ID = $ncbi_hit_info[0];
		    my $strand = $ncbi_hit_info[1];
		    my $mRNA_start = $ncbi_hit_info[2];
		    my $mRNA_end = $ncbi_hit_info[3];
		    my $ncbi_locus_name = $ncbi_hit_info[6];
		    my $ncbi_mRNA_ID = $ncbi_hit_info[7];
		    my $ncbi_CDS_ID = $ncbi_hit_info[8];

		    # Only print to file if longest isoform
		    if($ncbi_mRNA_ID ~~ @longest_mrnas || $ncbi_CDS_ID ~~ @longest_isoforms){
			
			
			# Contig length
			my $contig_length_value = $contig_lengths{$contig_ID};

			# CDS length
			my $cds_length_value = $cds_length_hash{$ncbi_CDS_ID};

			# CDS coords
			my $cds_ncbi_coordinate_vals = $ncbi_cds_coords{$ncbi_CDS_ID};
			my @coord_split = split(/\|/, $cds_ncbi_coordinate_vals);
			my $ncbi_CDS_start = $coord_split[0];
			my $ncbi_CDS_end = $coord_split[1];
			
			# Pseudogene status
			my $functional_value = "NA";
			if($pseudogene_check =~ m/^yes$/){
			    $functional_value = $cds_functional_hash{$ncbi_CDS_ID};
			}
			
			# HMMER details (protein)
			my $protein_full_seq_eval = "";
			my $protein_domain_eval = "";
			if(exists $protein_hmm_key{$ncbi_CDS_ID}){
			    my @hmmer_protein_details = @{$protein_hmm_key{$ncbi_CDS_ID}};
			    $protein_full_seq_eval = $hmmer_protein_details[0];
			    $protein_domain_eval = $hmmer_protein_details[1];
			}
			else{
			    $protein_full_seq_eval = "NA";
			    $protein_domain_eval = "NA";
			}
			
			# nHMMER details (mRNA)
			my $nhmmer_mrna_eval;
			my $nhmmer_mrna_start;
			my $nhmmer_mrna_end;
			if(exists $mrna_hmm_key{$ncbi_mRNA_ID}){
			    my @nhmmer_mrna_details = @{$mrna_hmm_key{$ncbi_mRNA_ID}};
			    $nhmmer_mrna_eval = $nhmmer_mrna_details[0];
			    $nhmmer_mrna_start = $nhmmer_mrna_details[1];
			    $nhmmer_mrna_end = $nhmmer_mrna_details[2];
			}
			else{
			    $nhmmer_mrna_eval = "NA";
			    $nhmmer_mrna_start = "NA";
			    $nhmmer_mrna_end = "NA";
			}
			
			#nHMMER details (CDS)
			my $nhmmer_cds_eval;
			my $nhmmer_cds_start;
			my $nhmmer_cds_end;
			if(exists $cds_hmm_key{$ncbi_CDS_ID}){
			    my @nhmmer_cds_details = @{$cds_hmm_key{$ncbi_CDS_ID}};
			    $nhmmer_cds_eval = $nhmmer_cds_details[0];
			    $nhmmer_cds_start = $nhmmer_cds_details[1];
			    $nhmmer_cds_end = $nhmmer_cds_details[2];
			}
			else{
			    $nhmmer_cds_eval = "NA";
			    $nhmmer_cds_start = "NA";
			    $nhmmer_cds_end = "NA";
			}

			# Prepare line for TSV:

			$tsv_entry.= $ncbi_locus_name."\t";
			$tsv_entry.= "Yes"."\t";
			$tsv_entry.= "No"."\t";
			$tsv_entry.= "NA"."\t";
			$tsv_entry.= "NA"."\t";
			$tsv_entry.= $contig_ID."\t";
			$tsv_entry.= $contig_length_value."\t";
			$tsv_entry.= $strand."\t";
			$tsv_entry.= $ncbi_locus_name."\t";
			$tsv_entry.= $ncbi_mRNA_ID."\t";
			$tsv_entry.= $mRNA_start."\t";
			$tsv_entry.= $mRNA_end."\t";
			$tsv_entry.= $ncbi_CDS_ID."\t";
			$tsv_entry.= $ncbi_CDS_start."\t";
			$tsv_entry.= $ncbi_CDS_end."\t";
			$tsv_entry.= $cds_length_value."\t";
			$tsv_entry.= $functional_value."\t";
			$tsv_entry.= $protein_full_seq_eval."\t";
			$tsv_entry.= $protein_domain_eval."\t";
			$tsv_entry.= $nhmmer_mrna_eval."\t";
			$tsv_entry.= $nhmmer_cds_eval."\t";
			$tsv_entry.= "NA"."\t";
			$tsv_entry.= $nhmmer_mrna_start."\t";
			$tsv_entry.= $nhmmer_mrna_end."\t";
			$tsv_entry.= $nhmmer_cds_start."\t";
			$tsv_entry.= $nhmmer_cds_end."\t";
			

			my $genomic_Ds = "NA";
			my $genomic_De = "NA";

			if($nhmmer_mrna_start ne "NA" && $nhmmer_mrna_end ne "NA"){
			    if($strand eq "Forward"){
				$genomic_Ds = $mRNA_start + $nhmmer_mrna_start;
				$genomic_De = $mRNA_start + $nhmmer_mrna_end;
			    }
			    elsif($strand eq "Reverse"){
				$genomic_De = $mRNA_end - $nhmmer_mrna_start;
				$genomic_Ds = $mRNA_end - $nhmmer_mrna_end;
			    }
			    #print "Genomic domain coords ncbi: \n";
			    #print $genomic_Ds."\n";
			    #print $genomic_De."\n";
			    #if($genomic_Ds < $genomic_De){
				#print "TRUE \n";
			    #}
			    #else{
				#print "FALSE \n";
			    #}
			}
			$tsv_entry.=  $genomic_Ds."\t";
			$tsv_entry.=  $genomic_De."\n";
		    
			# Print line to TSV
			print TSV $tsv_entry;
		    }
		    
		}
		
		
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
		    my $nevalue = $nhmmdetails[0];
		    push(@nhmmer_evalues, $nevalue);
		    my $ncontig = $nhmmdetails[3]; #print $ncontig."\n";
		    my $nstart =  $nhmmdetails[4]; #print $nstart."\n";
		    my $nend =  $nhmmdetails[5]; #print $nend."\n";
		    my $nhmm_hit = $ncontig."|".$nstart."|".$nend;
		    push @nhmmer_coordinates, $nhmm_hit;
		}

		#####################################################
		# Infer forward/reverse strand based on coordinates #
		#####################################################
		
		foreach my $nhmm_locus(@nhmmer_coordinates){
		    my $overlap = "";
		    my $ntmp = "";
		    my $strand = "";
		    my @details = split(/\|/, $nhmm_locus);
		    my $contig_n = $details[0]; #need to add Xnts 
		    my $start_n = $details[1]; #need to add Xnts
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
			my $contig_length = $contig_lengths{$contig_n};
			my $nhmm_dets = $nhmm_locus."|".$contig_length."|".$strand;
			push (@domain_details, $nhmm_dets);
		
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
			    `$cmd`;
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
			    `$cmd`;
			}
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
		    my $nevalue = $nhmmdetails[0];
		    push(@nhmmer_evalues, $nevalue);
		    my $ncontig = $nhmmdetails[3]; #print $ncontig."\n";
		    my $nstart =  $nhmmdetails[4]; #print $nstart."\n";
		    my $nend =  $nhmmdetails[5]; #print $nend."\n";
		    my $nhmm_hit = $ncontig."|".$nstart."|".$nend;
		    push @nhmmer_coordinates, $nhmm_hit;
		}


		######################
		# Get contig lengths #
		######################
		
		my $contig_hash_ref =&get_contig_lengths($genome);
		my %contig_lengths = %$contig_hash_ref;
		
		#foreach my $contig_key(keys %contig_lengths){
		#    print $contig_key.":".$contig_lengths{$contig_key}."\n";
		#}
		
		#####################################################
		# Infer forward/reverse strand based on coordinates #
		#####################################################
		
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
		    
		    my $contig_length = $contig_lengths{$contig_n};
		    my $nhmm_dets = $nhmm_locus."|".$contig_length."|".$strand;
		    push(@domain_details, $nhmm_dets);
		   
		    ########################################################
		    # Reverse strand: Esl-sfetch range and output to files #
		    ########################################################
		    
		    if($strand eq "rev"){
			$start_n -= $nhmmer_plus; #3' end is $start (hence - plus)
			$end_n += $nhmmer_minus; #5' end is $end (hence + minus)
			if ($start_n < 1){
			    $start_n = 1;
			}
			if ($end_n > $contig_length){
			    #print "exceeds contig range, reverting to max contig bounds!\n";
			    $end_n = $contig_length;
			}
			my $range = $start_n."\.\.".$end_n;
			my $cmd = "esl-sfetch -c $range -r $genome $contig_n >> $nhmmer_nucleotide_sequences";
			`$cmd`;
		    }

		    ########################################################
		    # Forward strand: Esl-sfetch range and output to files #
		    ########################################################
		    
		    elsif($strand eq "pos"){
			$start_n -= $nhmmer_minus;
			$end_n += $nhmmer_plus;
			if ($start_n < 1){
			    $start_n = 1;
			}
			if ($end_n > $contig_length){
			   # print "exceeds contig range, reverting to max contig bounds!\n";
			    $end_n = $contig_length;
			}
			my $range = $start_n."\.\.".$end_n;
			my $cmd = "esl-sfetch -c $range $genome $contig_n >> $nhmmer_nucleotide_sequences";
			`$cmd`;
		    }
		}
	
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

	    ######################################
	    # Declare hash to store hit details  #
	    ######################################

	    my %Hit_Hash;

	    
	    if ($predict_new_hits eq "Yes" || $predict_new_hits eq "yes"){
		print "running augustus to predict novel hits ...\n\n";

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
		
	
		######################################################################
		# Iterate over each novel hit and feed into augustus for predictions #
		######################################################################
		
		foreach my $seq(@seqs){

		    #Declare coord variables
		    my $genomic_gene_start;
		    my $genomic_gene_end;
		    my $genomic_cds_start;
		    my $genomic_cds_end;

		    # Name novel hit
		    $hit_no +=1;
		    my $hit_annotation = $hit_prefix.$hit_no;

		    #Begin TSV append
		    my @tsv_details = ("NA") x 28;
		    $tsv_details[0] = $hit_annotation;
		    $tsv_details[8] = $hit_annotation;
		    $tsv_details[9] = $hit_annotation;
		    $tsv_details[12] = $hit_annotation;
		    $tsv_details[1] = "No";	    
		    
			
		    if($seq =~ m/\>(.*)\n([\S\n]+)/){
			my $seq_header = $1;
			my $seq_nt = $2;

			###########################
			# Get domain coordinates  #
			###########################

			# evalue
			my $nhmm_eval = $nhmmer_evalues[0];
			shift @nhmmer_evalues;
			$tsv_details[21] = $nhmm_eval;

			# coordinates
			my $domain_info = $domain_details[0];
			#print "These are domain details for $hit_no : $domain_details[0]\n";
			my @info = split(/\|/, $domain_info);
			my $contig_name = $info[0];
			my $domain_info_start = $info[1];
			my $domain_info_end = $info[2];
			my $contig_max = $info[3];
			my $strand_direction = $info[4];
			if($strand_direction eq "pos"){
			    $strand_direction = "Forward";
			}
			elsif($strand_direction eq "rev"){
			    $strand_direction = "Reverse";
			}
			if($domain_info_start > $domain_info_end){
			    my $bu = $domain_info_start;
			    $domain_info_start = $domain_info_end;
			    $domain_info_end = $bu;
			}
			shift(@domain_details);

			# Append to tsv details
			$tsv_details[5] = $contig_name;
			$tsv_details[6] = $contig_max;
			$tsv_details[7] = $strand_direction;
			$tsv_details[26] = $domain_info_start;
			$tsv_details[27] = $domain_info_end;
		     

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
			my $a1 = $domain_info_start - $nhmmer_minus; #position a
			#print "This is a: $a \n";
			if($a1 < 1){
			    my $x = 1 - $a1;
			   # print "This is x: $x \n";
			    my $y = $nhmmer_minus - $x;
			   # print "This is y $y \n";
			    $domain_scoord = 1 + $y;
			   # print "this is domain scoord \n";
			}
			elsif($a1 >= 1){
			    $domain_scoord = 1 + $nhmmer_minus;
			    #print "This is domain scoord: $domain_scoord \n";
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
			$seq_header = ">".$hit_prefix.$hit_no."_".$seq_header;
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
			my $psl = $hit_prefix.$hit_no."_ref.psl";
			my $hints = $hit_prefix.$hit_no."_hints.gff";
			my $reference_outseq = $hit_prefix.$hit_no."_ref.fa";
			my $prediction_gff = $hit_prefix.$hit_no."prediction_out.gff";

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
			

			my $prediction_found = 0;
			my $verified_cds_seq = "";
			my $verified_prot_seq = "";

			
			foreach my $pred(@predictions){
			    
			    ##########################
			    # Get prediction details #
			    ##########################
			    unless($prediction_found == 1){
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
				my $cds_start = "";
				my $cds_end = "";
				my $maximum_end = "";
				my $minimum_start = "";
				my $total_length3 = "";
				my $max_span3 = "";
				my $prediction_length = "";
				my @prediction_coords = ();
				my @cds_coords = ();
				
				foreach my $gd(@gene_dets){
				    $gd =~ s/[\s]+/\|/g;
				    if ($gd =~ m/\|gene\|([0-9]+)\|([0-9]+)/){
					$pred_start = $1; $pred_end = $2;
					push(@prediction_coords, $pred_start);
					push(@prediction_coords, $pred_end);
				    }
				    elsif($gd =~ m/\|CDS\|([0-9]+)\|([0-9]+)/){
					$cds_start =$1;
					$cds_end = $2;
					push(@cds_coords, $cds_start);
					push(@cds_coords, $cds_end);
				    }
				}
				@cds_coords = sort { $a <=> $b } @cds_coords;
				# Get the minimum and maximum values
				$cds_start = shift(@cds_coords);
				$cds_end  = pop(@cds_coords);    
				push(@prediction_coords, $cds_start);
				push(@prediction_coords, $cds_end);
				
				$pred_start = $prediction_coords[0];
				$pred_end = $prediction_coords[1];

				if($prediction_coords[2]){
				    $cds_start = $prediction_coords[2];
				}
				else{
				    $cds_start = "NA";
				}
				if($prediction_coords[3]){
				    $cds_end = $prediction_coords[3];
				}
				else{
				    $cds_end = "NA";
				}
				

				
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

				#############################
				# Check prediction overlap: #
				#############################
				
				if ($max_span3 < $total_length3){ 
				    ### NEW:
				    my $domain_not_covered = $max_span3 - $prediction_length;
				    my $percentage_domain_not_covered = $domain_not_covered / $domain_length;
				    my $percentage_domain_cover = 1 - $percentage_domain_not_covered;
				    #print "This is percentage-domain-cover: $percentage_domain_cover \n";


				    #######################################################
				    # Ensure overlap is greater than percentage threshold #
				    #######################################################
				    
				    if($percentage_domain_cover >= $domain_cover_threshold){


					# Test somthing
					#print "========================================================\n";
					#print "Coords for $hit_annotation \n";
					#print "This is prediction start: $pred_start \n";
					#print "This is prediction end: $pred_end \n";
					#print "This is CDS start: $cds_start \n";
					#print "This is CDS end: $cds_end \n";
					#print "This is domain start: $domain_start \n";
					#print "This is domain end: $domain_end \n";
					#print "This is genomic domain start: $domain_info_start \n";
					#print "This is genomic domain end: $domain_info_end \n";



					# Convert to TSV info coords:

					#################
					#1) mRNA coords #
					#################
					
					# 1.1 Within prediction coordinates (zero adjusted)

					my $zero_gene_domain_start;
					my $zero_gene_domain_end;
					#my $genomic_gene_start;
					#my $genomic_gene_end;
					my $zero_gene_end = $pred_end - $pred_start;

			
					#######################
					# FORWARD mRNA COORDS #
					#######################
					
					if($strand_direction eq "Forward"){
					    
					#Domain is fully captured within mRNA at 5' end
					    if($pred_start <= $domain_start){
						$zero_gene_domain_start = $domain_start - $pred_start;
						$genomic_gene_start = $domain_info_start - $zero_gene_domain_start;
						#print "mRNA Forward: Condition 1 \n";
						#print "mRNA Condition 1 start: $genomic_gene_start \n";
					    }
					    
					    #Domain is cut off at the 5' end
					    if($pred_start > $domain_start){
						$zero_gene_domain_start = 0;
						my $diff = $pred_start - $domain_start;
						$genomic_gene_start = $domain_info_start + $diff;
						#print "mRNA Forward: Condition 2 \n";
						#print "mRNA Condition 2 start: $genomic_gene_start \n";
					    }
					    
					    #Domain is fully captured within mRNA at 3' end
					    if($pred_end <= $pred_end){
						$zero_gene_domain_end = $domain_end - $pred_start;
						$genomic_gene_end = $genomic_gene_start + $zero_gene_end;
						#print "mRNA Forward: Condition 3 \n";
						#print "mRNA Condition 3 end: $genomic_gene_end \n";
					    }
					    
					    #Domain is cut off at the 3' end
					    if($pred_end > $pred_end){
						$zero_gene_domain_end = $zero_gene_end;
						#$zero_gene_domain_end = $pred_end;
						$genomic_gene_end = $genomic_gene_start + $zero_gene_end;
						#print "mRNA Forward: Condition 4 \n";
						#print "mRNA Condition 4 end: $genomic_gene_end \n";
					    }
					    
					    #Domain does not occur in mRNA
					    # This condition will never be met #
					    if($domain_start > $pred_end){
						$zero_gene_domain_start = 0;
						$zero_gene_domain_end = 0;
						my $x1 = $domain_start - $pred_end;
						$genomic_gene_end = $domain_info_start - $x1;
						my $x2 = $pred_end - $pred_start;
						$genomic_gene_start = $genomic_gene_end - $x2;
					    }
					    
					    #Domain does not occur in mRNA
					    # This condition will never be met #
					    if($domain_end < $pred_start){
						$zero_gene_domain_end = 0;
						$zero_gene_domain_start = 0;
						my $x1 = $pred_start - $domain_end;
						$genomic_gene_start = $domain_info_end + $x1;
						my $x2 = $pred_end - $pred_start;
						$genomic_gene_end = $genomic_gene_start + $x2;
					    }
					}

					#######################
					# REVERSE mRNA COORDS #
					#######################

					if($strand_direction eq "Reverse"){
					    #print "Reverse strand\n";
					    
					    #Domain is fully captured within mRNA at 5' end
					    if($pred_start <= $domain_start){
						$zero_gene_domain_start = $domain_start - $pred_start;
						$genomic_gene_end = $domain_info_end + $zero_gene_domain_start;
						#print "Reverse mRNA, condition 1 \n";
						#print "This is condition 1 end: $genomic_gene_end \n";
					    }

					    #Domain is cut off at the 5' end
					    if($pred_start > $domain_start){
						$zero_gene_domain_start = 0;
						my $diff = $pred_start - $domain_start;
						$genomic_gene_end = $domain_info_end - $diff;
						#print "Reverse mRNA, condition 2 \n";
						#print "This is condition 2 end: $genomic_gene_end \n";
						#print "I am here \n";
					    }
					    
					    #Domain is fully captured within mRNA at 3' end
					    if($domain_end <= $pred_end){
						$zero_gene_domain_end = $domain_end - $pred_start;
						my $x2 = $zero_gene_end - $zero_gene_domain_end;
						$genomic_gene_start = $domain_info_start - $x2;
						$genomic_gene_start -=1;
						#print "Reverse mRNA, condition 3 \n";
						#print "This is condition 3 start: $genomic_gene_start \n";
						
					    }
					    
					    #Domain is cut off at the 3' end
					    if($domain_end > $pred_end){
						#$zero_gene_domain_end = $pred_end;
						$zero_gene_domain_end = $zero_gene_end;
						$genomic_gene_start = $genomic_gene_end - $zero_gene_end;
						#print "Reverse mRNA, condition 4 \n";
						#print "This is condition 4 start: $genomic_gene_start \n";
					    }
					    
					    # Domain does not occur in mRNA
					    # This condition will never be met
					    if($domain_start > $pred_end){
						# same as condition 1
						$genomic_gene_end = $domain_info_end + $zero_gene_domain_start;

						# same as condition 4
						$genomic_gene_start = $genomic_gene_end - $zero_gene_end;

						# zero domain as 0 as not in mRNA
						$zero_gene_domain_start = 0;
						$zero_gene_domain_end = 0;	
					    }
					    
					    #Domain does not occur in mRNA
					    # This condition will never be met
					    if($domain_end < $pred_start){
						# same as condition 2
						my $diff = $pred_start - $domain_start;
						$genomic_gene_end = $domain_info_end - $diff;

						# same as condition 3
						my $x2 = $zero_gene_end - $zero_gene_domain_end;
						$genomic_gene_start = $domain_info_start - $x2;
						$genomic_gene_start -=1;

						# zer domain as 0 as not in mRNA
						$zero_gene_domain_end = 0;
						$zero_gene_domain_start = 0;
					    }
					    
					}
				
					#print "This is zero domain mRNA start: $zero_gene_domain_start \n";
					#print "This is zero domain mRNA end: $zero_gene_domain_end \n";					
					#print "This is genomic GENE start: $genomic_gene_start \n";
					#print "This is genomic GENE end: $genomic_gene_end \n";

					#################
					#2) CDS coords: #
					#################
					
					# 2.1 Within prediction coordinates (zero adjusted)

					my $zero_cds_domain_start;
					my $zero_cds_domain_end;
					#my $genomic_cds_start;
					#my $genomic_cds_end;
					my $zero_cds_end = $cds_end - $cds_start;


					######################
					# FORWARD CDS COORDS #
					######################
					
					if($strand_direction eq "Forward"){
					    #Domain is fully captured within CDS at 5' end
					    if($cds_start <= $domain_start){
						#print "CDS Forward: Condition 1\n";
						$zero_cds_domain_start = $domain_start - $cds_start;
						#print "zero domain start: $zero_cds_domain_start \n";
						$genomic_cds_start = $domain_info_start - $zero_cds_domain_start;
					    }
					    
					    #Domain is cut off at the 5' end
					    if($cds_start > $domain_start){
						#print "CDS Forward: Condition 2\n";
						$zero_cds_domain_start = 0;
						my $diff = $cds_start - $domain_start;
						$genomic_cds_start = $domain_info_start + $diff;
					    }
					    
					    #Domain is fully captured within CDS at 3' end
					    if($domain_end <= $cds_end){
						#print "CDS Forward: Condition 3 \n";
						$zero_cds_domain_end = $domain_end - $cds_start;
						$genomic_cds_end = $genomic_cds_start + $zero_cds_end;
						#print "CDS Condition 3, end: $genomic_cds_end \n";
					    }
					    
					    #Domain is cut off at the 3' end
					    if($domain_end > $cds_end){
						#print "CDS Forward: Condition 4 \n";
						$zero_cds_domain_end = $zero_cds_end;
						$genomic_cds_end = $genomic_cds_start + $zero_cds_end;
						#print "CDS Condition 4, end: $genomic_cds_end \n";
					    }
					    
					    #Domain does not occur in CDS
					    if($domain_start > $cds_end){
						$zero_cds_domain_start = 0;
						$zero_cds_domain_end = 0;
						my $x1 = $domain_start - $cds_end;
						$genomic_cds_end = $domain_info_start - $x1;
						my $x2 = $cds_end - $cds_start;
						$genomic_cds_start = $genomic_cds_end - $x2;
						#print "CDS Forward: Condition 5 \n";
						#print "CDS Condition 5 start: $genomic_cds_start \n";
						#print "CDS Condition 5 end: $genomic_cds_end \n";
					    }
					    
					    #Domain does not occur in CDS
					    if($domain_end < $cds_start){
						$zero_cds_domain_end = 0;
						$zero_cds_domain_start = 0;
						my $x1 = $cds_start - $domain_end;
						$genomic_cds_start = $domain_info_end + $x1;
						my $x2 = $cds_end - $cds_start;
						$genomic_cds_end = $genomic_cds_start + $x2;
						$genomic_cds_start +=1;
						$genomic_cds_end +=1;
						#print "CDS Forward: Condition 6 \n";
						#print "CDS Condition 6 start: $genomic_cds_start \n";
						#print "CDS Condition 6 end: $genomic_cds_end \n";
					    }
					
					    #################
					    # check mRNA
					    #################
					    #my $checker = $cds_start - $pred_start;
					    #if($genomic_cds_start - $checker == $genomic_gene_start){
					    #print "Good start!\n";
					    #}
					    #else{
					    #	print "BAD start\n";
					    #   }
					    #  
					    # #check mRNA
					    #my $checker = $pred_end - $cds_end;
					    #if($genomic_gene_end - $checker == $genomic_cds_end){
					    #	print "Good end!\n";
					    #   }
					    #  else{
					    #	print "BAD end\n";
					    #   }
					    #}
					    
					}
					
					######################
					# REVERSE CDS COORDS #
					######################
					
					if($strand_direction eq "Reverse"){
					    #print "Reverse strand\n";
					    
					    #Domain is fully captured within CDS at 5' end
					    if($cds_start <= $domain_start){
						$zero_cds_domain_start = $domain_start - $cds_start;
						$genomic_cds_end = $domain_info_end + $zero_cds_domain_start;
						#print "Reverse CDS, condition 1 \n";
						#print "This is condition 1 end: $genomic_cds_end \n";
					    }
					    
					    #Domain is cut off at the 5' end
					    if($cds_start > $domain_start){
						$zero_cds_domain_start = 0;
						my $diff = $cds_start - $domain_start;
						$genomic_cds_end = $domain_info_end - $diff;
						#print "Reverse CDS, condition 2 \n";
						#print "This is condition 2 end: $genomic_cds_end \n";
					    }
					    
					    #Domain is fully captured within CDS at 3' end
					    if($domain_end <= $cds_end){
						$zero_cds_domain_end = $domain_end - $cds_start;
						my $x2 = $zero_cds_end - $zero_cds_domain_end;
						$genomic_cds_start = $domain_info_start - $x2;
						$genomic_cds_start -=1;
						#print "Reverse CDS, condition 3 \n";
						#print "This is condition 3 start: $genomic_cds_start \n";
					    }
					    
					    #Domain is cut off at the 3' end
					    if($domain_end > $cds_end){
						$zero_cds_domain_end = $zero_cds_end;
						$genomic_cds_start = $genomic_cds_end - $zero_cds_end;
						#print "Reverse CDS, condition 4 \n";
						#print "This is condition 4 start: $genomic_cds_start \n";
					    }
					    
					    #Domain does not occur in CDS
					    if($domain_start > $cds_end){
						# same as condition 1
						$genomic_cds_end = $domain_info_end + $zero_cds_domain_start;

						# same as condition 4
						$genomic_cds_start = $genomic_cds_end - $zero_cds_end;
						
						# zero domain as 0 as not in CDS
						$zero_cds_domain_start = 0;
						$zero_cds_domain_end = 0;
						#print "This is condition 5 start: $genomic_cds_start \n";
						#print "This is condition 5 end: $genomic_cds_end \n";
					    }
					    
					    #Domain does not occur in CDS
					    if($domain_end < $cds_start){
						# same as condition 2
						my $diff = $cds_start - $domain_start;
						$genomic_cds_end = $domain_info_end - $diff;

						# same as condition 3
						my $x2 = $zero_cds_end - $zero_cds_domain_end;
						$genomic_cds_start = $domain_info_start - $x2;
						$genomic_cds_start -=1;

						# zer domain as 0 as not in CDS
						$zero_cds_domain_end = 0;
						$zero_cds_domain_start = 0;
						
						#print "This is condition 6 start: $genomic_cds_start \n";
						#print "This is condition 6 end: $genomic_cds_end \n";
					    }
					 
					    #check mRNA

					    #my $checker = $pred_end - $cds_end;
					    #if($genomic_cds_start - $checker == $genomic_gene_start){
					    #    print "Good start!\n";
					    #}
					    #else{
					    #    print "BAD start\n";
					    #}
					   # 
					    #my $checker = $cds_start - $pred_start;
					    #if($genomic_gene_end - $checker == $genomic_cds_end){
					    #    print "Good end!\n";
					    #}
					    #else{
					    #    print "BAD end\n";
					    #}

					}
					
					#print "This is zero domain cds start: $zero_cds_domain_start \n";
					#print "This is zero domain cds end: $zero_cds_domain_end \n";
					#print "This is genomic cds start: $genomic_cds_start \n";
					#print "This is genomic cds end: $genomic_cds_end \n";
					
					
					# Prepare CDS header and seq
					$protein_details =~ s/\n//g;
					$cds_details =~ s/\n//g;

					
					my $hit_header = ">Hit".$hit_no."_new_augustus_prediction\n";

					$tsv_details[2] = "Yes";

					
					#######################
					# CDS nucleotide      #
					#######################

					my $augustus_cds_length = "";
					
					my $cds_seq = "";
					if ($cds_details =~ m/\[([^\]]+)\]/){
					    $cds_seq = $1;
					    $cds_seq =~ s/\#//g;
					    $cds_seq =~ s/\s//g;
					    $augustus_cds_length = length($cds_seq);
					    $cds_seq =~ s/.{80}\K/\n/g;
					    $cds_seq = uc($cds_seq);
					}

					
					#######################
					# CDS protein         #
					#######################
					
					my $protein_seq = "";
					if ($protein_details =~ m/\[([^\]]+)\]/){
					    $protein_seq = $1;
					    $protein_seq =~ s/\#//g;
					    $protein_seq =~ s/\s//g;
					    $protein_seq =~ s/.{80}\K/\n/g;
					}
					

					#print $cds_seq."\n";
					#######################################################
					# Prep temporary file with protein sequence for hmmer #
					#######################################################
					
					my $tmp_prot_out = "tmp_prot_out.fa"; 
					open(TPO, ">$tmp_prot_out");
					print TPO $hit_header.$protein_seq."\n";
					close TPO;
					
					######################################
					# Prepare input files for subroutine #
					######################################
					
					my @prediction_details = ();
					push(@prediction_details, $phmm_profile);
					push(@prediction_details, $tmp_prot_out);
					
					####################################
					# Push evaule for hmmer subroutine #
					####################################
					
					if($default_phmmer_evalue =~ m/^yes$/i){
					    my $default = "default";
					    push(@prediction_details, $default);
					}
					else{
					    push(@prediction_details, $phmmer_evalue);
					}


					# Add mRNA coords to tsv
				
					$tsv_details[10] = $genomic_gene_start;
					$tsv_details[11] = $genomic_gene_end;
				
					$tsv_details[22] = $zero_gene_domain_start;
					$tsv_details[23] = $zero_gene_domain_end;
					
					##############
					# HMM filter #
					##############
					
					my $hmm_status =&hmm_filter(@prediction_details);
					
					#print "This is hmm status: $hmm_status \n";
					
					if($hmm_status == 1){
					    $prediction_found = 1;
					    my $hit_details = $hit_prefix.$hit_no." prediction:\n";
					    $verified_cds_seq = $cds_seq;
					    $verified_prot_seq = $protein_seq;
					    `echo \"-----------------------------\n\" >> predictions_log.gff`;
					    `echo \"$hit_details\" >> predictions_log.gff`;
					    `echo \"$pred\" >> predictions_log.gff`;
					    $tsv_details[2] = "Yes";
					    $tsv_details[3] = "Pass";


					    # Update these (HERE HERE HERE)
					    $tsv_details[24] = $zero_cds_domain_start;
					    $tsv_details[25] = $zero_cds_domain_end;

					    #Update these
					    $tsv_details[13] = $genomic_cds_start;
					    $tsv_details[14] = $genomic_cds_end;

					    # length
					    $tsv_details[15] = $augustus_cds_length;
					    
					    # Functional or Pseudogene status
					    #Get functional status
					    if($pseudogene_check =~ m/^yes$/i){
						my $copy_seq = $cds_seq;
						$copy_seq =~ s/\n//g;
						my @status=&checkframe($copy_seq);
						my $stat = $status[0];
						my $annotation = "";
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
						if ($augustus_cds_length < $pseudogene_length){
						    $annotation = $annotation_short;
						}
						$tsv_details[16] = $annotation;
						$cds_functional_hash{$hit_annotation} = $annotation;
					    }
					}
					else{
					    $prediction_found = 0;
					    $tsv_details[3] = "Fail";
					}
				    }
				    else{
					$prediction_found = 0;
				    }
				}
				else{
				    $prediction_found = 0;
				  
				}
			    }
	
			}
			
			#################################################################
			# If prediction overlaps with NHMMER coordinates, print to file #
			#################################################################
			
			if ($prediction_found  == 1){
			    my $locus_flag = " [Coordinates=$contig:$genomic_cds_start-$genomic_cds_end]";
			    my $prediction_flag =" [Augustus prediction]";
			    my $cds_header = ">".$hit_prefix.$hit_no.$locus_flag.$prediction_flag."\n";
			    #print $cds_header."\n";
			    #my $cds_header = ">".$hit_prefix.$hit_no."_new_augustus_prediction\n";
			    #$cds_header = ">".$hit_prefix.$hit_no."_new_augustus_prediction\n";
			    my $hit_annotation = $hit_prefix.$hit_no;

			    my $cds_copy = $verified_cds_seq;
			    $cds_copy =~ s/\n//g;
			    my $cds_length = length($cds_copy);
			    $tsv_details[12] = $cds_length;
			    
			    # Print to CDS nucleotide
			    open(CDS_NT, ">>$cds_final_nuc");
			    print CDS_NT $cds_header.$verified_cds_seq."\n";
			    close CDS_NT;
			    
			    # Print to CDS protein
			    open(CDS_PROT, ">>$cds_final_prot");
			    print CDS_PROT $cds_header.$verified_prot_seq."\n";
			    close CDS_PROT;
			}
			else{
			    if($tsv_details[3] eq "Fail"){
				$tsv_details[2] = "Yes";
			    }
			}
			if($tsv_details[2] eq "NA"){
			    $tsv_details[2] = "No";
			}
			open(TSV, ">>$tsv_summary");
			my $last_entry = pop(@tsv_details);
			print TSV join("\t", @tsv_details), "\t";
			print TSV $last_entry."\n";
			close TSV;
		    }
		}
	    
		##############################
		# Remove the temporary files #
		##############################

		if($hit_prefix){
		    `rm $hit_prefix*`;
		}
		`rm *ssi tmp.fa`;
		my @domain_seqs = "";
		
	    }

	    ###############################################################
	    # 4.8. Assign pseudogene/functional status to all predictions #
	    ###############################################################
	    
	    
	    #open(CDS, $cds_final_nuc);
	    my $out = "tmp_cds.fa";
	    my @cds_final_files = ();
	    my %prot_headers;
	    push(@cds_final_files, $cds_final_prot);
	    push(@cds_final_files, $cds_final_nuc);
	    
	    my $counter = 0;
	    foreach my $cds_final_file(@cds_final_files){
		my @ncbi_seqs_final = ();
		my @augustus_seqs_final = ();
		$counter ++;
		
		#print "This is cds final file: $cds_final_file\n";
		
		if(-e $out){
		    `rm $out`;
		}
		
		my @finalseqs =&parse_fasta($cds_final_file);
		#my @finalseqs = split (/\>/, $final_seqs);
		
		foreach my $cds(@finalseqs){
		    if($cds =~ m/(.*)\n([\S\n]+)/){
			my $cds_header = $1;
			my $cds_seq = $2;
			$cds_header =~ s/\n//g;
			$cds_seq =~ s/\n//g;
			$cds_seq =~ s/.{80}\K/\n/g;
			my $cds_ID_name;
			my $ncbi_prediction = 0;
			if($cds_header =~ m/NCBI/i){
			    $ncbi_prediction = 1;
			}
			if($cds_header =~ m/\[protein\_id\=([^\]]+)/){
			    $cds_ID_name = $1;
			}
			elsif($cds_header =~ m/\>([\S]+)/){
			    $cds_ID_name = $1;
			}

			my $functional_flag = "";
			if ($pseudogene_check eq "yes" || $pseudogene_check eq "Yes"){
			    #print "Annotating hits with functional or pseudogene status ...\n\n";
			    my $functional_status = $cds_functional_hash{$cds_ID_name};
			    $functional_flag = " [status=$functional_status]";
			}
			else{
			    $functional_flag = "NA";
			}
			if($counter ==1){
			    if($functional_flag ne "NA"){
				$cds_header.=$functional_flag;
			    }
			}
			if($counter == 1){ # protein
			    $prot_headers{$cds_ID_name} = $cds_header;
			}
			elsif($counter == 2){ #nucleotide
			    $cds_header = $prot_headers{$cds_ID_name};
			}
			
			my $seq_final = $cds_header."\n".$cds_seq."\n";
			
			if($ncbi_prediction ==1){
			    push(@ncbi_seqs_final, $seq_final);
			}
			else{
			    push(@augustus_seqs_final, $seq_final);
			}
		    }
		}
		my @sorted_ncbi_seqs_final = sort(@ncbi_seqs_final);
		foreach my $annotated_seq(@sorted_ncbi_seqs_final){	
		    open(OUT_FINAL, ">>$out");
		    print OUT_FINAL $annotated_seq;
		    close OUT_FINAL;
		}
		if(@augustus_seqs_final){
		    #my @sorted_augustus_seqs_final = sort(@augustus_seqs_final);
		    foreach my $annotated_seq(@augustus_seqs_final){	
			open(OUT_FINAL, ">>$out");
			print OUT_FINAL $annotated_seq;
			close OUT_FINAL;
		    }
		}
		`mv $out $cds_final_file`;
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
			    my $file_gff = $accession."_".$gff_suffix.".gz";
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
			    my $final_link_gff = $link."/".$file_gff;
			    `wget $final_link_gff`;
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



sub hmm_filter{
    my @hmminputs = @_;
    my $hmm_file = $hmminputs[0];
    my $prediction_file = $hmminputs[1];
    my $evalue_threshold = $hmminputs[2];
    my $hmm_out = "tmp_hmmfilter.out";
    
    if($evalue_threshold =~ m/^default$/i){
	`hmmsearch $hmm_file $prediction_file > $hmm_out`;
    }
    else{
	`hmmsearch --incE $evalue_threshold $hmm_file $prediction_file > $hmm_out`;
    }	
    
    my $results_file;
    open(IN, $hmm_out);
    {
        local $/;
        $results_file = <IN>;
    }
    close IN;
    my $regex = "\-";
    my $expression = $regex x 37;
    my @hmmhits = split(/$expression/, $results_file);
    my $summary = $hmmhits[1];
    my @stats = split(/\n/, $summary);
    my $target_seqs = $stats[2];
    my $passed_hits = $stats[8];
    my $targets = "";
    my $hits = "";
    if($target_seqs){
	if ($target_seqs =~ m/\:[\s]+([0-9]+)\s/){
	    $targets = $1;
	}
    }
    if($passed_hits){
	if ($passed_hits =~ m/\:[\s]+([0-9]+)\s/){
	    $hits = $1;
	}
    }
    my $hmm_status = 0;
    if($targets){
	if($targets == 1){
	    if ($hits ==1){
		$hmm_status = 1;
		return $hmm_status;
	    }
	    else{
		$hmm_status = 0;
		return $hmm_status;
	    }
	}
	else{
	    $hmm_status = 0;
	    return $hmm_status;
	}
    }
    else{
	$hmm_status = 0;
	return $hmm_status;   
    }
}


sub parse_gff{
    #Inputs:
    my ($gff_file, $identifier_type, $identifier_names_ref) = @_;
    my @identifiers = @{$identifier_names_ref}; #dereference identifiers array

    #create array for gene details with NA values
    my @gene_details = ();
    for my $i (1..11) {
	push @gene_details, "NA"; 
    }

    #Declare variables
    my %gene_info;
    my @rna_IDs = ();
    
    #Read in GFF file
    open my $fh, '<', $gff_file or die "Cannot open GFF file: $!";
    {
	local $/ = "gene\t"; #read in file in 'gene' chunks
	while (my $chunk = <$fh>) {
	    chomp $chunk;  # Remove newline character
	    
	    # Skip comment lines (lines starting with '#')
	    next if $chunk =~ /^#/;
		
	    #Split gene chunk into features (mRNA, CDS, exon etc ..)
	    my @features = split("\n", $chunk);
	    shift(@features);
	    foreach my $line(@features){
		my $rna_ID;
		my $protein_ID;
		my $gene_ID;
		
		# Split the line into fields using tab as the delimiter
		my @fields = split /\t/, $line;

		if(scalar(@fields) >= 7){
		# Access information about the feature
		    my $seqname = $fields[0];
		    my $source = $fields[1];
		    my $type = $fields[2];
		    my $start = $fields[3];
		    my $end = $fields[4];
		    my $score = $fields[5];
		    my $strand = $fields[6];
		    my $phase = $fields[7];
		    my $attributes = $fields[8];

		    if($type eq "CDS"){
			my %attribute_hash;
			foreach my $attribute (split /;/, $attributes) {
			    my ($key, $value) = split /=/, $attribute;
			    $attribute_hash{$key} = $value;
			}

			#Check that ID is target
			if($attribute_hash{$identifier_type} ~~ @identifiers){

			    $gene_details[0] = $seqname; #contig name
			    #$gene_details[1] = $type;
			    
			    #Account for reverse strand
			    if($end < $start){
				my $tmp_val = $start;
				$start = $end;
				$end = $tmp_val;
			    }

			    if($strand eq "+"){
				$strand = "Forward";
			    }
			    elsif($strand eq "-"){
				$strand = "Reverse";
			    }
			    
			    $gene_details[1] = $strand; #Forward or Reverse
			    #$gene_details[4] = $start;	#CDS start		    
			    #$gene_details[5] = $end; #CDS end
			    
			    #Gather attributes
			    foreach my $key (keys %attribute_hash) {
				my $value = $attribute_hash{$key};
				
				if($key eq "Parent"){
				    my $adjust_value = $value;
				    $adjust_value =~ s/rna\-//g;
				    $rna_ID = $adjust_value;
				    push(@rna_IDs, $value);
				}
				elsif($key eq "protein_id"){
				    $protein_ID = $value;
				}
			    }

			    #Populate @gene_details
			    if($rna_ID){
				$gene_details[7] = $rna_ID;
			    }
			    if($protein_ID){
				$gene_details[8] = $protein_ID;
			    }
			    $gene_info{$rna_ID} = join("|", @gene_details);
			}
		    }
		}

	    }
	    #Get GENE NAME from mRNA features.
	    foreach my $line(@features){
		my $gene_ID;

		# Split the line into fields using tab as the delimiter
		my @fields = split /\t/, $line;

		if(scalar(@fields) >= 7){
		    # Access information about the feature
		    my $seqname = $fields[0];
		    my $source = $fields[1];
		    my $type = $fields[2];
		    my $start = $fields[3];
		    my $end = $fields[4];
		    my $score = $fields[5];
		    my $strand = $fields[6];
		    my $phase = $fields[7];
		    my $attributes = $fields[8];

		    if($end < $start){
			my $tmp_val = $start;
			$start = $end;
			$end = $tmp_val;
		    }
		    my $mrna_length = ($end - $start) +1;
		    
		    if($type eq "mRNA"){
			my %attribute_hash;
			foreach my $attribute (split /;/, $attributes) {
			    my ($key, $value) = split /=/, $attribute;
			    $attribute_hash{$key} = $value;
			}
			my $id = "ID";
			if($attribute_hash{$id} ~~ @rna_IDs){
			    my $gene_id = $attribute_hash{"Parent"};
			    $gene_id =~ s/gene-//g;
			    my $rna_name = $attribute_hash{$id};
			    my $rna = $rna_name;
			    $rna =~ s/rna-//g;	    
			    my @info_array = split(/\|/, $gene_info{$rna});
			    $info_array[6] = $gene_id; #gene ID
			    $info_array[2] = $start; #mRNA start
			    $info_array[3] = $end; #mRNA end
			    $info_array[9] = $mrna_length; #mRNA length
			    $gene_info{$rna} = join("|", @info_array);
			}
		    }		
		}
	    }
	}
	close $fh;
	return(\%gene_info);
    }
}

sub parse_fasta_hash{ #returns sequences stored in an array
    my $seqfile = $_[0];
    my %seqs = ();
    open(SEQS, "$seqfile");
    my $header;   
    while(<SEQS>){
	my $line = $_;
	if($line =~ m/^>/){
	    $header = $line;
	    chomp $header;
	}
	elsif($line =~ m/[\S]+/){
	    if($seqs{$header}){
		my $seq = $seqs{$header};
		$seq.=$line;
		$seqs{$header} = $seq;
	    }
	    else{
		$seqs{$header} = $line;
	    }
	}
    }
    close SEQS;
    return %seqs;
}

sub parse_fasta{ #returns sequences stored in an array
    my $seqfile = $_[0];
    my @seqs = ();
    open(SEQFILE, "$seqfile");
    {
	local $/ = ">";
	while(<SEQFILE>){
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
    close SEQFILE;
    shift @seqs;
    return @seqs;
}


sub mine_seqs{
    my ($in_file, $out_file, $ID_array_ref, $seq_type) =  @_;
    my @ID_array = @{$ID_array_ref}; #dereference identifiers array

    my $ncbi_flag = " [NCBI Prediction]";
    #Declare return array
    my @headers = ();
    
    # Open the output file
    open(FILEOUT, ">>$out_file");
    
    #Loop through file and search for IDs
    my %seq_hash =&parse_fasta_hash($in_file);
    foreach my $key(keys %seq_hash){
	my $header = $key;
	$header =~ s/\n//g;
	my $seq = $seq_hash{$key};

	my $seqID;
	if($seq_type eq "mRNA" || $seq_type eq "prot_CDS"){
	    if($header =~ m/\>([\S]+)/){
		$seqID = $1;
	    }
	}
	elsif($seq_type eq "nuc_CDS"){
	    if($header =~ m/\[protein\_id\=([^\]]+)/){
		$seqID = $1;
	    }
	}   
	if($seqID ~~ @ID_array){
	    #print "SEQID: $seqID \n";
	    unless($header =~ m/\Q$ncbi_flag\E/){
		$header.=$ncbi_flag;
	    }
	    print FILEOUT $header."\n".$seq."\n";
	    push(@headers, $header);
	}
    }
    close FILEOUT;
    return @headers;
}

sub get_contig_lengths{
    my $genome_file = $_[0];
    open(GENOME_IN, "$genome_file");
    my %contig_length_hash;
    my $contig;
    while(<GENOME_IN>){
	my $line = $_;
	if($line =~ m/^>([\S]+)/){
	    $contig = $1;
	    chomp $contig;
	}
	elsif($line =~ m/[\S]+/){
	    if($contig_length_hash{$contig}){
		my $contig_length = $contig_length_hash{$contig};
		$line =~ s/\n//g;
		my $add_length = length($line);
		my $new_length = $contig_length + $add_length;
		$contig_length_hash{$contig} = $new_length;
	    }
	    else{
		$line =~ s/\n//g;
		my $contig_length = length($line);
		$contig_length_hash{$contig} = $contig_length;
	    }
	}
    }
    close GENOME_IN;
    return \%contig_length_hash;
}
	

