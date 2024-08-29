#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::SeqIO;
use File::Basename;

# Define the reference sequence
my $reference = 'CTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCGGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGGAAAAAACAGCGATTACAATTGGGACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCA';
my $reference_aa = 'GASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARGKNSDYNWDFQHWGQGTLVTVSS';

# Define the patterns to search for
my $pattern1 = 'tctggggctgaggtgaagatgc';
my @pattern2_options = (
    'GAGAGTCAGTCCTTCCCAAATGTCTTCC',
    'GCCAAAACGACACCCCCATCTG',
    'GCCAAAACAACAGCCCCATCGG',
    'GCCAAAACAACACCCCCATCAGT',
    'GCCAAAACAACAGCCCCATCGG',
    'GCTACAACAACAGCCCCATCTG',
);

# Prompt user for the input folder path containing FASTA files
print "Enter the path to the folder containing FASTA files: ";
my $input_folder_path = <STDIN>;
chomp $input_folder_path;

# Create a folder to store the output files
print "\nEnter path to the folder for output files: ";
my $output_folder = <STDIN>;
chomp $output_folder;

# Open the input directory
opendir(my $dir_fh, $input_folder_path) or die "Could not open directory $input_folder_path: $!";

while (my $file = readdir($dir_fh)) {
	
	next unless $file =~ /\.fasta$/; # Process only files with .fasta extension

	my $input_file_path = "$input_folder_path/$file";
	
	# Read sequences from the FASTA file
	my $seqio = Bio::SeqIO->new(-file => $input_file_path, -format => 'Fasta');

	# Initialize variables for counting
	my $total_sequences_analyzed = 0;
	my $total_pattern_matches = 0;
	my $total_sequences_ignored = 0;

	# Initialize an array to store mismatch counts at each position
	my @mismatch_counts = (0) x length($reference);
	my @aa_mismatch_counts = (0) x length($reference_aa);

	# Initialize an array to store amino acid sequences
	my @proteins;

	# Iterate through each sequence in the input file
	while (my $seq = $seqio->next_seq) {
		$total_sequences_analyzed++;

		my $sequence = $seq->seq;
		my $header = $seq->display_id;

		# Extract the relevant substring based on the patterns
		if ($sequence =~ /$pattern1(.{323})($pattern2_options[0]|$pattern2_options[1]|$pattern2_options[2]|$pattern2_options[3]|$pattern2_options[4]|$pattern2_options[5])/i) {
			$total_pattern_matches++;

			my $relevant_sequence = $1;  # Extract the 323 base pairs between the patterns

			# Ensure the relevant sequence length matches the reference length
			unless (length($relevant_sequence) == length($reference)) {
				print "Skipped: Length mismatch for sequence:\n$relevant_sequence\n";
				$total_sequences_ignored++;
				next;
			}

			# Count the number of mismatches within the sequence
			my $mismatches_within_sequence = 0;
			for my $position (0 .. length($reference) - 1) {
				my $ref_base = uc substr($reference, $position, 1);
				my $seq_base = uc substr($relevant_sequence, $position, 1);
				$mismatches_within_sequence++ if $seq_base ne $ref_base;
				$mismatch_counts[$position]++ if $seq_base ne $ref_base;
			}

			# Translate the sequence in reading frame 3 and add the sequence to the protein list
			my $relevant_sequence_rf3 = substr($relevant_sequence, 2);
			my $seq_obj = Bio::Seq->new(-seq => $relevant_sequence_rf3, -alphabet => 'dna');
			my $relevant_aa_sequence = $seq_obj->translate->seq;
		
			# Count the number of amino acid substitutions
			my $mismatches_within_aa_sequence = 0;
			my $aa_description = "GL";
			for my $position_aa (0 .. length($reference_aa) - 1) {
				my $ref_aa = uc substr($reference_aa, $position_aa, 1);
				my $seq_aa = uc substr($relevant_aa_sequence, $position_aa, 1);
				$mismatches_within_aa_sequence++ if $seq_aa ne $ref_aa;
				$aa_mismatch_counts[$position_aa]++ if $seq_aa ne $ref_aa;
				my $position_aa_abbott = $position_aa + 15;
				if ($seq_aa ne $ref_aa) {
					if ($aa_description eq "GL") {
						$aa_description = "$ref_aa$position_aa_abbott$seq_aa";
					} else {
						$aa_description = "$aa_description,$ref_aa$position_aa_abbott$seq_aa";
					}
				}
			}
		
			push(@proteins , ("$relevant_aa_sequence\t$mismatches_within_aa_sequence\t$aa_description"));

			# Print the relevant sequence, its length, the number of mismatches, and the header for debugging
			print "Header: $header\n";
			print "Relevant Sequence:\n$relevant_sequence\n";
			print "Length: " . length($relevant_sequence) . "\n";
			print "Mismatches Within Sequence: $mismatches_within_sequence\n";
			print "Amino acid subsitutions: $mismatches_within_aa_sequence\n";
		}
		else {
			print "Skipped: Sequence does not match the pattern:\n$sequence\n";
			$total_sequences_ignored++;
		}
	}

	# Count unique amino acid sequences and print the list with frequencies
	my %clones;
	my @unique = grep !$clones{$_}++, @proteins;
	foreach my $key (sort {$clones{$b} <=> $clones{$a}} keys %clones) {
		print "$clones{$key}\t$key\n";
	}
	my $protein_count = keys %clones;
	print "\n\nTotal number of unique amino acid sequences: $protein_count\n\n";

	# Print the results
	print "Position\tMismatch Counts\tMismatchFrequency\n";
	for my $position (0 .. $#mismatch_counts) {
		my $positionplusone = $position + 1;
		print "$positionplusone\t$mismatch_counts[$position]\t";
		my $mismatch_frequency = $mismatch_counts[$position] / $total_pattern_matches;
		print "$mismatch_frequency\n";
	}

	print "Position\tSubstitutions\tSubstitutionFrequency\n";
	for my $position_aa (0.. $#aa_mismatch_counts) {
		my $position_aaplusone = $position_aa + 1;
		print "$position_aaplusone\t$aa_mismatch_counts[$position_aa]\t";
		my $substitution_frequency = $aa_mismatch_counts[$position_aa] / $total_pattern_matches;
		print "$substitution_frequency\n";
	}

	# Print the summary
	print "\nSummary:\n";
	print "Total Sequences Analyzed: $total_sequences_analyzed\n";
	print "Total Pattern Matches: $total_pattern_matches\n";
	print "Total Sequences Ignored: $total_sequences_ignored\n";

	### Print the results to an output file

	# Modify the output file path based on the input file name
    my $output_file_name = fileparse($file, qr/\.[^.]*/);
    my $output_file_path = "$output_folder/$output_file_name.txt";


	# Open the file for writing
	open my $output_fh, '>', $output_file_path or die "Could not open file '$output_file_path' for writing: $!\n";

	# Print a header line

	my $short_file_name = fileparse($input_file_path);
	print $output_fh "$short_file_name\n\n";

	# Print the results to the text file

	print $output_fh "Nucleotides\nPosition_Abbott\tMismatch Counts\tMismatchFrequency\n";
	for my $position (0 .. $#mismatch_counts) {
		my $positionabbott = $position + 46;
		print $output_fh "$positionabbott\t$mismatch_counts[$position]\t";
		my $mismatch_frequency = $mismatch_counts[$position] / $total_pattern_matches;
		print $output_fh "$mismatch_frequency\n";
	}


	print $output_fh "\nAmino Acids\nPosition_Abbott\tSubstitutions\tSubstitutionFrequency\n";
	for my $position_aa (0.. $#aa_mismatch_counts) {
		my $position_aa_abbott = $position_aa + 15;
		print $output_fh "$position_aa_abbott\t$aa_mismatch_counts[$position_aa]\t";
		my $substitution_frequency = $aa_mismatch_counts[$position_aa] / $total_pattern_matches;
		print $output_fh "$substitution_frequency\n";
	}

	# Print the summary to the text file
	print $output_fh "\nSummary:\n";
	print $output_fh "Total Sequences Analyzed: $total_sequences_analyzed\n";
	print $output_fh "Total Pattern Matches: $total_pattern_matches\n";
	print $output_fh "Total Sequences Ignored: $total_sequences_ignored\n";

	#print unique amino acid sequences, counts, and distance from reference
	print $output_fh "\n\nUnique Amino Acid Sequences:\nCount\tSequence\tDistanceFromGL\n";
	foreach my $key (sort {$clones{$b} <=> $clones{$a}} keys %clones) {
		print $output_fh "$clones{$key}\t$key\n";
	}

	# Close the file
	close $output_fh;

	# Inform the user that the results have been saved
	print "Results have been saved to: $output_file_path\n";
}

closedir($dir_fh);

# Script for analysis of SHM and amino acid substitutions to the germline VRC01 BCR heavy chain sequence
# John M. Lindner (john.m.lindner@gmail.com)
# This script generated with the aid of ChatGPT 3.5
# 2024_01_10_v2 Added translation to amino acid sequence and comparison of each translation to the reference
# 2024_01_11_v3 Added frequency of amino acid, counting, and sorting. "Clones" are reported with frequency and distance from reference
# 2024_01_16_v4 Added naming of clonotypes based on amino acid substitutions from germline
# This script counts mismatches between the reference VRC01 heavy chain sequence
# and reports mutation frequencies from NGS data in FASTA format
# In v3, this sequence is specified by the script and is a flexible protein linker
# Note the following CDR sequences and coordinates for reference (NOT from IgBLAST):
# CDR1: GYTFTGYYMH , positions 12-21 of the translated reference (add 14 for Abbott)
# CDR2: INPNSGGTNYAQKFQG , positions 37-52 of the translated reference (add 14 for Abbott)
# CDR3: ARGKNSDYNWDFQH , positions 83-96 of the translated reference (add 14 for Abbott)

