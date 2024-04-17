#!/usr/bin/perl

use strict;

my $fasta_file = $ARGV[0];

open F, "< $fasta_file" or die "Cannot open $fasta_file : $!";
while(<F>){

		  chomp;
		  if($_=~/^>/){
					 print "$_\n";
		  }else{
					 my $positive_strand_seq = $_;
					 my $reverse_complement_seq = reverse $positive_strand_seq;
					 $reverse_complement_seq =~ tr/ATGCatgc/TACGtacg/;
					 print "$reverse_complement_seq\n";
		 }
}
close F;

				
