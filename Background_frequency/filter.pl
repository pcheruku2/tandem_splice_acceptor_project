#!/usr/bin/perl

use strict;

# filter out any sequences with "N" bases 
# convert all lower case bases to upper case

my $fa_file = $ARGV[0];

my $c=0;
my $h=();
open F, "< $fa_file" or die "Cannot open $fa_file : $!";
while(<F>){
	chomp;
	$c++;
	my $m = $c%2;
	if($m == 1){
		$h=$_;
	}elsif($m == 0){
		my $seq=$_;
		$seq=uc($seq);
		if($seq=~/N/){ next; }
		else{
			print "$h\n$seq\n";
		}
	}
}
close F;

		
