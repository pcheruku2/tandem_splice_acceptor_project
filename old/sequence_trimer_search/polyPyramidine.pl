#!/usr/bin/perl

use strict;

my $file = $ARGV[0];

open F, "< $file " or die "Cannot open $file : $!";
while(<F>){
	chomp;
	if($_=~/^>/){
		print "$_\n";
	}else{
		my $polyString=$_;
		my $nuclString=$_;
		$polyString=~s/A/R/g;
		$polyString=~s/G/R/g;
		$polyString=~s/C/Y/g;
		$polyString=~s/T/Y/g;

		#print "$nuclString\n";
		print "$polyString\n";
	}
}
close F;
