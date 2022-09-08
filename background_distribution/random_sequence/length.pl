#!/usr/bin/perl

use strict;

my $size_file = $ARGV[0];
my $chrom_file = $ARGV[1];

my %size_hash=();
my $chr=();
my $total_length=0;

open F, "< $size_file" or die "Cannot open $size_file: $!";
while(<F>){
	chomp;
	my(@l)=split;
	$size_hash{$l[0]}=$l[1];
}
close F;


open(IN, "gunzip -c $chrom_file |") or die "Cannot gunzip $chrom_file: $!";
while(<IN>){
	chomp;
	if($_=~/>/){
		$_=~/>(.*)/;
		$chr=$1;		
	}else{
		my $l=length($_);
		$total_length+=$l;
	}
}
close IN;


print "$chr\t";
if($size_hash{$chr}){
	print "$size_hash{$chr}\t";
}else{
	print "NA\t";
}
print "$total_length\n";




