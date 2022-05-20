#!/usr/bin/perl 

use strict; 

my $trimer_hash_file = $ARGV[0];
my $trimer_bed_file  = $ARGV[1];
my $trimer_vcf_file  = $ARGV[2];


my %trimer_lookup_hash =();
my %trimer_loci_hash   =();

# Lookup table - unique trimers
# AAG	C:loss;G:neutral;T:loss
open F, "< $trimer_hash_file" or die "Cannot open $trimer_hash_file : $!";
while(<F>){
	chomp;
	my(@l)=split;
	my $t = $l[0];
	my (@c)=split(/;/, $l[1]);
	foreach my $k(@c){
		my($a,$b)=split(/:/, $k);
		$trimer_lookup_hash{$t}{$a}=$b;
	}
}
close F;

# Trimer bed file 
# 19	1218398	1218399	CTG
open F, "< $trimer_bed_file" or die "Cannot open $trimer_bed_file : $!";
while(<F>){
	chomp;
	my(@l)=split;
	$trimer_loci_hash{$l[0]}{$l[2]}=$l[3];
}
close F;


# VCF file subset for trimer loci

open F, "< $trimer_vcf_file " or die "Cannot open $trimer_vcf_file : $!";
while(<F>){
	chomp;
	if($_=~/^#/){ next; }
	my (@l)=split;
	if($trimer_loci_hash{$l[0]}{$l[1]}){
		my $trimer_data=$trimer_loci_hash{$l[0]}{$l[1]};
		my ($trimer,$local_coord)=split(/:/,$trimer_data);
		
		my($a,$b,$c)=split(//,$trimer_loci_hash{$l[0]}{$l[1]});
		if(($a eq $l[3]) &&  $trimer_lookup_hash{$trimer}{$l[4]}){
			print "$trimer\t$a\t$l[4]\t$trimer_lookup_hash{$trimer}{$l[4]}\t$local_coord\t$_\n";
		}
	}
}
close F; 

	

