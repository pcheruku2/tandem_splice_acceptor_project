#!/usr/bin/perl 

use strict; 

my $trimer_hash_file = $ARGV[0];
my $trimer_bed_file  = $ARGV[1];
my $trimer_vcf_file  = $ARGV[2];


my %trimer_lookup_hash =();
my %trimer_lookup_hash_rc =();
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

		# reverse complement

		my $trc = reverse $t;
		my $arc = reverse $a;
	
		$trc =~ tr/ATGCatgc/TACGtacg/; 
		$arc =~ tr/ATGCatgc/TACGtacg/; 

		$trimer_lookup_hash_rc{$trc}{$arc}=$b;
		
	}
}
close F;

# Trimer bed file 
# 21	30402916	30402917	AAC:2:+

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
		my ($trimer,$local_coord,$strand)=split(/:/,$trimer_data);

		if($strand eq "+"){
			my($a,$b,$c)=split(//,$trimer); 

			# make sure the second base is checked against reference allele
			# $l[3] = ref allele ; $l[4] = alt allele

			if(($b eq $l[3]) &&  $trimer_lookup_hash{$trimer}{$l[4]}){
				print "$trimer\t$b\t$l[4]\t$strand\t$trimer\t$b\t$l[4]\t$trimer_lookup_hash{$trimer}{$l[4]}\t$local_coord\t$_\n";
			}
		}elsif($strand eq "-"){
			my $trimer_rc = reverse $trimer;
			$trimer_rc =~ tr/ATGCatgc/TACGtacg/; 

			my($a,$b,$c)=split(//,$trimer_rc); 

			my $b_rc=$b;
			$b_rc=~ tr/ATGCatgc/TACGtacg/;
			my $alt_rc=$l[4];
			$alt_rc=~ tr/ATGCatgc/TACGtacg/;
			

			if(($b eq $l[3]) &&  $trimer_lookup_hash_rc{$trimer_rc}{$l[4]}){
				print "$trimer\t$b_rc\t$alt_rc\t$strand\t$trimer_rc\t$b\t$l[4]\t$trimer_lookup_hash_rc{$trimer_rc}{$l[4]}\t$local_coord\t$_\n";
			}
		}
	}
}
close F; 

	

