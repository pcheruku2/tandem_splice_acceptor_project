#!/usr/bin/perl

use strict;

my $trimer_file = $ARGV[0];
my $fasta_file = $ARGV[1];
my $chrom      = $ARGV[2];

my $frequency_file = "counts_$chrom.out";
my $trimer_sites = "trimer_sites_$chrom.out";


open C, "> $frequency_file" or die "Cannot open $frequency_file : $!";
open T, "> $trimer_sites" or die "Cannot open $trimer_sites : $!";


my %trimer_hash =();
my %trimer_count=();
my $sequence_count=();

my $c=0;
my $chr=();
my $from_base1=();
my $to_base1=();

open F, "< $trimer_file" or die "Cannot open $trimer_file :$!";
while(<F>){
	chomp;
	my(@l)=split;
	$trimer_hash{$l[0]}=1;
}
close F;



open F, "< $fasta_file" or die "Cannot open $fasta_file : $!";
while(<F>){
	chomp;
	$c++;
	my $mod=$c%2;
	if($_=~/^>/ && $mod==1){
		$_=~/^>(.*):(\d+)-(\d+)/;		
		$chr=$1;
		$from_base1=$2;
		$to_base1=$3;
		$sequence_count++;
		
	}elsif($mod==0){
		my $seq=$_;

		for(my $a=0;$a<=length($seq)-3; $a++){
                      	my $trimer=substr($seq, $a, 3);
			my $coord=$from_base1+$a;
			if($trimer_hash{$trimer}){
				print T "$chr\t$coord\t$trimer\n";
				$trimer_count{$trimer}++;
			}
                }

		$chr=();
		$from_base1=();
		$to_base1=();
	}		
}
close F;


foreach my $k(sort keys %trimer_hash){
	if($trimer_count{$k}){
		print C "$k\t$trimer_count{$k}\n";
	}else{
		print C "$k\t0\n";
	}
}

close T;
close C;
