#!/usr/bin/perl

use strict;

my $fa_file = $ARGV[0];

my %ag_count_pos=();
my %ag_count_seq=();

my $total_sequence_count=0;
my $canonical_AG=0;
my $non_canonical=0;

open F, "< $fa_file" or die "Cannot open $fa_file : $!";
while(<F>){
	chomp;
	if($_=~/^>/){ $total_sequence_count++; next; }

	my($seq)=$_;
	$seq=uc($seq);
	my $len=length($seq);
	my (@nuc)=split(//,$seq);
	my $AG_count=0;
	my $flag=0;

	#print "$seq\t$len\n";

	for(my $i=@nuc-1;$i>0;$i--){
	 my $pos=$i-@nuc+1;

		if(($pos==-1) && ($nuc[$i] eq "G") && ($nuc[$i-1] eq "A")){ $canonical_AG++; $flag=1;  }elsif(($pos==-1) && ($nuc[$i] ne "G") && ($nuc[$i-1] ne "A")) { $non_canonical++;} 
	 	if(($pos<0) && ($nuc[$i] eq "G") && ($nuc[$i-1] eq "A") && $flag==1){
			if($pos==-2){ 
				print "$seq\t$len\n";	
	 			print "$nuc[$i-1]$nuc[$i]\t$pos\n";
			}
			$ag_count_pos{$pos}++;
			$AG_count++;
			if($pos<-1){ last; }
	 	}

	}
	#print "AG count : $AG_count\n";
	$ag_count_seq{$AG_count}++;
}

print "Canonical AG sequences = $canonical_AG\n";
print "Non-canonical sequences = $non_canonical\n";
print "Total sequences = $total_sequence_count\n";

print "\n\n";

foreach my $k(sort {$a<=>$b} keys %ag_count_seq){
	#print "$k\t$ag_count_seq{$k}\n";
} 

print "\n\n";

foreach my $k(sort {$b<=>$a} keys %ag_count_pos){
	print "$k\t$ag_count_pos{$k}\n";
} 


