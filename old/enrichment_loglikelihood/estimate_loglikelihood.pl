#!/bin/usr/perl

use strict;

my $nucl_prob_file = $ARGV[0];
my $trimer_count_file = $ARGV[1];

my %nucl_prob=();
my %trimer_count=();
my $trimer_sum=0;
my %trimer_prob=();
my %trimer_log_likelihood=();

#A	0.295269605
open F, "< $nucl_prob_file" or die "Cannot open $nucl_prob_file : $!";
while(<F>){
	chomp;
	my(@l)=split;
	$nucl_prob{$l[0]}=$l[1];
}
close F;


#405 TCG
open F, "< $trimer_count_file" or die "Cannot open trimer_count_file : $!";
while(<F>){
	chomp;
	my(@l)=split;
	$trimer_count{$l[1]}=$l[0];
	$trimer_sum+=$l[0];
}
close F;


foreach my $t(sort keys %trimer_count){
	$trimer_prob{$t}=$trimer_count{$t}/$trimer_sum;
	my($a,$b,$c)=split(//,$t);
	$trimer_log_likelihood{$t}=log(($trimer_prob{$t})/($nucl_prob{$a}*$nucl_prob{$b}*$nucl_prob{$c}));
	print "$t\t$trimer_count{$t}\t$trimer_prob{$t}\t$trimer_log_likelihood{$t}\n";
}

print "\ntrimer sum : $trimer_sum\n";


	
