#!/usr/bin/perl

use strict;

my $size_file = $ARGV[0];
my $window_size = $ARGV[1];
my $number_of_seq = $ARGV[2];

my $outfile_name = "random_w$window_size"."_n$number_of_seq";


open B, "> $outfile_name.bed" or die "Cannot open $outfile_name.bed : $!";
open S, "> $outfile_name.sam.txt" or die "Cannot open $outfile_name.sam.txt : $!";

my @array=();
my $c=0;

open F, "< $size_file " or die "Cannot open $size_file : $!";
while(<F>){
	chomp;
	my(@l)=split;
	$array[$c]="$l[0]:$l[1]";
	$c++;
}
close F;



for (my $i=0; $i<$number_of_seq; $i++){
	my $chr_random = int(rand(@array));
	#print "$chr_random\t$array[$chr_random]\n";
	my($chr,$len)=split(/:/,$array[$chr_random]);
 	$len=$len-$window_size;
	my $random_position = int( rand( $len ) ) + 1;
	my $bed_from_0=$random_position;
	my $bed_from_1=$random_position+1;
	my $bed_to=$random_position+$window_size;
	print B "$chr\t$bed_from_0\t$bed_to\trandom_$window_size\n";
	print S "$chr:$bed_from_1-$bed_to\n";
}

close B;
close S;

