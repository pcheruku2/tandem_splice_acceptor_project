#!/usr/bin/perl

#############################################################################################################
# script name : trimer_count.pl
#
# purpose : Given a trimers of interest look up file and fasta sequences 
# of fixed length (31 bases in exon, 33 bases in intron; with one base overlap at junction),
# the goal of the script is to search the primary input fasta sequence for trimers of interest
# and return a bedfile with genomic coordinates, the trimer, local coordinates relative to the junction, 
# and strand infromation; it also returns a file with a summary of counts of trimers found. 
#############################################################################################################



use strict;

#########
# Input #
#########

my $trimer_file = $ARGV[0];
my $fasta_file = $ARGV[1];
my $suffix      = $ARGV[2];


# Output file names
my $frequency_file = "counts_$suffix.out";
my $trimer_sites_bedfile = "trimer_sites_$suffix.bed";


# Open output files 

open C, "> $frequency_file" or die "Cannot open $frequency_file : $!";
open B, "> $trimer_sites_bedfile" or die "Cannot open $trimer_sites_bedfile : $!";

# Hashes & variables 

my %trimer_hash =();
my %trimer_count=();
my $sequence_count=();

my $c=0;
my $chr=();
my $cn=();
my $from_base1=();
my $to_base1=();
my $seq_class=();
my $strand=();


#######################################################
# Read in the trimer look up table
# Table format: 
# trimer and consequence of variant at 2nd base
# AAA	G:gain 
#######################################################

open F, "< $trimer_file" or die "Cannot open $trimer_file :$!";
while(<F>){
	chomp;
	my(@l)=split;
	$trimer_hash{$l[0]}=1;
}
close F;


#######################################################################################
# Read in fasta sequence : 
#    intronic sequences = 33 bases; 1 base into exon at junction
#    exonic sequences = 31 bases; 1 base into intron at junction
#######################################################################################

open F, "< $fasta_file" or die "Cannot open $fasta_file : $!";
while(<F>){
	chomp;
	$c++;
	my $mod=$c%2;
	if($_=~/^>/ && $mod==1){
		if($_=~/rc/){ 			
			$_=~/^>(.*):(\d+)-(\d+)\/rc/;		
			$chr=$1;
			$from_base1=$2;
			$to_base1=$3;
			$strand="-";
		}else{ 
			$_=~/^>(.*):(\d+)-(\d+)/;		
			$chr=$1;
			$from_base1=$2;
			$to_base1=$3;
			$strand ="+";
		}
		$sequence_count++;
		
	}elsif($mod==0){
		my $seq=$_;
		my $len = length($seq);
		my $local_coordn = ();

		if($len == 31){ 
			$seq_class="exon";
		}elsif($len == 33){
			$seq_class="intron";
		}
			
		for(my $a=0;$a<=length($seq)-3; $a++){
                      	my $trimer=substr($seq, $a, 3);
			my $coord=();
			my $local_coordinate="N";

			if($seq_class eq "exon"  && $strand eq "+"){
				$coord=$from_base1+$a; # zero-base for bedfile;
				$local_coordinate=$a+1;

				#print "seq class : exon\n";
				#print "seq : $seq\n";
				#print "strand : $strand\n";
				#print "trimer : $trimer\n";
				#print "a   : $a\n";
				#print "from : $from_base1\n";
				#print "to   : $to_base1\n";
				#print "coordinate : $coord\n";
				#print "local coordinate : $local_coordinate\n\n";

				if($trimer_hash{$trimer}){
					my $lcd=$coord+1;
					my $cn=$chr;
					$cn=~s/chr//g;

					# Write to bedfile
					print B "$cn\t$coord\t$lcd\t$trimer:$local_coordinate:$strand\n";
					$trimer_count{$trimer}++;
				}
	
			}elsif($seq_class eq "intron" && $strand eq "+"){
				$coord=$from_base1+$a;
				$local_coordinate=$a-$len+2;

				#print "seq class : intron\n";
				#print "seq : $seq\n";
				#print "strand : $strand\n";
				#print "trimer : $trimer\n";
				#print "a   : $a\n";
				#print "from : $from_base1\n";
				#print "to   : $to_base1\n";
				#print "coordinate : $coord\n";
				#print "local coordinate : $local_coordinate\n\n";

				if($trimer_hash{$trimer}){
					my $lcd=$coord+1;
					my $cn=$chr;
					$cn=~s/chr//g;

					# Write to bedfile
					print B "$cn\t$coord\t$lcd\t$trimer:$local_coordinate:$strand\n";
					$trimer_count{$trimer}++;
				}

			}elsif($seq_class eq "exon"  && $strand eq "-"){
				$coord=$to_base1-$a-2; # zero-base for bedfile;
				$local_coordinate=$a+1;

				#print "seq class : exon\n";
				#print "seq : $seq\n";
				#print "strand : $strand\n";
				#print "trimer : $trimer\n";
				#print "a   : $a\n";
				#print "from : $from_base1\n";
				#print "to   : $to_base1\n";
				#print "coordinate : $coord\n";
				#print "local coordinate : $local_coordinate\n\n";

				if($trimer_hash{$trimer}){
					my $lcd=$coord+1;
					my $cn=$chr;
					$cn=~s/chr//g;

					# Write to bedfile
					print B "$cn\t$coord\t$lcd\t$trimer:$local_coordinate:$strand\n";
					$trimer_count{$trimer}++;
				}
			}elsif($seq_class eq "intron" && $strand eq "-"){
				$coord=$to_base1-$a-2;
				$local_coordinate=$a-$len+2;

				#print "seq class : intron\n";
				#print "seq : $seq\n";
				#print "strand : $strand\n";
				#print "trimer : $trimer\n";
				#print "a   : $a\n";
				#print "from : $from_base1\n";
				#print "to   : $to_base1\n";
				#print "coordinate : $coord\n";
				#print "local coordinate : $local_coordinate\n\n";

				if($trimer_hash{$trimer}){
					my $lcd=$coord+1;
					my $cn=$chr;
					$cn=~s/chr//g;

					# Write to bedfile
					print B "$cn\t$coord\t$lcd\t$trimer:$local_coordinate:$strand\n";
					$trimer_count{$trimer}++;
				}
			}

	

                }

		$chr=();
		$from_base1=();
		$to_base1=();
		$strand=();
	}		
}
close F;
########################################################
# Write out basic trimer count statistics to output file
########################################################

foreach my $k(sort keys %trimer_hash){
	if($trimer_count{$k}){
		print C "$k\t$trimer_count{$k}\n";
	}else{
		print C "$k\t0\n";
	}
}

close C;
close B;
