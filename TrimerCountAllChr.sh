#!/bin/sh

homeDir='/home/frostfg/NAGNAG/genome_wide_trimer_code_V2'
fastaDir='/home/frostfg/NAGNAG/ChrFastas'
refDir='/home/frostfg/NAGNAG/'

cd $refDir
AllChroms=`\ls UCSCfastas/*fa | sed -e "s/UCSCfastas\\///g" | sed -e "s/\\.fa//g"`

module load samtools
cd $homeDir
for chr in $AllChroms;
do
	echo Working on "$chr" exonic windows...
	# Make region.fa file for this chromosome, also convert from 0- to 1-base
	# Positive strand gene regions
	cat "$refDir"/GRCh37_MANE_Transcripts_Splice_Acceptor_Exonic_Windows.bed | awk -v chr="$chr" '$1 == chr && $6=="+" {print $1":"$2+1"-"$3}' > "$homeDir"/OutputRegions/"$chr"_+_Exonic_regions.txt
	
	# Negative strand gene regions
        cat "$refDir"/GRCh37_MANE_Transcripts_Splice_Acceptor_Exonic_Windows.bed | awk -v chr="$chr" '$1 == chr && $6=="-" {print $1":"$2+1"-"$3}' > "$homeDir"/OutputRegions/"$chr"_-_Exonic_regions.txt

	# Not entirely sure what's happening here, but shit is being indexed, so that's good
	# Positive strand gene regions first
	samtools faidx -n 100 "$refDir"UCSCfastas/"$chr".fa -r "$homeDir"/OutputRegions/"$chr"_+_Exonic_regions.txt > "$fastaDir"/"$chr"_+_Exonic_regions.fa
	# Now negative strand gene regions
        samtools faidx -n 100 -i "$refDir"UCSCfastas/"$chr".fa -r "$homeDir"/OutputRegions/"$chr"_-_Exonic_regions.txt > "$fastaDir"/"$chr"_-_Exonic_regions.fa

	# Slap both strands into one .fa file
	cat "$fastaDir"/"$chr"_+_Exonic_regions.fa "$fastaDir"/"$chr"_-_Exonic_regions.fa >  "$fastaDir"/"$chr"_all_Exonic_regions.fa
	rm "$fastaDir"/"$chr"_+_Exonic_regions.fa
	rm "$fastaDir"/"$chr"_-_Exonic_regions.fa

	# On to the trimer count script
	cd PerlTrimerOut/ExonicOut
	perl "$homeDir"/sequence_trimer_search/trimer_count.pl "$homeDir"/sequence_trimer_search/uniq_trimers_AG_gain_loss_neutral.txt "$fastaDir"/"$chr"_all_Exonic_regions.fa ExonicOut_"$chr"
	cd $homeDir
	echo Working on "$chr" exonic windows...... Done
done

cd $homeDir
for chr in $AllChroms;
do
        echo Working on "$chr" intronic windows...
        # Make region.fa file for this chromosome, also convert from 0- to 1-base
        # Positive strand gene regions
        cat "$refDir"/GRCh37_MANE_Transcripts_Splice_Acceptor_Intronic_Windows.bed | awk -v chr="$chr" '$1 == chr && $6=="+" {print $1":"$2+1"-"$3}' > "$homeDir"/OutputRegions/"$chr"_+_Intronic_regions.txt

        # Negative strand gene regions
        cat "$refDir"/GRCh37_MANE_Transcripts_Splice_Acceptor_Intronic_Windows.bed | awk -v chr="$chr" '$1 == chr && $6=="-" {print $1":"$2+1"-"$3}' > "$homeDir"/OutputRegions/"$chr"_-_Intronic_regions.txt

        # Not entirely sure what's happening here, but shit is being indexed, so that's good
        # Positive strand gene regions first
        samtools faidx -n 100 "$refDir"UCSCfastas/"$chr".fa -r "$homeDir"/OutputRegions/"$chr"_+_Intronic_regions.txt > "$fastaDir"/"$chr"_+_Intronic_regions.fa
        # Now negative strand gene regions
        samtools faidx -n 100 -i "$refDir"UCSCfastas/"$chr".fa -r "$homeDir"/OutputRegions/"$chr"_-_Intronic_regions.txt > "$fastaDir"/"$chr"_-_Intronic_regions.fa

        # Slap both strands into one .fa file
        cat "$fastaDir"/"$chr"_+_Intronic_regions.fa "$fastaDir"/"$chr"_-_Intronic_regions.fa >  "$fastaDir"/"$chr"_all_Intronic_regions.fa
        rm "$fastaDir"/"$chr"_+_Intronic_regions.fa
        rm "$fastaDir"/"$chr"_-_Intronic_regions.fa

        # On to the trimer count script
        cd PerlTrimerOut/IntronicOut
        perl "$homeDir"/sequence_trimer_search/trimer_count.pl "$homeDir"/sequence_trimer_search/uniq_trimers_AG_gain_loss_neutral.txt "$fastaDir"/"$chr"_all_Intronic_regions.fa IntronicOut_"$chr"
        cd $homeDir
        echo Working on "$chr" intronic windows...... Done
done
