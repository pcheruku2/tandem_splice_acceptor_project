#!/bin/sh

homeDir='/home/frostfg/NAGNAG/genome_wide_trimer_code_V2'
refDir='/home/frostfg/NAGNAG/'

cd $refDir
AllChroms=`\ls UCSCfastas/*fa | sed -e "s/UCSCfastas\\///g" | sed -e "s/\\.fa//g"`

cd $homeDir
for chr in $AllChroms;
do
	echo Working on "$chr"...
	
	# Make region.fa file for this chromosome, also convert from 0- to 1-base
	cat "$homeDir"/PerlTrimerOut/ExonicOut/trimer_sites_ExonicOut_"$chr".out | awk '{print $1 "\t" $2 "\t"  $2+1 "\t" $3}' | sed -e 's/chr//g' > "$homeDir"/tmp_ExonicOut_"$chr".txt
        cat "$homeDir"/PerlTrimerOut/IntronicOut/trimer_sites_IntronicOut_"$chr".out | awk '{print $1 "\t" $2 "\t"  $2+1 "\t" $3}' | sed -e 's/chr//g' > "$homeDir"/tmp_IntronicOut_"$chr".txt

	echo Working on "$chr"... Done
done

cat tmp_ExonicOut_*txt > AllChr_NAGNAG_Exonic_Windows.bed
rm tmp_ExonicOut_*txt

cat tmp_IntronicOut_*txt > AllChr_NAGNAG_Intronic_Windows.bed
rm tmp_IntronicOut_*txt


