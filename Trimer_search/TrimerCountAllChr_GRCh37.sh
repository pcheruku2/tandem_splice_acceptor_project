#!/bin/sh

homeDir='/path/to/analysis directory'
fastaDir='/path/to/splice/region/fastAs'
refDir='/path/to/reference/chromosome/fastAs'


cd $refDir
AllChroms=`\ls UCSCfastas/*fa | sed -e "s/UCSCfastas\\///g" | sed -e "s/\\.fa//g"`

module load samtools

# Ensure output directory structure exists:
mkdir -p ${homeDir}/OutputRegions

for chr in $AllChroms;
do
	cd ${homeDir}
	echo Working on "$chr" exonic windows...
	# Make region.fa file for this chromosome, also convert from 0- to 1-base
	# Positive strand gene regions
	cat ${refDir}/GRCh37_MANE_Transcripts_Splice_Acceptor_Exonic_Windows.bed | \
		awk -v chr=${chr} '$1 == chr && $6=="+" {print $1":"$2+1"-"$3}' >
			${homeDir}/OutputRegions/${chr}_+_Exonic_regions.txt
	
	# Negative strand gene regions
        cat ${refDir}/GRCh37_MANE_Transcripts_Splice_Acceptor_Exonic_Windows.bed | \
		awk -v chr=${chr} '$1 == chr && $6=="-" {print $1":"$2+1"-"$3}' > \
			${homeDir}/OutputRegions/${chr}_-_Exonic_regions.txt

	# Index chromosomal fastAs
	# Positive strand gene regions first
	samtools faidx -n 100 \
		${refDir}/${chr}.fa \
		-r ${homeDir}/OutputRegions/${chr}_+_Exonic_regions.txt > \
			${fastaDir}/${chr}_+_Exonic_regions.fa

	# Now negative strand gene regions
        samtools faidx -n 100 \
		-i \
		${refDir}/${chr}.fa \
		-r ${homeDir}/OutputRegions/${chr}_-_Exonic_regions.txt > \
			${fastaDir}/${chr}_-_Exonic_regions.fa

	# On to the trimer count script
        cd PerlTrimerOut/ExonicOut
	perl ${homeDir}/sequence_trimer_search/trimer_count.pl \
		${homeDir}/sequence_trimer_search/uniq_trimers_AG_gain_loss_neutral.txt \
		${fastaDir}/${chr}_+_Exonic_regions.fa \
		 ${chr}_+_Exonic_region

	perl  ${homeDir}/sequence_trimer_search/trimer_count.pl \
                 ${homeDir}/sequence_trimer_search/uniq_trimers_AG_gain_loss_neutral.txt \
                 ${fastaDir}/${chr}_-_Exonic_regions.fa \
                 ${chr}_-_Exonic_region
	
	echo Working on "$chr" exonic windows...... Done
done

for chr in $AllChroms;
do
	cd ${homeDir}
        echo Working on "$chr" intronic windows...
        # Make region.fa file for this chromosome, also convert from 0- to 1-base
        # Positive strand gene regions
        cat ${refDir}/GRCh37_MANE_Transcripts_Splice_Acceptor_Intronic_Windows.bed | \
		awk -v chr=${chr} '$1 == chr && $6=="+" {print $1":"$2+1"-"$3}' > \
			${homeDir}/OutputRegions/${chr}_+_Intronic_regions.txt

        # Negative strand gene regions
        cat ${refDir}/GRCh37_MANE_Transcripts_Splice_Acceptor_Intronic_Windows.bed | \
		awk -v chr=${chr} '$1 == chr && $6=="-" {print $1":"$2+1"-"$3}' > \
			${homeDir}/OutputRegions/${chr}_-_Intronic_regions.txt

        # Positive strand gene regions first
        samtools faidx -n 100 \
		${refDir}/${chr}.fa \
		-r ${homeDir}/OutputRegions/${chr}_+_Intronic_regions.txt > \
			${fastaDir}/${chr}_+_Intronic_regions.fa

        # Now negative strand gene regions
        samtools faidx -n 100 \
		-i \
		${refDir}/${chr}.fa \
		-r ${homeDir}/OutputRegions/${chr}_-_Intronic_regions.txt > \
			${fastaDir}/${chr}_-_Intronic_regions.fa

        # On to the trimer count script
        cd PerlTrimerOut/IntronicOut

        perl ${homeDir}/sequence_trimer_search/trimer_count.pl \
 		${homeDir}/sequence_trimer_search/uniq_trimers_AG_gain_loss_neutral.txt \
                ${fastaDir}/${chr}_+_Intronic_regions.fa \
                ${chr}_+_Intronic_region

        perl ${homeDir}/sequence_trimer_search/trimer_count.pl \
                ${homeDir}/sequence_trimer_search/uniq_trimers_AG_gain_loss_neutral.txt \
                ${fastaDir}/${chr}_-_Intronic_regions.fa \
                ${chr}_-_Intronic_region

        cd $homeDir
        echo Working on "$chr" intronic windows...... Done
done
