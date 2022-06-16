#!/bin/sh

homeDir='/home/frostfg/NAGNAG/'

cd $homeDir
module load bcftools
# Testing Chr13
# awk '$1 == 13 {print $1 "\t" $2 "\t" $3 "\t" $4}' genome_wide_trimer_code/AllChr_NAGNAG_Exonic_Windows.bed > chr13exonicTest.bed
# bcftools view --regions-file chr13exonicTest.bed ClinVar/clinvar.vcf.gz > chr13Exonictest.vcf
# perl subset_trimer_vcf.pl uniq_trimers_AG_gain_loss_neutral.txt chr13exonicTest.bed chr13Exonictest.vcf > chr13_Exonic_Test_Clinvar_Intersection.vcf



# Extract all variation that intersects the second-base (middle base) in the trimer
echo Working on intersection of exonic window
bcftools view --regions-file genome_wide_trimer_code_V2/AllChr_NAGNAG_Exonic_Windows.bed  ClinVar/clinvar.vcf.gz > All_trimer_sites_Exonic_clinvar.vcf
echo Working on intersection of intronic window
bcftools view --regions-file genome_wide_trimer_code_V2/AllChr_NAGNAG_Intronic_Windows.bed ClinVar/clinvar.vcf.gz > All_trimer_sites_Intronic_clinvar.vcf

# Ensure the variant allele results either in gain, loss, or neutral change (since not all variation at the second base in the trimer lead to AG site gain, loss, or neutral change.
echo Determining AG gains and losses for exonic window 
perl genome_wide_trimer_code_V2/vcf_AG_variation/subset_trimer_vcf.pl uniq_trimers_AG_gain_loss_neutral.txt genome_wide_trimer_code_V2/AllChr_NAGNAG_Exonic_Windows.bed All_trimer_sites_Exonic_clinvar.vcf > AllChr_Exonic_Clinvar_Intersection.vcf
echo Determining AG gains and losses for intronic window 
perl genome_wide_trimer_code_V2/vcf_AG_variation/subset_trimer_vcf.pl uniq_trimers_AG_gain_loss_neutral.txt genome_wide_trimer_code_V2/AllChr_NAGNAG_Intronic_Windows.bed All_trimer_sites_Intronic_clinvar.vcf > AllChr_Intronic_Clinvar_Intersection.vcf


