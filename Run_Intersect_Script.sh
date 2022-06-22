#!/bin/sh

# Useful sidenote.
# Any command with multiple options/flags can be run in 2 ways in bash:
# Method 1:
# command --option=variable --output=/some/output/directory --method=SomethingOvercomplicated input.txt
# However, method one gets hard to read (and therefore hard to spot typos) if you have long file paths as inputs, or many flags.
# Method 2 offers an alternative:
# command \
# --option=variable \
# --output=/some/output/directory \
# --method=SomethingOvercomplicated \
# input.txt
# The backslashes offer a lot better organization and legibility in my opinion, which is why I almost always use them for functions with several inputs.
# Another note here: the tab indent is not necessary, I just think it looks better/is easier to read



./VCF-NAGNAG_BED_Window_Intersection.sh  --BedDir='/home/frostfg/NAGNAG/genome_wide_trimer_code_V2/PerlTrimerOut' \
	--TrimerRef='/home/frostfg/NAGNAG/genome_wide_trimer_code_V2/sequence_trimer_search/uniq_trimers_AG_gain_loss_neutral.txt' \
	--ScriptDir='/home/frostfg/NAGNAG/genome_wide_trimer_code_V2/vcf_AG_variation/subset_trimer_vcf.pl' \
	--InputVCFdir='/home/frostfg/NAGNAG/ClinVar' \
	--InputVCFfile='clinvar.vcf' \
	--OutDir='/home/frostfg/NAGNAG/IntersectTest' \
	--OutName='ClinVar_Example'

