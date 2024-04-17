#!/bin/sh

##############################################
# Option values for testing purposed:
# BedDir='/home/frostfg/NAGNAG/genome_wide_trimer_code_V2/PerlTrimerOut'
# TrimerRef='/home/frostfg/NAGNAG/genome_wide_trimer_code_V2/sequence_trimer_search/uniq_trimers_AG_gain_loss_neutral.txt'
# ScriptDir='/home/frostfg/NAGNAG/genome_wide_trimer_code_V2/vcf_AG_variation/subset_trimer_vcf.pl'
# InputVCFdir='/home/frostfg/NAGNAG/ClinVar'
# InputVCFfile='clinvar.vcf'
# OutDir='/home/frostfg/NAGNAG/IntersectTest'
# OutName='ClinVar_Example'
##############################################


# Initialize variables that are getting passed in from options
BedDir='' # This is where the Trimer Site BED files are (The output from running trimer_count.pl)
TrimerRef='' # The location of uniq_trimers_AG_gain_loss_neutral.txt
ScriptDir='' # The location of subset_trimer_vcf.pl
InputVCFdir='' # The directory containing your input VCF. NOTE: do not put a slash after the directory name
InputVCFfile='' # The name of the inout VCF file
OutDir='' # The output directory where the intersections will go
OutName='' # The sample name you want included in the  output file names


# Here's how those options get passed into the bash script
# Shoutout to https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
# This code is taken from the above link, and isn't actually mine. Someone who actually knows what they're doing wrote it, so hopefully it works.
# Author's Note: This code does actually work quite well. Thanks internet!

for i in "$@"
do
case $i in
    -b=*|--BedDir=*)
    BedDir="${i#*=}"
    shift # past argument=value
    ;;
    -T=*|--TrimerRef=*)
    TrimerRef="${i#*=}"
    shift # past argument=value
    ;;
    -s=*|--ScriptDir=*)
    ScriptDir="${i#*=}"
    shift # past argument=value
    ;;
    -i=*|--InputVCFdir=*)
    InputVCFdir="${i#*=}"
    shift # past argument=value
    ;;
    -I=*|--InputVCFfile=*)
    InputVCFfile="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--OutDir=*)
    OutDir="${i#*=}"
    shift # past argument=value
    ;;
    -O=*|--OutName=*)
    OutName="${i#*=}"
    shift # past argument=value
    ;;
esac
done

# Need to make a variable with all the chromosomes we need to loop through
AllChroms=`echo "chr"{1..22}`
AllChroms=`echo $AllChroms "chrX" "chrY"`

# I need to make two temporary directories for the chromosome-by-chromsome results
cd $OutDir
mkdir $OutDir"/tmp_VCFexonic"
mkdir $OutDir"/tmp_VCFintronic"

for chr in $AllChroms;
do

	echo working on exonic windows for "$chr"...
	
	perl ${ScriptDir} \
                ${TrimerRef} \
                ${BedDir}/ExonicOut/trimer_sites_${chr}_+_Exonic_region.bed \
                ${InputVCFdir}/${InputVCFfile} >  ${OutDir}/tmp_VCFexonic/${OutName}_${chr}_+_exonic_trimer_vcf_AG_variation.tsv
	
	perl ${ScriptDir}\
                ${TrimerRef}\
                ${BedDir}/ExonicOut/trimer_sites_${chr}_-_Exonic_region.bed\
                $InputVCFdir/${InputVCFfile} >  ${OutDir}/tmp_VCFexonic/${OutName}_${chr}_-_exonic_trimer_vcf_AG_variation.tsv

	echo working on intronic windows for "$chr"...
        perl ${ScriptDir}\
                ${TrimerRef}\
                ${BedDir}/IntronicOut/trimer_sites_${chr}_+_Intronic_region.bed\
                ${InputVCFdir}/${InputVCFfile} >  ${OutDir}/tmp_VCFintronic/${OutName}_${chr}_+_intronic_trimer_vcf_AG_variation.tsv

        perl ${ScriptDir}\
                ${TrimerRef}\
                ${BedDir}/IntronicOut/trimer_sites_${chr}_-_Intronic_region.bed\
                ${InputVCFdir}/${InputVCFfile} >  ${OutDir}/tmp_VCFintronic/${OutName}_${chr}_-_intronic_trimer_vcf_AG_variation.tsv

done

# Paste all the chromosome-wise results into one big .tsv file
cat ${OutDir}/tmp_VCFexonic/* > ${OutDir}/${OutName}_Exonic_Window_AG_Variants.tsv
cat ${OutDir}/tmp_VCFintronic/* > ${OutDir}/${OutName}_Intronic_Window_AG_Variants.tsv

rm -r "$OutDir"/tmp_VCFexonic
rm -r "$OutDir"/tmp_VCFintronic
