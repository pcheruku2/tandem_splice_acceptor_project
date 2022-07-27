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

# Firt 10 files
FileList=`\ls *vcf.gz | head`

BedDir="/home/frostfg/NAGNAG/genome_wide_trimer_code_V2/PerlTrimerOut"
TrimerRef="/home/frostfg/NAGNAG/genome_wide_trimer_code_V2/sequence_trimer_search/uniq_trimers_AG_gain_loss_neutral.txt"
scriptDir="/home/frostfg/NAGNAG/genome_wide_trimer_code_V2/vcf_AG_variation/subset_trimer_vcf.pl"
InputVCFdir="/home/frostfg/NAGNAG/ClinVar"
OutDir="/home/frostfg/NAGNAG/IntersectTest"


n=0
nTotal=`echo $FileList | wc`
for file in $FileList;
do

n=$((n+1))
echo "Working on ${file}... (${n}/${nTotal})"

# if file ends in "gz", unzip 

if echo $file | grep -q "gz$";
then
	FileStatus="Compressed"
	echo "File is compressed, unzipping..."
	gzip -d ${InputVCFDir}/${file}
	file=`echo $file | sed -e "s/\.vcf\.gz//g"`
else
	file=`echo $file | sed -e "s/\.vcf//g"`
fi

echo "Running intersection on ${file}..."
./VCF-NAGNAG_BED_Window_Intersection.sh  \
	--BedDir=${BedDir} \
	--TrimerRef=${TrimerRef} \
	--ScriptDir=${ScriptDir} \
	--InputVCFdir=${InputVCFdir} \
	--InputVCFfile=${file}.vcf \
	--OutDir=${OutDir} \
	--OutName=${file}

if echo $FileStatus | grep -q "^Compressed$";
        then
        
        echo "Rezipping file..."
        gzip  ${InputVCFDir}/${file}.vcf

fi

done



