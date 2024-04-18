#!/bin/sh
homeDir="/path/to/output/"
bamDir="/path/to/input/RNAseq/BAMs"
cd ${bamDir}

module load samtools

cd ${bamDir}
allBams=`\ls *bam`
cd ${homeDir}
mkdir AG-Gain_BAM_Files


n=0
nFinal=`echo ${allBams} | wc -w`
for bamFile in ${allBams}
do
	n=$((n+1))
	# Testing: bamFile="UDN515393_FB.hg19.bam"
	# Get the nagnag Variant information for all intronic gain variants this case
	awk -v bam=${bamFile} 'NR > 1 && $23 == bam && $18 == "gain" && $19 < 0 {print $0}' \
		All_NAGNAG_Vars_With_NewExonStart.tsv > \
			caseBamNagNag.tmp
	individualID=`cut -f10 caseBamNagNag.tmp | sort | uniq`
	echo "Checking for evidence of usage of intronic AG-Gain variants in ${bamFile} AKA ${individualID} (${n}/${nFinal})"

	# Check for alignments between the reference exon start and the predicted new exon start in the BAM file
	# This is being done for each intronic AG-gain variant individually
	# Make a destination directory for bam outputs of a given variant
	mkdir -p "AG-Gain_BAM_Files/${individualID}_${bamFile}_AG-gain-supporting_reads"
	caseBamOut="AG-Gain_BAM_Files/${individualID}_${bamFile}_AG-gain-supporting_reads"

	# This while loop takes the variant info and calculates the positions of the new exon and reference exon starts.
	# The loop then takes those two genomic positions and uses samtools view to pull reads that map within that span with MQ > 30
	while read line
	do
		echo ${line} | sed "s/ /\t/g" > caseBamNagNag_OneVariant.tmp
                variantName=`echo ${line} | sed "s/ /\t/g" | awk '{print $1"_"$2"_"$4"_"$5}'`

		caseBamFile=`cut -f23 caseBamNagNag_OneVariant.tmp`
		variantChr=`cut -f1 caseBamNagNag_OneVariant.tmp | sed "s/chr//g"`
		variantPos=`cut -f2 caseBamNagNag_OneVariant.tmp`
		variantStrand=`cut -f14 caseBamNagNag_OneVariant.tmp`
		variantNewExonBoundaryPosition=`cut -f22 caseBamNagNag_OneVariant.tmp`
		variantBoundaryPosition=`cut -f19 caseBamNagNag_OneVariant.tmp`
		
		# If the strand is + or - changes whether you need to add or subtract exon boundary positions (e.g. -13)
		if [ ${variantStrand} == "+" ]
		then
	        	variantOldExonStart=`expr ${variantPos} - ${variantBoundaryPosition}`
	        	variantNewExonStart=`expr ${variantOldExonStart} + ${variantNewExonBoundaryPosition}`
		else
                	variantOldExonStart=`expr ${variantPos} + ${variantBoundaryPosition}`
                	variantNewExonStart=`expr ${variantOldExonStart} - ${variantNewExonBoundaryPosition}`
		fi

		# Now Actually look for reads that map between the old, new exon starts
		# Filter these results map quality > 30 (0.1%), if needed go back to 20 (1%)
		if [ ${variantStrand} == "+" ]
		then
			samtools view \
				--min-MQ 30 \
				--bam \
				${bamDir}/${caseBamFile} \
				"${variantChr}:${variantNewExonStart}-${variantOldExonStart}" > \
					${caseBamOut}/${variantName}_${variantChr}_${variantStrand}_${variantNewExonStart}-${variantOldExonStart}_AG-Gain_MQ_30_Usage_Reads.bam
		else
			samtools view \
				--min-MQ 30 \
                                --bam \
				${bamDir}/${caseBamFile} \
				"${variantChr}:${variantNewExonStart}-${variantOldExonStart}" > \
                                        ${caseBamOut}/${variantName}_${variantChr}_${variantStrand}_${variantNewExonStart}-${variantOldExonStart}_AG-Gain_MQ_30_Usage_Reads.bam

		fi
	done < caseBamNagNag.tmp
	# Cleanup
	rm caseBamNagNag.tmp
	rm caseBamNagNag_OneVariant.tmp
done

