#!/bin/sh
#SBATCH -J NAGNAG_CIGAR
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=graeme.frost@nih.gov
#SBATCH --time=01:00:00
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=144
#SBATCH --mem=750g
#SBATCH --gres=lscratch:50

homeDir="/path/to/output"
bamOutDir="/path/to/AG-Gain/BAM/fles"
cigarScriptDir="/path/to/cigarRscript"

module load samtools
module load python
module load R
module load parallel

# Make an output directory for the CIGAR decoding
cd ${homeDir}
mkdir AG-Gain_CIGAR_Decoding
cigarOutDir="${homeDir}/AG-Gain_CIGAR_Decoding/AllSamples_Intronic_AG_Gain_Variant_Usage"

# This is where the for loop to go through all cases would start
nFinal=`\ls ${bamOutDir} | wc -l`
n=0
loopStartTime=`date +%s`
for caseName in `\ls -1d ${bamOutDir}/U* | xargs -n1 basename`
do
	cd ${bamOutDir}
	n=$( (n + 1) )
	echo "Evaluating intronic AG-gain variants in ${caseName} (${n}/${nFinal})"
	startTime=`date +%s`
	cd ${caseName}

	# Index all bam files with GNU parallel
	# NOTE: empty bam files are 3.1k, so only indexing and searching files bigger than that will save some time (hopefully)
	echo "Indexing variant-specific BAM files"
	find . -maxdepth 1 -type f -size +3k | \
		egrep "bam$" | \
		parallel 'samtools index {}'

	# Initiate header for file R writes to 
	caseID=`echo ${caseName} | grep -oe "UDP[0-9]*"`

	# Search for exon boundaries and check them against the variant
	export TMPDIR=/lscratch/${SLURM_JOB_ID}
	find . -maxdepth 1 -type f -size +3k | \
                egrep "bam$" | \
		sed "s/\.\///g" | \
		parallel  Rscript ${cigarScriptDir}/all_AG-Gains_bam_CIGAR_extract_and_intersect.R `pwd`/{} ${caseID}

        endTime=`date +%s`
	runTime=$((endTime - startTime))
	nLeft=$((nFinal - n))
	etaSec=$(($(($((endTime - loopStartTime))/n))*nLeft))
	etaMin=$((etaSec / 60))
	etaRemainder=$((etaSec %  60))
	echo "completed in ${runTime} seconds, ETA: ${etaMin}min, ${etaRemainder}sec"
	echo
done

