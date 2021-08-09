#!/bin/bash
#SBATCH --output=./slurm-out/slurm-%A_%a.out
#SBATCH --job-name=count
#SBATCH --mem=32000
#SBATCH -c 1

#JOB LOG HEADER
perl -E 'say"="x80'; echo "JOB STARTED: `date`"; echo "NODE: `hostname`"; echo "SCRIPT ${0}:"; echo "JOB ID: ${SLURM_JOB_ID}"; cat $0; perl -E 'say"="x80'

#SOFTWARE REQUIREMENTS

#VARIABLES
PREFIXES=($(ls -1 ${FQDIR} | grep $GREP | sed -r $SED | uniq))
INPUT=${PREFIXES[$SLURM_ARRAY_TASK_ID]}_sorted_rf_alignment.bam

#COMMAND(s) TO RUN
cd ${ALIGNDIR}

echo "PROCESSING FILE: ${INPUT}"
htseq-count \
	-i gene_id \
	-r pos \
	-f bam \
	-s reverse \
	-m union \
	--type exon \
	${INPUT} \
	${GTFPATH} > ${COUNTDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}_rf_GeneCount.txt

#JOB LOG FOOTER
perl -E 'say"="x80'; echo "JOB COMPLETED: `date`"; perl -E 'say"="x80'
