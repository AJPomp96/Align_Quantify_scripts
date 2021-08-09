#!/bin/bash
#SBATCH --output=./slurm-out/slurm-%A_%a.out
#SBATCH --job-name=sort
#SBATCH --mem=64000
#SBATCH -c 16

#JOB LOG HEADER
perl -E 'say"="x80'; echo "JOB STARTED: `date`"; echo "NODE: `hostname`"; echo "SCRIPT ${0}:"; echo "JOB ID: ${SLURM_JOB_ID}"; cat $0; perl -E 'say"="x80'

#SOFTWARE REQUIREMENTS

#VARIABLES
PREFIXES=($(ls -1 ${FQDIR} | grep $GREP | sed -r $SED | uniq))
FILES=($(ls -1 | grep "${PREFIXES[$SLURM_ARRAY_TASK_ID]}.*\.sam"))

#COMMAND(s) TO RUN
cd ${ALIGNDIR}
echo $FILES

#Generate BAM file
FILE1=${PREFIXES[$SLURM_ARRAY_TASK_ID]}_aligned_reads.sam
echo "PROCESSING FILE: ${FILE1}"
samtools view \
	-@ 8 \
	-bS ${FILE1} > ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_aligned_reads.bam

#Sort generated BAM file
samtools sort \
	-@ 8 \
	-m 12G \
	-o ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_sorted_rf_alignment.bam ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_aligned_reads.bam

rm ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_aligned_reads.bam ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_aligned_reads.sam

#Index sorted BAM file
samtools index ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_sorted_rf_alignment.bam

#JOB LOG FOOTER
perl -E 'say"="x80'; echo "JOB COMPLETED: `date`"; perl -E 'say"="x80'
