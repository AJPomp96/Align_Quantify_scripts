#!/bin/bash
#SBATCH --output=./slurm-out/map-%A_%a.out
#SBATCH --job-name=map
#SBATCH --mem=32000
#SBATCH -c 16

#JOB LOG HEADER
perl -E 'say"="x80'; echo "JOB STARTED: `date`"; echo "NODE: `hostname`"; echo "SCRIPT ${0}:"; echo "JOB ID: ${SLURM_JOB_ID}"; cat $0; perl -E 'say"="x80'

#SOFTWARE REQUIREMENTS

#VARIABLES
PREFIXES=($(ls -1 ${FQDIR} | grep $GREP | sed -r $SED | uniq))

#COMMAND(s) TO RUN
cd ${TRIMDIR}

FILES=($(ls -1 | grep "${PREFIXES[$SLURM_ARRAY_TASK_ID]}.*ribotrim"))

if ls -1 $FQDIR | grep -q "R[12]"
then
  echo "fastq is paired end"
  FILE1=${FILES[0]}
  FILE2=${FILES[1]}
  echo "PROCESSING ${FILE1} and ${FILE2}"
  hisat2 -p 8 \
        --verbose \
        --phred33 \
        --dta \
        --fr \
        --summary-file ${ALIGNDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}_AlignStat.txt \
        -x ${HISAT2_INDEXES}/${HISAT2_PREFIX} \
        -1 ${FILE1} \
        -2 ${FILE2} \
	--rna-strandness RF \
        -S ${ALIGNDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}_aligned_reads.sam
else
  echo "fastq is single end"
  FILE1=${FILES[0]}
  echo "PROCESSING ${FILE1}"
  hisat2 -p 8 \
        --verbose \
        --phred33 \
        --dta \
        --fr \
        --summary-file ${ALIGNDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}_AlignStat.txt \
        -x ${HISAT2_INDEXES}/${HISAT2_PREFIX} \
        -U ${FILE1} \
	--rna-strandness RF \
        -S ${ALIGNDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}_aligned_reads.sam
fi

#JOB LOG FOOTER
perl -E 'say"="x80'; echo "JOB COMPLETED: `date`"; perl -E 'say"="x80'
