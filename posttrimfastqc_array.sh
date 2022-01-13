#!/bin/bash
#SBATCH --output=./slurm-out/pofqc-%A_%a.out
#SBATCH --job-name=postfastqc
#SBATCH --mem=8000
#SBATCH -c 8

#JOB LOG HEADER
perl -E 'say"="x80'; echo "JOB STARTED: `date`"; echo "NODE: `hostname`"; echo "SCRIPT ${0}:"; echo "JOB ID: ${SLURM_JOB_ID}"; cat $0; perl -E 'say"="x80'

#SOFTWARE REQUIREMENTS

#VARIABLES
PREFIXES=($(ls -1 ${FQDIR} | grep $GREP | sed -r $SED | uniq))

#COMMAND(s) TO RUN
cd ${TRIMDIR}

#get quality report of fastq files posttrim
if ls -1 $FQDIR | grep -q "R[12]"
then
  echo "fastq is paired end"
  fastqc ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_ribotrim_R1.fq -o ${POSTTRIM_QC}
  fastqc ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_ribotrim_R2.fq -o ${POSTTRIM_QC}
else
  echo "fastq is single end"
  fastqc ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_ribotrim.fq -o ${POSTTRIM_QC}
fi

#JOB LOG FOOTER
perl -E 'say"="x80'; echo "JOB COMPLETED: `date`"; perl -E 'say"="x80'
