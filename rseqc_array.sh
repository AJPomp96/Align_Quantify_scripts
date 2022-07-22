#!/bin/bash
#SBATCH --output=./slurm-out/rseqc-%A_%a.out
#SBATCH --job-name=RSEQC
#SBATCH --mem=32000
#SBATCH -c 2

#JOB LOG HEADER
perl -E 'say"="x80'; echo "JOB STARTED: `date`"; echo "NODE: `hostname`"; echo "SCRIPT ${0}:"; echo "JOB ID: ${SLURM_JOB_ID}"; cat $0; perl -E 'say"="x80'

#SOFTWARE REQUIREMENTS
cd ${ALIGNDIR}

#VARIABLES
PREFIXES=($(ls -1 ${FQDIR} | grep $GREP | sed -r $SED | uniq))

INPUT=${PREFIXES[$SLURM_ARRAY_TASK_ID]}_sorted_rf_alignment.bam

#COMMAND(s) TO RUN
echo "PROCESSING FILE: ${INPUT}"

if ls -1 $FQDIR | grep -q "R[2]"
then
  echo "fastq is paired end"

else
  echo "fastq is single end"
  read_distribution.py\
  -i $INPUT\
  -r ${BEDPATH}\
  > ${RSEQCDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}_read_dists.txt

  read_GC.py\
  -i $INPUT\
  -o ${RSEQCDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}

  bam_stat.py\
  -i $INPUT\
  > ${RSEQCDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}_bam_stats.txt

  infer_experiment.py\
  -i $INPUT\
  -r $BEDPATH\
  > ${RSEQCDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}_inf_exp.txt

fi

#JOB LOG FOOTER
perl -E 'say"="x80'; echo "JOB COMPLETED: `date`"; perl -E 'say"="x80'
