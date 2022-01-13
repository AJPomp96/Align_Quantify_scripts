#!/bin/bash
#SBATCH --output=./slurm-out/trim-%A_%a.out
#SBATCH --job-name=trim
#SBATCH --mem=8000
#SBATCH -c 8

#JOB LOG HEADER
perl -E 'say"="x80'; echo "JOB STARTED: `date`"; echo "NODE: `hostname`"; echo "SCRIPT ${0}:"; echo "JOB ID: ${SLURM_JOB_ID}"; cat $0; perl -E 'say"="x80'

#SOFTWARE REQUIREMENTS

#VARIABLES
PREFIXES=($(ls -1 ${FQDIR} | grep $GREP | sed -r $SED | uniq))

#COMMAND(s) TO RUN
cd ${FQDIR}
FILES=($(ls -1 | grep "${PREFIXES[$SLURM_ARRAY_TASK_ID]}.*\.fastq.gz"))

if ls -1 $FQDIR | grep -q "R[12]"
then
  echo "fastq is paired end"
  FILE1=${FILES[0]}
  FILE2=${FILES[1]}
  echo "PROCESSING ${FILE1} and ${FILE2}"
  trim_galore -o ${TRIMDIR} \
        --length 35 \
        --quality 28 \
        --paired \
        --phred33 \
        --clip_R2 3 \
        --cores 8 \
        ${FILE1}\
        ${FILE2}
else
  echo "fastq is single end"
  FILE1=${FILES[0]}
  echo "PROCESSING ${FILE1}"
  trim_galore -o ${TRIMDIR} \
        --length 35 \
        --quality 28 \
        --phred33 \
        --cores 8 \
        ${FILE1}
fi

#JOB LOG FOOTER
perl -E 'say"="x80'; echo "JOB COMPLETED: `date`"; perl -E 'say"="x80'
