#!/bin/bash
#SBATCH --output=./slurm-out/slurm-%A_%a.out
#SBATCH --job-name=prefastqc
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

#get quality report of fastq files pretrim
if ls -1 $FQDIR | grep -q "R[12]"
then
  echo "fastq is paired end"
  FILE1=${FILES[0]}
  FILE2=${FILES[0]}
  echo $FILE1,$FILE2
  fastqc ${FILE1} -o ${PRETRIM_QC}
  fastqc ${FILE2} -o ${PRETRIM_QC}
else
  echo "fastq is single end"
  FILE1=${FILES[0]}
  echo $FILE1
  fastqc ${FILE1} -o ${PRETRIM_QC}
fi

#JOB LOG FOOTER
perl -E 'say"="x80'; echo "JOB COMPLETED: `date`"; perl -E 'say"="x80'
