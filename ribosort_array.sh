#!/bin/bash
#SBATCH --output=./slurm-out/ribos-%A_%a.out
#SBATCH --job-name=ribosort
#SBATCH --mem=64000
#SBATCH -c 8

#JOB LOG HEADER
perl -E 'say"="x80'; echo "JOB STARTED: `date`"; echo "NODE: `hostname`"; echo "SCRIPT ${0}:"; echo "JOB ID: ${SLURM_JOB_ID}"; cat $0; perl -E 'say"="x80'

#SOFTWARE REQUIREMENTS

#VARIABLES
PREFIXES=($(ls -1 ${FQDIR} | grep $GREP | sed -r $SED | uniq))

#COMMAND(s) TO RUN
cd ${ALIGNDIR}

FILES=($(ls -1 | grep "${PREFIXES[$SLURM_ARRAY_TASK_ID]}.*\.sam"))

FILE1=${FILES[0]}

echo "PROCESSING FILE: ${FILE1}"

#Generate BAM file from SAM file
samtools view \
	-@ 8 \
	-bS ${FILE1} > ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_ribo_reads.bam

#Sort the generated BAM file
samtools sort \
	-n \
	-@ 8 \
	-m 12G \
	-o ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_byname_ribo.bam ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_ribo_reads.bam

#Remove the ribo_reads bam/sam file
rm ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_ribo_reads.bam ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_ribo_reads.sam

#Generate Fastq files that results in Fastq files without rRNA reads
if ls -1 $FQDIR | grep -q "R[12]"
then
  echo "fastq is paired end"
  bedtools bamtofastq -i <(samtools view -h -f 13 ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_byname_ribo.bam) \
  -fq ${TRIMDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}_ribotrim_R1.fq \
  -fq2 ${TRIMDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}_ribotrim_R2.fq
else
  echo "fastq is single end"
  bedtools bamtofastq -i <(samtools view -h -f 4 ${PREFIXES[$SLURM_ARRAY_TASK_ID]}_byname_ribo.bam) \
  -fq ${TRIMDIR}/${PREFIXES[$SLURM_ARRAY_TASK_ID]}_ribotrim.fq
fi

#JOB LOG FOOTER
perl -E 'say"="x80'; echo "JOB COMPLETED: `date`"; perl -E 'say"="x80'
