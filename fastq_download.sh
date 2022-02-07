#!/bin/bash
#
#File:
#Author:
#Purpose:
#
#
#Created:

#SBATCH --output=./slurm-out/slurm-%j.out
#SBATCH --job-name=fastqdownload
#SBATCH --mem=10000
#SBATCH -c 8

#JOB LOG HEADER
perl -E 'say"="x80'; echo "JOB STARTED: `date`"; echo "NODE: `hostname`"; echo "SCRIPT ${0}:"; echo "JOB ID: ${SLURM_JOB_ID}"; cat $0; perl -E 'say"="x80'

#SOFTWARE REQUIREMENTS

#VARIABLES

#COMMAND(s) TO RUN
python fastq_download.py
#JOB LOG FOOTER
perl -E 'say"="x80'; echo "JOB COMPLETED: `date`"; perl -E 'say"="x80'
