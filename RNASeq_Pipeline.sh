#!/bin/bash
#
#File:Paired_End_Pipeline.sh
#Author: Anthony Pompetti
#Purpose:
#
#
#Created:

#SBATCH --output=./slurm-out/main-%j.out
#SBATCH --job-name=RNA_Seq_Pipeline
#SBATCH --mem=8000
#SBATCH -c 1

#JOB LOG HEADER
perl -E 'say"="x80'; echo "JOB STARTED: `date`"; echo "NODE: `hostname`"; echo "SCRIPT ${0}:"; echo "JOB ID: ${SLURM_JOB_ID}"; cat $0; perl -E 'say"="x80'

#SOFTWARE REQUIREMENTS

#VARIABLES
#Define Main Directories
export PATH=${PATH}:$(pwd)/scripts/Align_Quantify_scripts  # Add this pipeline to the executeable path
export FQDIR=$(pwd)/fastq
export TRIMDIR=$(pwd)/Trimmed
export ALIGNDIR=$(pwd)/Alignments
export COUNTDIR=$(pwd)/Counts
export PRETRIM_QC=$(pwd)/PreTrim_FastQC
export POSTTRIM_QC=$(pwd)/PostTrim_FastQC
export KLSTODIR=$(pwd)/Kallisto
export RSEQCDIR=$(pwd)/RSeQC_Results
#export KLSTOIDX=/work/ajpompet/EnsMm_grc39_104/EnsMm_104_Kallisto_total
export GTFPATH=/work/ajpompet/EnsMm_grc39_104/Mus_musculus.GRCm39.104.gtf
export DEXPATH=/work/ajpompet/EnsMm_grc39_104/EnsMm_104_DEXSeq_Annot.gff
export BEDPATH=/work/ajpompet/EnsMm_grc39_104/rseqc_gene_models.bed
export HISAT2_INDEXES=/work/ajpompet/EnsMm_grc39_104
export HISAT2_PREFIX=EnsMm_grc39_104
#Detect if fastq files are paired end or single end reads
if ls -1 $FQDIR | grep -q "R[12]"
then
  echo "fastq is paired end"
  export GREP="L[0-9]\{3\}_R1_[0-9]\{3\}\.fastq\.gz"
  export SED='s/_R[12]_001[.]fastq.gz//'
  export PREFIXES=($(ls -1 ${FQDIR} | grep $GREP | sed -r $SED | uniq))
else
  echo "fastq is single end"
  export GREP=".fastq.gz"
  export SED='s/[.]fastq.gz//'
  export PREFIXES=($(ls -1 $FQDIR | grep $GREP | sed -r $SED | uniq))
fi
export ENDINDX=$((${#PREFIXES[@]} - 1))

#Mark DELOLD as 1 if you wish to delete old data in folders
export DELOLD=0

#COMMAND(s) TO RUN
#Create Directories
DIRARR=($TRIMDIR\
	 $ALIGNDIR \
	 $COUNTDIR \
	 $PRETRIM_QC \
	 $POSTTRIM_QC \
	 $KLSTODIR \
	 $RSEQCDIR)

for dir in ${DIRARR[@]}
do
  if [ ! -d $dir ]
  then
    mkdir $dir
  
  elif [ $DELOLD == 1 ];
  then
    rm -rf $dir
    mkdir $dir
  fi
done

echo ${PREFIXES[@]}
echo $GREP
echo $SED
echo $ENDINDX

#Trimgalore step
trimJB=$(sbatch --array [0-$ENDINDX] trim_array.sh | gawk '{print $4}')

#Pretrim FastQC step
prefastqcJB=$(sbatch --array [0-$ENDINDX]%3 pretrimfastqc_array.sh | gawk '{print $4}')

#Align to rRNA step
alignRiboJB=$(sbatch --array [0-$ENDINDX] --dependency=afterok:$trimJB ribomap_array.sh | gawk '{print $4}')

#rRNA sort step
sortRiboJB=$(sbatch --array [0-$ENDINDX] --dependency=afterok:$alignRiboJB ribosort_array.sh | gawk '{print $4}')

#Align rRNA-Trimmed reads step
alignJB=$(sbatch --array [0-$ENDINDX] --dependency=afterok:$sortRiboJB map_array.sh | gawk '{print $4}')

#Posttrim FastQC step
postfastqcJB=$(sbatch --array [0-$ENDINDX]%3 --dependency=afterok:$sortRiboJB posttrimfastqc_array.sh | gawk '{print $4}')

#Sort rRNA-Trimmed reads step
sortJB=$(sbatch --array [0-$ENDINDX] --dependency=afterok:$alignJB sort_array.sh | gawk '{print $4}')

#Perform RSeQC modules
rseqcJB=$(sbatch --array [0-$ENDINDX]%3 --dependency=afterok:$sortJB rseqc_array.sh | gawk '{print $4}')

#Gene-level Counts Htseq-count step
geneCountJB=$(sbatch --array [0-$ENDINDX]%3 --dependency=afterok:$sortJB count_array.sh)

#JOB LOG FOOTER
perl -E 'say"="x80'; echo "JOB COMPLETED: `date`"; perl -E 'say"="x80'
