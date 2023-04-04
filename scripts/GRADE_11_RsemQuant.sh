
#!/bin/sh

############################################
############################################
#                                          #
#         **GRADE in PBS Cluster**         #
#        General RNAseq Analysis for       #
#         Differential Expression          #
#                                          #
############################################
############################################
# #        Version 4.0 - APR 2023        # #
############################################


############################################
##             PIPELINE STEPS             ##
############################################

step="RSEM QUANTIFICATION"

echo "This script will run STEP 11 of the GRADE pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Quantify counts with Rsem"

#  1 FASTQC PRE TRIM
#  2 TRIMMOMATIC
#  3 FASTQC POST TRIM
#  4 KALLISTO INDEXING (if needed)
#  5 KALLISTO QUANTIFICATION
#  6 Create KALLISTO count tables
#  7 STAR INDEXING (if needed)
#  8 STAR ALIGNMENT
#  9 SAMTOOLS
# 10 NOVOSORT
# 11 RSEM INDEXING (if needed)
# 12 RSEM QUANTIFICATION <---
# 13 Create STAR/RSEM count tables


############################################
##            EXAMPLES OF DATA            ##
############################################

# Data should be named as: ID1_ID2_R1.fq.gz (forward) and ID1_ID2_R2.fq.gz (reverse)
# For example:
# Cell1_Sample1_R1.fq.gz  Cell1_Sample2_R1.fq.gz  Cell2_Sample1_R1.fq.gz  Cell2_Sample2_R1.fq.gz  Cell3_Sample1_R1.fq.gz  Cell3_Sample2_R1.fq.gz
# Cell1_Sample1_R2.fq.gz  Cell1_Sample2_R2.fq.gz  Cell2_Sample1_R2.fq.gz  Cell2_Sample2_R2.fq.gz  Cell3_Sample1_R2.fq.gz  Cell3_Sample2_R2.fq.gz

############################################
##            DECLARE  MODULES            ##
############################################
rsem=RSEM/1.2.30

############################################
##          REFERENCE  DIRECTORY          ##
############################################
RefDirectory=/working/lab_julietF/mainaB/ReferenceGenomes


############################################
##              SET USER VARs             ##
############################################
# Edit the following line(s) to reflect your environment

# Email to which messages will be sent when each PBS job is initiated/concluded/terminated
email=Maina.Bitar@qimrberghofer.edu.au

# User name within your cluster environment
user=`whoami`

############################################
##         SET MAIN WORK DIRECTORY        ##
############################################
# Define directory:

if [ "$#" -gt 1 ]; then MainDirectory=$2; fi
if [ "$#" -eq 1 ]; then MainDirectory=`pwd`; fi

echo ""; echo "WARNING: The main directory for this run was set to ${MainDirectory}."; echo ""

############################################
##        CREATE OUTPUT DIRECTORIES       ##
############################################

# Go inside Main Directory (containing ONLY original fastq files)
cd $MainDirectory

# Create Star directory:
mkdir $MainDirectory/RsemQuant


# 11 RSEM QUANT
############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Enter directory:
cd $MainDirectory/SamtoolsNovosort

# Memory required (in GB):
memreq=8

# CPUs required (integer):
cpureq=12


############################################
##       CREATE PBS FILES FOR RSEM        ##
############################################
# List all fastq files and create correspondent PBS files:
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do cp $MainDirectory/Header.pbs Rsem_${line}.pbs; done

# Populate PBS files (name of job):
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "#PBS -N ${line}_Rsem" >> Rsem_${line}.pbs; done

# Populate PBS files (number of threads):
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "#PBS -r n" >> Rsem_${line}.pbs; done

# Populate PBS files (time and resources for running the job):
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "#PBS -l mem=${memreq}GB,walltime=08:00:00,ncpus=${cpureq}" >> Rsem_${line}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "#PBS -m ae" >> Rsem_${line}.pbs; done

# Populate PBS files (e-mail address for correspondence):
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "#PBS -M $email" >> Rsem_${line}.pbs; done

# Populate PBS files (empty line):
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "" >> Rsem_${line}.pbs; done

# Populate PBS files (load module(s)):
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "module load $rsem" >> Rsem_${line}.pbs; done

# Populate PBS files (empty line):
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "" >> Rsem_${line}.pbs; done


############################################
##          COMMAND LINE FOR RSEM         ##
############################################
# ENSEMBL:
if [ "$1" == ensembl ] || [ "$1" == E ] ; then
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "rsem-calculate-expression --paired-end --bam --forward-prob 0 --no-bam-output -p $cpureq $MainDirectory/SamtoolsNovosort/${line}.novosort.bam $RefDirectory/RsemIndex/EnsemblTranscriptome $MainDirectory/RsemQuant/${line}" >> Rsem_${line}.pbs; done
fi

# GENCODE
if [ "$1" == gencode ] || [ "$1" == G ] ; then
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "rsem-calculate-expression --paired-end --bam --forward-prob 0 --no-bam-output -p $cpureq $MainDirectory/SamtoolsNovosort/${line}.novosort.bam $RefDirectory/RsemIndex/GencodeTranscriptome $MainDirectory/RsemQuant/${line}" >> Rsem_${line}.pbs; done
fi

############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "==> Rsem Quantification" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- Command Lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n1 Rsem_*.pbs >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "-- Job Submission IDs --" >> $MainDirectory/Pipeline_Log

############################################
##             SUBMIT Rsem JOBS           ##
############################################
# Submit jobs and populate Log file (submission ID numbers):
ls Trimmed_*.novosort.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do qsub Rsem_${line}.pbs >> $MainDirectory/Pipeline_Log; done

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user"


############################################
##                   END                  ##
############################################


