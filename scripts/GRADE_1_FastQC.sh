
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

step="FASTQC PRE TRIM"

echo "This script will run STEP 1 of the GRADE pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Pre-trimming quality control using FastQC"

#  1 FASTQC PRE TRIM <---
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
# 12 RSEM QUANTIFICATION
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
python=python/2.7.10
python3=python/3.6.1
fastqc=fastqc/0.11.8

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

if [ "$#" -gt 0 ]; then MainDirectory=$1; fi
if [ "$#" -eq 0 ]; then MainDirectory=`pwd`; fi

echo ""; echo "WARNING: The main directory for this run was set to ${MainDirectory}."; echo ""

############################################
##        CREATE OUTPUT DIRECTORIES       ##
############################################

# Go inside Main Directory (containing ONLY original fastq files)
cd $MainDirectory

# Create FastQC-PreTrim directory:
mkdir $MainDirectory/FastQCPreTrim


# 1 FASTQC PRE TRIM
############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Enter directory:
cd $MainDirectory

# Memory required (in GB):
memreq=20

# CPUs required (integer):
cpureq=12

############################################
##      CREATE PBS FILES FOR FASTQC       ##
############################################

# List all fastq files and create correspondent PBS files:
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do cp $MainDirectory/Header.pbs preTrimFastQC_${line}.pbs; done

# Populate PBS files (name of job):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "#PBS -N ${line}_1FastQC" >> preTrimFastQC_${line}.pbs; done

# Populate PBS files (number of threads):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "#PBS -r n" >> preTrimFastQC_${line}.pbs; done

# Populate PBS files (time and resources for running the job):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "#PBS -l mem=${memreq}GB,walltime=08:00:00,ncpus=${cpureq}" >> preTrimFastQC_${line}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "#PBS -m ae" >> preTrimFastQC_${line}.pbs; done

# Populate PBS files (e-mail address for correspondence):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "#PBS -M $email" >> preTrimFastQC_${line}.pbs; done

# Populate PBS files (empty line):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "" >> preTrimFastQC_${line}.pbs; done

# Populate PBS files (load module(s)):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "module load $fastqc" >> preTrimFastQC_${line}.pbs; done

# Populate PBS files (empty line):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "" >> preTrimFastQC_${line}.pbs; done

############################################
##         COMMAND LINE FOR FASTQC        ##
############################################
# Populate PBS files (actual command line):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "fastqc -t $cpureq --outdir ${MainDirectory}/FastQCPreTrim $MainDirectory/${line}_R1.fq.gz $MainDirectory/${line}_R2.fq.gz" >> preTrimFastQC_${line}.pbs ; done

############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> FastQC Pre Trimming" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n1 preTrimFastQC_*.pbs >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "-- Job Submission IDs --" >> $MainDirectory/Pipeline_Log

############################################
##            SUBMIT FASTQC JOBS          ##
############################################
# Submit jobs and populate Log file (submission ID numbers):
ls preTrimFastQC_*.pbs | while read line; do echo "$line" >> $MainDirectory/Pipeline_Log ; qsub "$line"; done >> $MainDirectory/Pipeline_Log

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user"


############################################
##                   END                  ##
############################################

