
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

step="NOVOSORT"

echo "This script will run STEP 10 of the GRADE pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Sort Star alignments with Novosort"

#  1 FASTQC PRE TRIM
#  2 TRIMMOMATIC
#  3 FASTQC POST TRIM
#  4 KALLISTO INDEXING (if needed)
#  5 KALLISTO QUANTIFICATION
#  6 Create KALLISTO count tables
#  7 STAR INDEXING (if needed)
#  8 STAR ALIGNMENT
#  9 SAMTOOLS
# 10 NOVOSORT <---
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
novosort=novoalign/3.05.01

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


# 10 NOVOSORT
############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Enter directory:
cd $MainDirectory/SamtoolsNovosort

# Memory required (in GB):
memreq=20

# CPUs required (integer):
cpureq=8

# Internal memory for sort buffer:
novomem=16


############################################
##     CREATE PBS FILES FOR NOVOSORT      ##
############################################
# List all fastq files and create correspondent PBS files:
ls Trimmed_*.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do cp $MainDirectory/Header.pbs Novosort_${line}.pbs; done

# Populate PBS files (name of job):
ls Trimmed_*.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "#PBS -N ${line}_Novosort" >> Novosort_${line}.pbs; done

# Populate PBS files (number of threads):
ls Trimmed_*.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "#PBS -r n" >> Novosort_${line}.pbs; done

# Populate PBS files (time and resources for running the job):
ls Trimmed_*.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "#PBS -l mem=${memreq}GB,walltime=08:00:00,ncpus=${cpureq}" >> Novosort_${line}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
ls Trimmed_*.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "#PBS -m ae" >> Novosort_${line}.pbs; done

# Populate PBS files (e-mail address for correspondence):
ls Trimmed_*.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "#PBS -M $email" >> Novosort_${line}.pbs; done

# Populate PBS files (empty line):
ls Trimmed_*.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "" >> Novosort_${line}.pbs; done

# Populate PBS files (load module(s)):
ls Trimmed_*.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "module load $novosort" >> Novosort_${line}.pbs; done

# Populate PBS files (empty line):
ls Trimmed_*.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "" >> Novosort_${line}.pbs; done

############################################
##        COMMAND LINE FOR NOVOSORT       ##
############################################
ls Trimmed_*.bam | cut -d"_" -f1-3 | cut -d"." -f1 | uniq | while read line; do echo "novosort -n -m ${novomem}G -c 6 $MainDirectory/SamtoolsNovosort/${line}.samtools.bam > $MainDirectory/SamtoolsNovosort/${line}.novosort.bam" >> Novosort_${line}.pbs; done

############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "==> Novosort sorting" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- Command Lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n1 Novosort_*.pbs >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "-- Job Submission IDs --" >> $MainDirectory/Pipeline_Log

############################################
##          SUBMIT NOVOSORT JOBS          ##
############################################
# Submit jobs and populate Log file (submission ID numbers):
ls Novosort_*.pbs | while read line; do echo "$line" >> $MainDirectory/Pipeline_Log ; qsub "$line"; done >> $MainDirectory/Pipeline_Log

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user"


############################################
##                   END                  ##
############################################

