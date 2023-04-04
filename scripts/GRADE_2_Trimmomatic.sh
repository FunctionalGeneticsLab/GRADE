
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

step="TRIMMING"

echo "This script will run STEP 2 of the GRADE pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Trim reads using Trimmomatic"

#  1 FASTQC PRE TRIM
#  2 TRIMMOMATIC <---
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
trimmomatic=trimmomatic/0.36


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
mkdir $MainDirectory/Trimmomatic


# 2 TRIMMING
############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Enter directory:
cd $MainDirectory

# Memory required (in GB):
memreq=80

# CPUs required (integer):
cpureq=12

############################################
##    CREATE PBS FILES FOR TRIMMOMATIC    ##
############################################
# PBS files for Trimmomatic are named Trimmomatic_LINE.pbs
# Cluster job names for Trimmomatic are named LINE_Trimmomatic

# List all fastq files and create correspondent PBS files:
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do cp $MainDirectory/Header.pbs Trimmomatic_${line}.pbs; done

# Populate PBS files (name of job):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "#PBS -N ${line}_Trimmomatic" >> Trimmomatic_${line}.pbs; done

# Populate PBS files (number of threads):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "#PBS -r n" >> Trimmomatic_${line}.pbs; done

# Populate PBS files (time and resources for running the job):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "#PBS -l mem=${memreq}GB,walltime=08:00:00,ncpus=${cpureq}" >> Trimmomatic_${line}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "#PBS -m ae" >> Trimmomatic_${line}.pbs; done

# Populate PBS files (e-mail address for correspondence):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "#PBS -M $email" >> Trimmomatic_${line}.pbs; done

# Populate PBS files (empty line):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "" >> Trimmomatic_${line}.pbs; done

# Populate PBS files (load module(s)):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "module load $trimmomatic" >> Trimmomatic_${line}.pbs; done

# Populate PBS files (empty line):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "" >> Trimmomatic_${line}.pbs; done

############################################
##       COMMAND LINE FOR TRIMMOMATIC     ##
############################################
# Populate PBS files (actual command line):
ls *.f* | cut -d"_" -f1-2 | uniq | while read line; do echo "trimmomatic PE -phred33 -threads $cpureq ${MainDirectory}/${line}_R1.fq.gz ${MainDirectory}/${line}_R2.fq.gz ${MainDirectory}/Trimmomatic/Trimmed_${line}_R1.fq.gz ${MainDirectory}/Trimmomatic/Unpaired_Trimmed_${line}_R1.fq.gz ${MainDirectory}/Trimmomatic/Trimmed_${line}_R2.fq.gz ${MainDirectory}/Trimmomatic/Unpaired_Trimmed_${line}_R2.fq.gz ILLUMINACLIP:/software/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:40" >> Trimmomatic_${line}.pbs ; done

# Other parameters to consider - HEADCROP:12 (the number of bases to remove from the start of the read)


############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "==> Trimmomatic Command Lines:" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- Command Lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n1 Trimmomatic_*.pbs >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "-- Job Submission IDs --" >> $MainDirectory/Pipeline_Log


############################################
##          SUBMIT TRIMMOMATIC JOBS       ##
############################################
# Submit jobs and populate Log file (submission ID numbers):
ls Trimmomatic_*.pbs | while read line; do echo "$line" >> $MainDirectory/Pipeline_Log ; qsub "$line"; done >> $MainDirectory/Pipeline_Log

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user"


############################################
##                   END                  ##
############################################

