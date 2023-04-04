
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

step="KALLISTO QUANTIFICATION"

echo "This script will run STEP 5 of the GRADE pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Quantify pseudo-counts with Kallisto"

#  1 FASTQC PRE TRIM
#  2 TRIMMOMATIC
#  3 FASTQC POST TRIM
#  4 KALLISTO INDEXING (if needed)
#  5 KALLISTO QUANTIFICATION <---
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
kallisto=kallisto/0.43.0

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

# Create Kallisto directory:
mkdir $MainDirectory/KallistoQuant


# 5 KALLISTO QUANTIFICATION
############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Enter directory:
cd $MainDirectory/Trimmomatic

# Memory required (in GB):
memreq=25

# CPUs required (integer):
cpureq=12


############################################
##     CREATE PBS FILES FOR KALLISTO      ##
############################################
# List all fastq files and create correspondent PBS files:
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do cp $MainDirectory/Header.pbs KallistoQuant_${line}.pbs; done

# Populate PBS files (name of job):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "#PBS -N ${line}_KallistoQuant" >> KallistoQuant_${line}.pbs; done

# Populate PBS files (number of threads):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "#PBS -r n" >> KallistoQuant_${line}.pbs; done

# Populate PBS files (time and resources for running the job):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "#PBS -l mem=${memreq}GB,walltime=08:00:00,ncpus=${cpureq}" >> KallistoQuant_${line}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "#PBS -m ae" >> KallistoQuant_${line}.pbs; done

# Populate PBS files (e-mail address for correspondence):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "#PBS -M $email" >> KallistoQuant_${line}.pbs; done

# Populate PBS files (empty line):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "" >> KallistoQuant_${line}.pbs; done

# Populate PBS files (load module(s)):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "module load $kallisto" >> KallistoQuant_${line}.pbs; done

# Populate PBS files (empty line):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "" >> KallistoQuant_${line}.pbs; done

############################################
##        COMMAND LINE FOR KALLISTO       ##
############################################

# ENSEMBL:
if [ "$1" == ensembl ] || [ "$1" == E ] ; then
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "kallisto quant -t $cpureq -i $RefDirectory/KallistoIndex/EnsemblTranscriptome -o $MainDirectory/KallistoQuant/KallistoQuant_${line} -b 100 <(zcat $MainDirectory/Trimmomatic/${line}_R1.fq.gz) <(zcat $MainDirectory/Trimmomatic/${line}_R2.fq.gz)" >> KallistoQuant_${line}.pbs; done
fi

# GENCODE:
if [ "$1" == gencode ] || [ "$1" == G ] ; then
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "kallisto quant -t $cpureq -i $RefDirectory/KallistoIndex/GencodeTranscriptome -o $MainDirectory/KallistoQuant/KallistoQuant_${line} -b 100 <(zcat $MainDirectory/Trimmomatic/${line}_R1.fq.gz) <(zcat $MainDirectory/Trimmomatic/${line}_R2.fq.gz)" >> KallistoQuant_${line}.pbs; done
fi

# Other parameters to consider -l (estimated average fragment length) and -s (estimated standard deviation of fragment length) --fr-stranded (strand specific reads, first read forward)

############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "==> Kallisto Quantification" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- Command Lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n1 KallistoQuant_*.pbs >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "-- Job Submission IDs --" >> $MainDirectory/Pipeline_Log


############################################
##           SUBMIT KALLISTO JOBS         ##
############################################
# Submit jobs and populate Log file (submission ID numbers):
ls KallistoQuant_*.pbs | while read line; do echo "$line" >> $MainDirectory/Pipeline_Log ; qsub "$line"; done >> $MainDirectory/Pipeline_Log

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user"


############################################
##                   END                  ##
############################################

