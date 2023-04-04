
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
mkdir $MainDirectory/PBSin
mkdir $MainDirectory/PBSout

############################################
##           CREATE HEADER FILE           ##
############################################
# Populate header file:
echo "###########################################################################################" >> $MainDirectory/Header.pbs
echo "#" >> $MainDirectory/Header.pbs
echo "#  Script:    Script for GRADE (General RNAseq Analysis for Differential Expression) in PBS" >> $MainDirectory/Header.pbs
echo "#  Author:    Maina Bitar" >> $MainDirectory/Header.pbs
echo "#  Created:   2016 at QIMR Berghofer (Brisbane, Australia)" >> $MainDirectory/Header.pbs
echo "#  Updated:   2023 at QIMR Berghofer (Brisbane, Australia)" >> $MainDirectory/Header.pbs
echo "#  Email:     Maina.Bitar@qimrberghofer.edu.au" >> $MainDirectory/Header.pbs
echo "#" >> $MainDirectory/Header.pbs
echo "###########################################################################################" >> $MainDirectory/Header.pbs


############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""
