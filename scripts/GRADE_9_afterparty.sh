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
##         CHECK STEP CONCLUSION!         ##
############################################

############################################
##         SET MAIN WORK DIRECTORY        ##
############################################
# Define directory:

if [ "$#" -gt 0 ]; then MainDirectory=$1; fi
if [ "$#" -eq 0 ]; then MainDirectory=`pwd`; fi

echo ""; echo "WARNING: The main directory for this run was set to ${MainDirectory}."; echo ""

############################################
##          MOVE FILES TO OUTPUT          ##
############################################
cd $MainDirectory/StarAlignment
mv Samtools_*.pbs $MainDirectory/PBSin
mv *Samtools*.e* $MainDirectory/PBSout
mv *Samtools*.o* $MainDirectory/PBSout

############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""
