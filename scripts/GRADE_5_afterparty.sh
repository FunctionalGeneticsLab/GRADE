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
cd $MainDirectory/Trimmomatic
mv KallistoQuant*.pbs $MainDirectory/PBSin
mv *KallistoQuant*.e* $MainDirectory/PBSout
mv *KallistoQuant*.o* $MainDirectory/PBSout

############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""
