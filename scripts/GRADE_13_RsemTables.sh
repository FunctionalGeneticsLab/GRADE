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
##         MAKE INTEGRATED TABLES         ##
############################################

cd $MainDirectory/RsemQuant

# Expected Counts:

ls *genes.results | cut -d"." -f1 | while read line; do echo "$line"; sample=`echo "$line" | cut -d"_" -f2-3 | cut -c1-14`; echo -e "Gene\t$sample" >> RsemQuant_${line}_ExpCounts; cut -f1,5 ${line}.genes.results | grep "ENST" >> RsemQuant_${line}_ExpCounts; done

first=`ls RsemQuant_*_ExpCounts | head -n2 | head -n1`; second=`ls RsemQuant_*_ExpCounts | head -n2 | tail -n1`; n=`ls RsemQuant_*_ExpCounts | wc -l`; ns=`echo "${n}-2" | bc`

join -t$'\t' ${first} ${second} >> Seed_ExpCounts

ls RsemQuant_*_ExpCounts | tail -n${ns} | while read line; do mv Seed_ExpCounts Active_ExpCounts; join -t$'\t' Active_ExpCounts ${line} >> Seed_ExpCounts; done

mv Seed_ExpCounts Rsem_All_ExpCounts; rm -rf Active_ExpCounts

echo -e "\n==> Your ExpCounts table is ready!\n Here is a sneak peak:"
head  $MainDirectory/RsemQuant/Rsem_All_ExpCounts



echo ""; echo "---> FINISHED"; echo ""
