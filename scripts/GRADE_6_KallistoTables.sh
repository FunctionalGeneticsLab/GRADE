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

cd $MainDirectory/KallistoQuant

# Estimated Counts:

ls -d */ | while read line; do seed=`echo $line | cut -d'_' -f3-` ; cut -f1,4 ${line}abundance.tsv | sed 's/target_id/Transcript/g' | sed "s|est_counts|${seed}|g" | tr -d '/' >> ${line}EstCounts; done

first=`ls -d */ | head -n2 | head -n1`; second=`ls -d */ | head -n2 | tail -n1`; n=`ls -d */ | wc -l`; ns=`echo "${n}-2" | bc`; np=`echo "${n}+1" | bc`

join -t$'\t' ${first}EstCounts ${second}EstCounts >> Seed_EstCounts

ls -d */ | tail -n${ns} | while read line; do mv Seed_EstCounts Active_EstCounts; join -t$'\t' Active_EstCounts ${line}EstCounts >> Seed_EstCounts; done

cat Seed_EstCounts | sed 's/Transcript/FullID/g' >> Kallisto_All_EstCounts_FullIDs
cut -d'|' -f1 Seed_EstCounts >> temp_Col1_EstCounts
cut -f2- Seed_EstCounts >> temp_Col2_EstCounts
paste temp_Col1_EstCounts temp_Col2_EstCounts | cut -f1-${np} >> Kallisto_All_EstCounts
rm -rf Active_EstCounts Seed_EstCounts temp_Col*_EstCounts

echo -e "\n==> Your Estimated Counts table is ready!\n Here is a sneak peak:"
head  $MainDirectory/KallistoQuant/Kallisto_All_EstCounts


# TPM:

ls -d */ | while read line; do seed=`echo $line | cut -d'_' -f3-` ; cut -f1,5 ${line}abundance.tsv | sed 's/target_id/Transcript/g' | sed "s|tpm|${seed}|g" | tr -d '/' >> ${line}TPM; done

first=`ls -d */ | head -n2 | head -n1`; second=`ls -d */ | head -n2 | tail -n1`; n=`ls -d */ | wc -l`; ns=`echo "${n}-2" | bc`

join -t$'\t' ${first}TPM ${second}TPM >> Seed_TPM

ls -d */ | tail -n${ns} | while read line; do mv Seed_TPM Active_TPM; join -t$'\t' Active_TPM ${line}TPM >> Seed_TPM; done

cat Seed_TPM | sed 's/Transcript/FullID/g' >> Kallisto_All_TPM_FullIDs 
cut -d'|' -f1 Seed_TPM >> temp_Col1_TPM
cut -f2- Seed_TPM >> temp_Col2_TPM
paste temp_Col1_TPM temp_Col2_TPM | cut -f1-${np} >> Kallisto_All_TPM
rm -rf Active_TPM Seed_TPM temp_Col*_TPM


echo -e "\n==> Your TPM counts table is ready!\n Here is a sneak peak:"
head  $MainDirectory/KallistoQuant/Kallisto_All_TPM
echo ""

############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""
