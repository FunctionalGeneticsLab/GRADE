

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
##          REFERENCE  DIRECTORY          ##
############################################
RefDirectory=/working/lab_julietF/mainaB/ReferenceGenomes

############################################
##         SET MAIN WORK DIRECTORY        ##
############################################
# Define directory:

if [ "$#" -gt 1 ]; then MainDirectory=$2; fi
if [ "$#" -eq 1 ]; then MainDirectory=`pwd`; fi

echo ""; echo "WARNING: The main directory for this run was set to ${MainDirectory}."; echo ""

############################################
##         MAKE GENE-LEVEL TABLES         ##
############################################

## GENCODE
if [ "$1" == gencode ] || [ "$1" == G ] ; then
t2g=${RefDirectory}/gencode.v36.gene2transcript
fi

## ENSEMBL
if [ "$1" == ensembl ] || [ "$1" == E ] ; then
t2g=${RefDirectory}/ensembl.GRCh38rel79.gene2transcript
fi


## KALLISTO
cd $MainDirectory/KallistoQuant

# Estimated Counts (Kallisto_All_EstCounts):
sort -k1,1 -r Kallisto_All_EstCounts >> TempFile
join -t$'\t' ${t2g} TempFile | cut -f2- | sort -k1,1 -r | awk 'NR==1{print;next} {for (i=2;i<=NF;i++) {a[$1][i]+=$i}}END{for (j in a) {s=j; for (i=2;i<=NF;i++) {s=s"\t"a[j][i]}; print s}}' >> Kallisto_All_GeneLevel_EstCounts
rm -rf TempFile

# TPM (Kallisto_All_TPM):
sort -k1,1 -r Kallisto_All_TPM >> TempFile
join -t$'\t' ${t2g} TempFile | cut -f2- | sort -k1,1 -r | awk 'NR==1{print;next} {for (i=2;i<=NF;i++) {a[$1][i]+=$i}}END{for (j in a) {s=j; for (i=2;i<=NF;i++) {s=s"\t"a[j][i]}; print s}}' >> Kallisto_All_GeneLevel_TPM
rm -rf TempFile

echo -e "\n==> Your Kallisto Gene-Level counts tables are ready!\n Here is a sneak peak:"
head  $MainDirectory/KallistoQuant/Kallisto_All_GeneLevel_*
echo ""



## RSEM
cd $MainDirectory/RsemQuant

# Estimated Counts (Rsem_All_ExpCounts):
sort -k1,1 -r Rsem_All_ExpCounts >> TempFile
join -t$'\t' ${t2g} TempFile | cut -f2- | sort -k1,1 -r | awk 'NR==1{print;next} {for (i=2;i<=NF;i++) {a[$1][i]+=$i}}END{for (j in a) {s=j; for (i=2;i<=NF;i++) {s=s"\t"a[j][i]}; print s}}' >> Rsem_All_GeneLevel_ExpCounts
rm -rf TempFile

# TPM (Rsem_All_TPM):
sort -k1,1 -r Rsem_All_TPM >> TempFile
join -t$'\t' ${t2g} TempFile | cut -f2- | sort -k1,1 -r | awk 'NR==1{print;next} {for (i=2;i<=NF;i++) {a[$1][i]+=$i}}END{for (j in a) {s=j; for (i=2;i<=NF;i++) {s=s"\t"a[j][i]}; print s}}' >> Rsem_All_GeneLevel_TPM
rm -rf TempFile

# FPKM (Rsem_All_FPKM):
sort -k1,1 -r Rsem_All_FPKM >> TempFile
join -t$'\t' ${t2g} TempFile | cut -f2- | sort -k1,1 -r | awk 'NR==1{print;next} {for (i=2;i<=NF;i++) {a[$1][i]+=$i}}END{for (j in a) {s=j; for (i=2;i<=NF;i++) {s=s"\t"a[j][i]}; print s}}' >> Rsem_All_GeneLevel_FPKM
rm -rf TempFile

echo -e "\n==> Your Rsem Gene-Level counts tables are ready!\n Here is a sneak peak:"
head  $MainDirectory/RsemQuant/Rsem_All_GeneLevel_*
echo ""



############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""
