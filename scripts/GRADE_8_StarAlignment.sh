
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

step="STAR ALIGNMENT"

echo "This script will run STEP 8 of the GRADE pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Align to transcriptome with Star"

#  1 FASTQC PRE TRIM
#  2 TRIMMOMATIC
#  3 FASTQC POST TRIM
#  4 KALLISTO INDEXING (if needed)
#  5 KALLISTO QUANTIFICATION
#  6 Create KALLISTO count tables
#  7 STAR INDEXING (if needed)
#  8 STAR ALIGNMENT <---
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
star=STAR/2.7.1a

############################################
##          REFERENCE  DIRECTORY          ##
############################################
RefDirectory=/working/lab_julietF/mainaB/ReferenceGenomes

############################################
##           SET REFERENCE VARs           ##
############################################
# Genome sequence and transcriptome annotation files for STAR.

# ENSEMBL:
# Genome build GRCh38.p2, updated in 2015
# Supplemented with 92 ERCC entries.
# 65309 annotated genes and 198370 transcripts.

if [ "$1" == ensembl ] || [ "$1" == E ] ; then
transcriptomeseq=${RefDirectory}/TranscriptomeGRCh38p7_ERCC.fa
transcriptomegtf=${RefDirectory}/TranscriptomeGRCh38rel79_ERCC.gtf
GenDirectory=${RefDirectory}/StarIndex/Ensembl
fi

# GENCODE:

if [ "$1" == gencode ] || [ "$1" == G ] ; then
transcriptomeseq=${RefDirectory}/gencode.v36.transcripts.fa
transcriptomegtf=${RefDirectory}/gencode.v36.annotation.gtf
GenDirectory=${RefDirectory}/StarIndex/Gencode
fi


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

# Create Star directory:
mkdir $MainDirectory/StarAlignment


# 8 STAR ALIGNMENT
############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Enter directory:
cd $MainDirectory/Trimmomatic

# Memory required (in GB):
memreq=40

# CPUs required (integer):
cpureq=12


############################################
##       CREATE PBS FILES FOR STAR        ##
############################################
# List all fastq files and create correspondent PBS files:
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do cp $MainDirectory/Header.pbs StarAlignment_${line}.pbs; done

# Populate PBS files (name of job):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "#PBS -N ${line}_Star" >> StarAlignment_${line}.pbs; done

# Populate PBS files (number of threads):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "#PBS -r n" >> StarAlignment_${line}.pbs; done

# Populate PBS files (time and resources for running the job):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "#PBS -l mem=${memreq}GB,walltime=08:00:00,ncpus=${cpureq}" >> StarAlignment_${line}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "#PBS -m ae" >> StarAlignment_${line}.pbs; done

# Populate PBS files (e-mail address for correspondence):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "#PBS -M $email" >> StarAlignment_${line}.pbs; done

# Populate PBS files (empty line):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "" >> StarAlignment_${line}.pbs; done

# Populate PBS files (load module(s)):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "module load $star" >> StarAlignment_${line}.pbs; done

# Populate PBS files (empty line):
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do echo "" >> StarAlignment_${line}.pbs; done


############################################
##          COMMAND LINE FOR STAR         ##
############################################
ls Trimmed_*.f* | cut -d"_" -f1-3 | uniq | while read line; do mkdir $MainDirectory/StarAlignment/$line; echo "STAR --runMode alignReads --readFilesIn $MainDirectory/Trimmomatic/${line}_R1.fq.gz $MainDirectory/Trimmomatic/${line}_R2.fq.gz --readFilesCommand zcat --outFileNamePrefix $MainDirectory/StarAlignment/$line/$line --genomeDir $GenDirectory --sjdbGTFfile $transcriptomegtf --outSJfilterReads Unique --sjdbOverhang 100 --runThreadN ${cpureq} --genomeLoad NoSharedMemory --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 --chimSegmentMin 20 --outSAMattributes All --outSAMstrandField intronMotif --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate" >> StarAlignment_${line}.pbs; done


############################################
##          STAR COMMAND EXPLAINED        ##
############################################

# genomeDir: path to the directory where genome files are stored (if runMode!=generateGenome) or will be generated (if runMode==generateGenome)

#		--readFilesCommand zcat \			string(s): command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout
#		--sjdbGTFfile ${gtfDir} \			path to the GTF file with annotations
#		--outSJfilterReads Unique \			string: which reads to consider for collapsed splice junctions output (uniquely mapping reads only)
#		--sjdbOverhang 114 \				int>0: length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1)
#		--twopassMode Basic \				string: 2-pass mapping mode (basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly)
#		--outFilterType BySJout \			string: type of filtering (BySJout ... keep only those reads that contain junctions that passed filtering into SJ.out.tab)
#		--outFilterMultimapNmax 100 \			int: maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as "mapped to too many loci" in the Log.final.out .
#		--outFilterMismatchNmax 33 \			int: alignment will be output only if it has no more mismatches than this value.
#		--outFilterMatchNminOverLread 0 \		float: same as outFilterMatchNmin, but normalized to the read length (sum of mates' lengths for paired-end reads). // int: alignment will be output only if the number of matched bases is higher than or equal to this value.
#		--outFilterMismatchNoverLmax 0.3 \		float: alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value.
#		--outFilterScoreMinOverLread 0.3 \		float: same as outFilterScoreMin, but  normalized to read length (sum of mates' lengths for paired-end reads) // int: alignment will be output only if its score is higher than or equal to this value.
#		--limitOutSJcollapsed 1000000 \			int>0: max number of collapsed junctions
#		--limitSjdbInsertNsj 1000000 \			int>=0: maximum number of junction to be inserted to the genome on the fly at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run
#		--alignSJoverhangMin 8 \			int>0: minimum overhang (i.e. block size) for spliced alignments
#		--alignEndsType EndToEnd \			string: type of read ends alignment (force end-to-end read alignment, do not soft-clip)
#		--alignSJDBoverhangMin 3  \			int>0: minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments
#		--alignIntronMin 20 \				maximum intron size (if 0, max intron size will be determined by (2^winBinNbits)*winAnchorDistNbins)
#		--winAnchorMultimapNmax 50 \			int>0: max number of loci anchors are allowed to map to
#		--seedSearchStartLmax 12 \			int>0: defines the search start point through the read - the read is split into pieces no longer than this value
#		--chimSegmentMin 20 \				int>=0: minimum length of chimeric segment length (if ==0, no chimeric output)
#		--outSAMattributes All \			string: a string of desired SAM attributes, in the order desired for the output SAM (All = NH HI AS nM NM MD jM jI)
#		--outSAMstrandField intronMotif \		string: Cufflinks-like strand field flag (intronMotif = strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.)
#		--quantMode TranscriptomeSAM \			string(s): types of quantification requested (TranscriptomeSAM = output SAM/BAM alignments to transcriptome into a separate file)
#		--outSAMattrIHstart 0 \				int>=0: start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie.
#		--outSAMunmapped Within \			string(s): output of unmapped reads in the SAM format (Within = output unmapped reads within the main SAM file (i.e. Aligned.out.sam))

############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "==> Star Alignment" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- Command Lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n1 StarAlignment_*.pbs >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "-- Job Submission IDs --" >> $MainDirectory/Pipeline_Log


############################################
##             SUBMIT STAR JOBS           ##
############################################
# Submit jobs and populate Log file (submission ID numbers):
ls StarAlignment_*.pbs | while read line; do echo "$line" >> $MainDirectory/Pipeline_Log ; qsub "$line"; done >> $MainDirectory/Pipeline_Log

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user"


############################################
##                   END                  ##
############################################


