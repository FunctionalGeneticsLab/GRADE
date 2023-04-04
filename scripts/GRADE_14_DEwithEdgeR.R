
#!/usr/bin/env Rscript

#module load R/4.0.2

# Use arguments:
args<-commandArgs(TRUE)

# args[1] = input table of counts
# args[2] = minimum CPM to consider expression
# args[3] = minimum samples to contain gene so it is considered for analysis
# args[4] = number of genes or transcripts in the table
# args[5] = sample IDs separated by commas (e.g. 1,1,1,2,2,2)
# args[6] = Gene or Transcript
# args[7] = Output files seed name


# This is to allow printing via sink to write as many lines as needed to an output file.
options(max.print = .Machine$integer.max)

# This is to load the edgeR library
library(edgeR)

# This is to read a tab-delimited table of raw counts, which is the edgeR input.
# |--> Provide input file:
inputfile <- args[1]
x <- read.delim(inputfile, row.names=args[6])

# This is to create groups and define replicates within the samples.
# |--> Provide groups:
vector <- as.numeric(unlist(strsplit(args[5],",")))
group_v <- vector
group <- factor(group_v)

# Example of group:
#group_v <- c(1,1,1,2,2,2,3,3,3)

# This is to initialize the object that will store information on DE genes. 
y <- DGEList
y <- DGEList(counts=x,group=group)

# This is to print the library size of each sample and help define the threshold for keeping reads.
# A reasonable cutoff is in between 5 and 10 counts, which are frequently around 1 CPM.
y$samples

# This is to discard reads with too low counts.
# A requirement for expression in 2 or more libraries is used when the minimum number of samples in each group is 3.
# |--> Provide cutoffs:
keep <- rowSums(cpm(y)>args[2]) >= args[3]
y <- y[keep, , keep.lib.sizes=FALSE]

# This is to normalize for RNA composition using a trimmed mean of M-values (TMM) between each pair of samples.
y <- calcNormFactors(y)

# This is to create the design matrix to be used for dispersion estimation.
design <- model.matrix(~group)

# Use the Cox-Reid profile-adjusted likelihood (CR) method
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

# -- Alternatively:
#y <- estimateDisp(y, design)

# This is to fit the negative binomial GLM for each tag.
fit <- glmQLFit(y, design)

# This is for dispersion estimation and hypothesis testing.
# |--> CHANGE coefficients if pairwise analysis is desired:
qlf <- glmQLFTest(fit, coef=2:3)

# This is to get the number of lines (transcripts) in qlf.
nr <- length(rownames(qlf))

# This is to generate FDRs for each transcript.
qlftt <- topTags(qlf,adjust.method="BH",n=nr)

# The topTags is usually used to extract the top DE tags ranked by p-value or absolute log-fold change.
# We used it to add the FDR column to the statistics table.

# This is to save the generated data.
qlft <- qlftt$table

# This is to get the number of columns in qlftt.
i <- ncol(qlft)

# This is to get the calculated p-values and adjust them using the Benjamini Hochberg method for multiple comparison adjustment.
# |--> Provide number of genes/transcripts:
p <- qlft[,i-1]
padjust <- p.adjust(p, method = "BH", n = as.numeric(args[4]))
# For ENSEMBL, n=X or n=198278, respectively the number of annotated genes and transcripts.
# For GENCODE, n=60660 or n=232117, respectively the number of annotated genes and transcripts.

# This is to add the column of adjusted p-values to the main table.
i <- i+1
qlft[,i] <- padjust
colnames(qlft)[i] <- "AdjPValue"

# This is to transform the data frame into a table.
qlfo <- cbind(rownames(qlft), qlft)
rownames(qlfo) <- rownames(qlft)
colnames(qlfo) <- c(args[6], c(colnames(qlft)))

# This is to automatically get column names for logFCs and other statistics separately.
i <- ncol(qlf$table)
cnst <- colnames(qlfo)[c(seq(i-3, i+1))]
cnfc <- colnames(qlfo)[c(seq(2, 3, i-4))]

# |--> CHANGE output tables names, if desired:
# This is to write these into tables.
write.table(qlfo[,c(args[6], c(cnfc))],file=sprintf("FoldChanges_%s.table",args[7]),quote=FALSE,row.names=FALSE,sep="\t",eol="\n")
write.table(qlfo[,c(args[6], c(cnst))],file=sprintf("Statistics_%s.table",args[7]),quote=FALSE,row.names=FALSE,sep="\t",eol="\n")
write.table(qlfo,file=sprintf("FullTable_%s.table",args[7]),quote=FALSE,row.names=FALSE,sep="\t",eol="\n")
q()
