setwd("/path/to/Bioinfo-pipeline-automatization") #<-- Set this to your wd.

####################
# Import libraries #
####################
suppressPackageStartupMessages({
    library(dplyr)
    library(DESeq2)}
)

#################
# Configuration #
#################
padj_threshold <- 0.05
fc_threshold <- 1.5
lfc_threshold <- log2(fc_threshold)


#############
# Load data #
#############
counts_matrix <- read.csv('Analysis/counts_matrix_filtered.csv', row.names=1, check.names=FALSE)
sample_metadata <- read.csv('Analysis/Metadata/metadata.csv', row.names = 1)
sample_metadata <- sample_metadata[colnames(counts_matrix),]

####################
# Prepare metadata #
####################
sample_metadata$sex <- factor(sample_metadata$sex)
sample_metadata$condition <- factor(sample_metadata$condition)
sample_metadata$condition <- relevel(
  factor(sample_metadata$condition),
  ref = "C"
)

################
# DESeq2 model #
################
dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix, 
    colData = sample_metadata, 
    design = ~sex + condition
)
dds <- DESeq(dds)

###################################
# Normalized / transformed counts #
###################################
vst_counts <- varianceStabilizingTransformation(dds)   

normalized_counts <- counts(
    dds, 
    normalized = TRUE
)          #Normalized counts

log2_normalized_counts <- log2(normalized_counts+1)                 #Log2 of normalized counts

###########
# Results #
###########
results_deseq <- results(
    dds,
    alpha = padj_threshold,
    lfcThreshold = lfc_threshold,
    altHypothesis = "greaterAbs"
)

results_deseq <- results_deseq[order(results_deseq$padj), ]

#########################
# Significant smallRNAs #
#########################
significant_results <- results_deseq[
    which(
        results_deseq$padj < padj_threshold 
        & 
        abs(results_deseq$log2FoldChange) >= lfc_threshold
    ), 
]

######################
# Significant miRNAs #
######################
significant_mirnas <- row.names(
    significant_results[
        grep("miR|let",row.names(significant_results)
        ),
    ]
)