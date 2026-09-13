####################
# Import libraries #
####################
suppressPackageStartupMessages({
    library(DESeq2)
    library(here)}
)

##############
# Load files #
##############
COUNTS_FILE <- here("Analysis","counts_matrix_filtered.csv")
METADATA_FILE <- here("Analysis","Metadata","metadata.csv")
OUTPUT_DIR <- here("Analysis","DESeq_results")

#################
# Configuration #
#################
padj_threshold <- 0.05                  # Set adjusted p-value threshold
fc_threshold <- 1.5                     # Set foldchange threshold
lfc_threshold <- log2(fc_threshold)     # Set log2(foldchange) threshold

test <- "standard"                      # Options: standard or threshold

                                        # 'standard'    Ho: |log2FC| = 0 
                                        #               Ha: |log2FC| != 0
                                        # 'threshold'   Ho: |log2FC| <= lfc_threshold
                                        #               Ha: |log2FC| > lfc_threshold

#############
# Load data #
#############
counts_matrix <- read.csv(COUNTS_FILE, 
    row.names=1, 
    check.names=FALSE
)

sample_metadata <- read.csv(METADATA_FILE, 
    row.names = 1
)
sample_metadata <- sample_metadata[colnames(counts_matrix),]

######################
# Create results dir #
######################
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

####################
# Prepare metadata #
####################
sample_metadata$sex <- factor(sample_metadata$sex)
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

if (test == "threshold") {
    results_deseq <- results(
        dds,
        alpha = padj_threshold,
        lfcThreshold = lfc_threshold,
        altHypothesis = "greaterAbs"
    )
} else if (test == "standard") {
    results_deseq <- results(
    dds,
    alpha = padj_threshold
    )
} else {
    stop("Invalid test: '", test, "'. Use 'standard' or 'threshold'.")
}
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

################
# Save results #
################
write.csv(
    results_deseq,
    file.path(
        OUTPUT_DIR, 
        paste0("deseq_results_", test, ".csv")
    )
)
write.csv(
    log2_normalized_counts,
    file.path(OUTPUT_DIR, "log2_normalized_reads.csv")
)
write.csv(
    assay(vst_counts),
    file.path(OUTPUT_DIR, "vst_counts.csv")
)