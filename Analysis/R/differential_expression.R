setwd("/path/to/directory/Bioinfo-pipeline-automatization/Analysis") #<-- Set this to your wd.
library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2)
counts <- read.csv('counts_matrix_filtered.csv', row.names=1)
coldata <- read.csv('Metadata/metadata.csv', row.names = 1)

coldata <- coldata[colnames(counts),]

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~Sex + Condition)
dds <- DESeq(dds)
vst <- varianceStabilizingTransformation(dds)   
NCounts <- counts(dds, normalized = T)          #Normalized counts
log2_NCounts <- log2(NCounts+1)                 #Log2 of normalized counts

res <- results(dds)
resOrdered <- res[order(res$padj), ]
sig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) >= log2(1.5)), ]
sig <- sig[order(sig$padj),]

miRNAs <- row.names(sig[grep("miR|let",row.names(sig)),])
miRNAs_res <- sig[grep("miR|let",row.names(sig)),]