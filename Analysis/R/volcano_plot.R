####################
# Import libraries #
####################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(here)}
)

#################
# Configuration #
#################
padj_threshold <- 0.05                  # Set adjusted p-value threshold
fc_threshold <- 1.5                     # Set foldchange threshold
lfc_threshold <- log2(fc_threshold)     # Set log2(foldchange) threshold
test <- "standard"                      # Options: standard or threshold
label_top_n_results <- 10               # Set how many miRNAs should be labeled by name


##############
# Load files #
##############
DE_RESULTS <- here("Analysis",
                   "DESeq_results", 
                   paste0("deseq_results_", test, ".csv")
)
OUTPUT_DIR <- here("Analysis",
                      "DESeq_results"
)

#############
# Load data #
#############
de_results <- read.csv(DE_RESULTS, 
                       row.names = 1
                       )

################
# Prepare data #
################

#Remove NAs
de_results <- de_results[complete.cases(de_results), , drop = FALSE] 

#Classify results
de_results$significance <- "Not significant"

de_results$significance[
  de_results$padj < padj_threshold
  & 
  de_results$log2FoldChange > lfc_threshold  
] <- "Up-regulated"

de_results$significance[
  de_results$padj < padj_threshold
  & 
    de_results$log2FoldChange < -lfc_threshold  
] <- "Down-regulated"

#labels
label_mirnas <- de_results[1:label_top_n_results,]

################
# Volcano plot #
################
plot <- ggplot(de_results, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = c("Up-regulated" = "firebrick", "Down-regulated" = "cornflowerblue", "Not significant" = "gray75")) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed") +
  geom_label_repel(
    data = label_mirnas,
    aes(label = row.names(label_mirnas)),
    size = 4,
    color = "black",
    fill = NA,
    box.padding = 0.3,
    segment.color = "black",
    segment.size = 0.6,
    segment.alpha = 1,
    min.segment.length = 0,
    point.padding = 0.3,
    label.size = 0.3,
    max.overlaps = Inf,
    max.iter = 10000000,
    max.time = 5,
    seed = 21
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black")) +
  coord_cartesian(xlim = c(min(de_results$log2FoldChange)-0.5, max(de_results$log2FoldChange)+0.5), 
        ylim = c(0, max(-log10(de_results$padj))+1)
    )

#############
# Save plot #
#############
ggsave(
  filename = file.path(OUTPUT_DIR, paste0("volcano_plot_", test, ".png")),
  plot = plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)