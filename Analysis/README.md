**Follow these steps to analyze count vectors**

1 - Move/Copy all count vectors into Analysis/Counts.

2 - Load your metadata.csv file to the Analysis/Metadata directory

3 - Run Analysis/post_processing.ipynb to generate the counts matrix and filter low expression/inconsistent miRNAs.

4 - With _counts_matrix_filtered.csv_ in the directory, run differential_expression.R and volcano_plot.R in an R environment 
(environment must contain the necesary libraries form the **# Import libraries #** section)
(rename metadata variables to match your metadata.csv columns)
