This repo contains a pipeline process raw reads obtained from an Illumina S2 cell NGS run prepared with a QIAGEN miRNA library prep kit with UMIs.

DIRECTORIES
- "Local" To run the pipeline locally con a BASH terminal.
- "NextFlow" To use Nextflow to run the pipeline (Locally or in a HPC).
- "SLURM" To run the pipeline using SLURM in an HPC context.
- "Reference" Contains the _human miRNome_ to use as a reference when mapping reads.
- "Samples" Samples to be processed through the pipeline.
