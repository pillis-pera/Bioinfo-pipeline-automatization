This repo contains a pipeline to process raw reads obtained from an Illumina S2 cell NGS run prepared with a QIAGEN miRNA library prep kit with UMIs.

**DIRECTORIES**
- "Local" To run the pipeline locally on a BASH terminal.
- "NextFlow" To use Nextflow to run the pipeline (Locally or in a HPC).
- "SLURM" To run the pipeline using SLURM in an HPC context.
- "Reference" Contains the _human miRNome_ to use as a reference when mapping reads.
- "Samples" Samples to be processed through the pipeline.
- "Analysis" To convert vector counts resulting from running the pipeline into a counts matrix and perform biostatistical analyses with R.

Note: The pipeline needs to be run in an environment containing _umi_tools_, _samtools_, _bowtie_ & _fastp_.
