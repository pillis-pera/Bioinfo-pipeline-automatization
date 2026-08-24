This repo contains a bash script that lets you process raw reads obtained from an Illumina S2 cell NGS run prepared with a QIAGEN miRNA library prep kit with UMIs.

FILES
- "FromFastqToBamBai.sh" is a script that contains just the alignment part of the full pipeline.
- "Sample to Results.sh" is script that contains the full pipeline (from merging L1 & L2 FASTQs to creating the .csv with deduplicated counts).
- "hsa_mature_DNA.fa" is the collection of mature Homo Sapiens miRNAs in FASTA using DNA annotation retrieved from miRBase.



CONSIDERATIONS
- This pipeline process all FASTQ samples inside the directory, keeping their base name to generate all subsequent file types.
- Each step puts the output inside a new directory.
- Samples are processed one after the other (I'm working on parallel processing implementation using HPC).
- Only "hsa_mature_DNA.fa" is needed to run the pipeline, since the build step is already included in the "Sample to Results.sh"
