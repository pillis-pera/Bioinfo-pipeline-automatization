Use this directory to run the pipeline locally in a BASH terminal

**DIRECTORIES**
- "r_*" Empty directories that will store the output of each process in the pipeline.
- "legacy" The original files used to create the pipeline (they are no longer needed).

**FIES**
- "Pipeline.sh" Contains the main miRNA processing pipeline. To skip processes of the pipeline or run a specific process, comment the other processes using "#". Remember to add the files needed for that process
- "merge_samples.sh" Script to merge two (or more?) files from different Flowcell lanes of the same sample. Unmerged samples will be taken from "Samples/unmerged" and the resulting merged files will be placed in "Samples/merged".

Note: Remember to give permissions to each file using "chmod +x filename.sh".
