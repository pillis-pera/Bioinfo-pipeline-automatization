#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
UNMERGED_DIR="${REPO_DIR}/Samples/unmerged"
MERGED_DIR="${REPO_DIR}/Samples/merged"
UMI_EXTRACTED_DIR="${REPO_DIR}/Local/r_umi_extracted"
ALIGNED_DIR="${REPO_DIR}/Local/r_aligned"
SORTED_DIR="${REPO_DIR}/Local/r_sorted"
DEDUPED_DIR="${REPO_DIR}/Local/r_deduped"
COUNTS_DIR="${REPO_DIR}/Local/r_counts"
REFERENCE_DIR="${REPO_DIR}/Reference"

# echo $UNMERGED_DIR
# echo $MERGED_DIR
# echo $UMI_EXTRACTED_DIR
# echo $ALIGNED_DIR
# echo $SORTED_DIR
# echo $DEDUPED_DIR
# echo $COUNTS_DIR

# cd $MERGED_DIR
# for file in *.fastq*; do
# 	base="${file%.*}"
# 	echo "Processing file: $file"
#     # echo $base
# 	umi_tools extract --extract-method=regex --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.*)$' -I "$file" -S "${UMI_EXTRACTED_DIR}/${base}.umi_extracted.fastq.gz"
# 	echo "Completed file: $file"
# 	echo ""
# done

cd $UMI_EXTRACTED_DIR
# bowtie-build "${REFERENCE_DIR}/hsa_mature_DNA.fa" "${REFERENCE_DIR}/hsa_mature_index"
> ${ALIGNED_DIR}/mapping_logs/mapping_stats.txt
for file in *.fastq*; do
	base="${file%.umi_extracted.fastq*}"
	#creates a variable with the name of the file without the extention to name all subsequent files created from this .fastq
	
	echo "Sample: $base" >> "${ALIGNED_DIR}/mapping_logs/mapping_stats.txt"
	#writes the name of the file (without extension) to a new line in the .txt file
	
	bowtie -q -S "${REFERENCE_DIR}/hsa_mature_index" "$file" > "${ALIGNED_DIR}/${base}.sam" 2>> "${ALIGNED_DIR}/mapping_logs/temporal_bowtie_output.txt"
	# -q for .fastq files (change it to the corresponding option accordingly, check bowtie), this line writes the termianls output into a temporal .txt file for later grepping and creates a .sam file from the mapping
	
	result=$(grep '^#' "${ALIGNED_DIR}/mapping_logs/temporal_bowtie_output.txt")
	#bowties terminal oputput gives a small summary of the results with each line of the summary starting with a '#', this code exploits this to grep only those lines and peg them to the variable 'results'
	
	#echo "$result"                 #UNCOMMENT IF YOU WANT TO SEE THE MAPPING RESULTS ON TERMINAL IN REAL TIME
	echo "$result" >> "${ALIGNED_DIR}/mapping_logs/mapping_stats.txt"
	echo "------------------------------" >> "${ALIGNED_DIR}/mapping_logs/mapping_stats.txt"
	#adds the summary and a separator to the .txt file for easy sample identification
	
	> ${ALIGNED_DIR}/mapping_logs/temporal_bowtie_output.txt

	echo "file: $file done"
	#feedbacks into the terminal to notify the user that the file has been processed
done


	# #clears the temporal .txt where the bowtie output is being stored so the new file in the for loop can write on it from cero
	# samtools view -bS "$outdiremapping/${base1}.sam" | samtools sort -o "$outdiremapping/${base1}.sorted.bam"
	# #takes the .bam that was just created with bowtie and creates its sorted.bam version
	
	# samtools index "$outdiremapping/${base1}.sorted.bam"
	# #creates the sorted.bam.bai version of the .bam file so it can be manipulated with something like umitools
	