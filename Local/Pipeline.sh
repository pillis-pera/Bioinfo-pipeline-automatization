#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
UNMERGED_DIR="${REPO_DIR}/Samples/unmerged"
MERGED_DIR="${REPO_DIR}/Samples/merged"
UMI_EXTRACTED_DIR="${REPO_DIR}/Local/r_umi_extracted"
FASTP_DIR="${REPO_DIR}/Local/r_fastp"
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

cd $MERGED_DIR
for file in *.fastq*; do
	base="${file%.*}"
	echo "Processing file: $file"
    # echo $base
	umi_tools extract \
		--extract-method=regex \
		--bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.*)$' \
		-I "$file" \
		--log=/dev/null |
		fastp \
		--stdin  \
		-l 17 \
		-o "${UMI_EXTRACTED_DIR}/${base}.umi_extracted.fastq.gz" \
		--html "${FASTP_DIR}/${base}.fastp.html" \
		--json "${FASTP_DIR}/${base}.fastp.json"
		
	echo "Completed file: $file"
	echo ""
done

cd $UMI_EXTRACTED_DIR
bowtie-build "${REFERENCE_DIR}/hsa_mature_DNA.fa" "${REFERENCE_DIR}/hsa_mature_index"
> ${ALIGNED_DIR}/mapping_logs/mapping_stats.txt
for file in *.fastq*; do
	base="${file%.umi_extracted.fastq*}"
	#creates a variable with the name of the file without the extention to name all subsequent files created from this .fastq
	
	echo "Sample: $base" >> "${ALIGNED_DIR}/mapping_logs/mapping_stats.txt"
	#writes the name of the file (without extension) to a new line in the .txt file
	
	bowtie -q -S "${REFERENCE_DIR}/hsa_mature_index" "$file" \
		2>> "${ALIGNED_DIR}/mapping_logs/temporal_bowtie_output.txt" |
		samtools sort -o "${ALIGNED_DIR}/${base}.sorted.bam" - 
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

cd $ALIGNED_DIR
for file in *.bam; do
	base="${file%.sorted.bam}"
	
	samtools index "${ALIGNED_DIR}/${base}.sorted.bam"
	#creates the sorted.bam.bai version of the .bam file so it can be manipulated with something like umitools
    # echo $base
    # echo $file
done

for file in *.sorted.bam; do
	base="${file%.sorted.bam}"
	echo "Processing file: $file"
    echo $base
	umi_tools dedup --stdin="$file" --stdout="${DEDUPED_DIR}/${base}.sorted.dedup.bam"
    samtools index "${DEDUPED_DIR}/${base}.sorted.dedup.bam"
done

cd $DEDUPED_DIR
for file in *.bam; do
	base="${file%.sorted.dedup.bam}"
	(echo "Reference,Length,Mapped,Unmapped"; samtools idxstats "$file" | tr '\t' ',') > "${COUNTS_DIR}/${base}.csv"
done
