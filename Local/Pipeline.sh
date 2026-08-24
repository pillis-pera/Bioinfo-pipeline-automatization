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
	umi_tools extract --extract-method=regex --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.*)$' -I "$file" -S "${UMI_EXTRACTED_DIR}/${base}.umi_extracted.fastq.gz"
	echo "Completed file: $file"
	echo ""
done

cd $UMI_EXTRACTED_DIR
