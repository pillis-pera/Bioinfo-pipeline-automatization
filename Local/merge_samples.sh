#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
UNMERGED_DIR="${REPO_DIR}/Samples/unmerged"
MERGED_DIR="${REPO_DIR}/Samples/merged"

cd ${UNMERGED_DIR}          #File structure     [sampleID]_[laneID].[Filetype]
L1ID="_L001"                #Set the laneID     (eg. L001)
filetype=".fastq.gz"        #Set                (eg. .fastq.gz)
                            #Important to include .  ^ 
for f in *"${L1ID}"*; do
    base="${f%${L1ID}*}"
    echo "Merging $base"
    cat $base* > "$MERGED_DIR/$base.merged$filetype"
done