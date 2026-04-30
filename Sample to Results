#!/bin/bash

themaindir="$PWD"

echo "What is the identifier for lane 1 samples eg. _L001_R1_001 (do not write .fastq.gz)"
read L1ID
echo "What is the identifier for lane 2 samples eg. _L002_R1_001 (do not write .fastq.gz)"
read L2ID

outdirmerged="merged_fastqs"
mkdir -p "$outdirmerged"

echo "The following files have duplicates:"
find . -mindepth 2 -maxdepth 2 -type f -printf '%f\n' | sort | uniq -d
echo ""
echo "If there are files listed above please press ctrl+c to stop the script and decide what file to keep"

echo "If there are no files listed continue"
read -r -p "Press Enter to continue..." _


find . -mindepth 2 -maxdepth 2 -type f -exec mv -i {} . \;

for f in *"${L1ID}".fastq.gz; do
	base="${f%${L1ID}.fastq.gz}"
	echo "Merging $base$L1ID.fastq.gz and $base$L2ID.fastq.gz"
	cat "${base}${L1ID}.fastq.gz" "${base}${L2ID}.fastq.gz" > "$outdirmerged/${base}.fastq.gz"
	echo "Done with: $base"
	echo ""
done

read -r -p "Proceed with UMI extraction? Press Enter to continue..." _

cd "$outdirmerged" || exit 1

outdirextracted="UMIextracted_fastqs"
mkdir -p "$outdirextracted"

#Important notice, change the pattern in --bc-pattern to match your UMI barcode pattern eg. AACTGTAGGCACCATCAAT when [fragment of interest][AACTGTAGGCACCATCAAT][12bp UMI] is the structure of the library

for file in *.fastq.gz; do
	base1=$(basename "$file" .fastq.gz)
	echo "Processing file: $file"
	umi_tools extract --extract-method=regex --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.*)$' -I "$file" -S "$outdirextracted/${base1}.umi_extracted.fastq.gz"
	echo "Completed file: $file"
	echo ""
done

read -r -p "Proceed with mapping using bowtie2? Press Enter to continue..." _

cd "$outdirextracted" || exit 1

cp -i "$themaindir/hsa_mature_DNA.fa" .

outdiremapping="Mapped"
mkdir -p "$outdiremapping"
mkdir -p mapping_logs
> mapping_logs/mapping_stats.txt
> mapping_logs/temporal_bowtie_output.txt

bowtie-build hsa_mature_DNA.fa hsa_mature_index

for f in *.fastq.gz; do
	base1=$(basename "$f" .fastq.gz)
	#creates a variable with the name of the file without the extention to name all subsequent files created from this .fastq
	
	echo "Sample: $base1" >> mapping_logs/mapping_stats.txt
	#wites the name of the file (without extension) to a new line in the .txt file
	
	bowtie -q -S hsa_mature_index "$f" > "$outdiremapping/${base1}.sam" 2>> mapping_logs/temporal_bowtie_output.txt
	# -q for .fastq files (change it to the corresponding option accordingly, check bowtie --h), this line writes the termianls output into a temporal .txt file for later grepping and creates a .sam file from the mapping
	
	result=$(grep '^#' mapping_logs/temporal_bowtie_output.txt)
	#bowties terminal oputput gives a small summary of the results with each line of the summary starting with a '#', this code exploits this to grep only those lines and peg them to the variable 'results'
	
	#echo "$result" #UNCOMMENT IF YOU WANT TO SEE THE MAPPING RESULTS ON TERMINAL IN REAL TIME
	echo "$result" >> mapping_logs/mapping_stats.txt
	echo "------------------------------" >> mapping_logs/mapping_stats.txt
	#adds the summary and a separator to the .txt file for easy sample identification
	
	> mapping_logs/temporal_bowtie_output.txt
	#clears the temporal .txt where the bowtie output is being stored so the new file in the for loop can write on it from cero
	samtools view -bS "$outdiremapping/${base1}.sam" | samtools sort -o "$outdiremapping/${base1}.sorted.bam"
	#takes the .bam that was just created with bowtie and creates its sorted.bam version
	
	samtools index "$outdiremapping/${base1}.sorted.bam"
	#creates the sorted.bam.bai version of the .bam file so it can be manipulated with something like umitools
	
	echo "file: $base1 done"
	#feedbacks into the terminal to notify the user that the file has been processed
done

cd "$outdiremapping" || exit 1

outdirdedup="Deduped"
mkdir -p "$outdirdedup"

for f in *.sorted.bam; do
	base2=$(basename "$f" .sorted.bam)
	echo "Processing file: $f"
	
	umi_tools dedup --stdin="$f" --stdout="$outdirdedup/${base2}.sorted.dedup.bam"
	(echo "Reference,Length,Mapped,Unmapped"; samtools idxstats "$outdirdedup/${base2}.sorted.dedup.bam" | tr '\t' ',') > "$outdirdedup/${base2}.dedup.idxstats.csv"


done


