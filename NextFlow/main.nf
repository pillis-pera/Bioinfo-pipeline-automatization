// UMI Extraction Process
process UMI_EXTRACTION {
    tag "${sample}"

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.extracted.fastq")

    script:
    """
    umi_tools extract \
        --extract-method=regex --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.*)\$' \
        -I "${fastq}" \
        -S "${sample}.extracted.fastq"
    """
}

// Bowtie alignment
process BOWTIE_ALIGNMENT{
    tag "${sample}"

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.sam")

    script:
    """
    bowtie \
        -p ${task.cpus} \
        -q \
        -S \
        ${params.bowtie_index} \
        ${fastq} \
        > ${sample}.sam
    """
}

// Sorting steps
process SORT_BAM {
    tag "${sample}"

    input:
    tuple val(sample),path(sam)

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")

    script:
    """
    samtools view -bS ${sam} | samtools sort -@ ${task.cpus} -o "${sample}.sorted.bam"

    samtools index "${sample}.sorted.bam"
    """
}

// UMI deduplication
process UMI_DEDUP {
    tag "${sample}"

    input:
    tuple val(sample),path(sorted_bam),path(indexed_bam)

    output:
    tuple val(sample), path("${sample}.sorted.dedup.bam"), path("${sample}.sorted.dedup.bam.bai")

    script:
    """
    umi_tools dedup --method=${task.ext.method} --stdin=${sorted_bam} --stdout="${sample}.sorted.dedup.bam"

    samtools index "${sample}.sorted.dedup.bam"
    """
}

// Final counts
process COUNTS_VECTOR {
    tag "${sample}"

    input:
    tuple val(sample),path(sorted_dedup_bam),path(indexed_dedup_bam)

    output:
    path("${sample}.csv")
    
    script:
    """
    (echo "Reference,Length,Mapped,Unmapped"; samtools idxstats ${sorted_dedup_bam} | tr '\t' ',') > "${sample}.csv"
    """
}

// Merging samples
process MERGING {
    tag "${base}"

    input:
    tuple val(base), val(type), path(samples)

    output:
    tuple val(base), path("${base}${type}")

    script:
    """
    cat ${samples} > ${base}${type}
    """
}

workflow MIRNA_PIPELINE {
    take:
    sample
    
    main:
    extracted_ch = UMI_EXTRACTION(sample)

    aligned_ch = BOWTIE_ALIGNMENT(extracted_ch)

    sorted_ch = SORT_BAM(aligned_ch)

    dedup_ch = UMI_DEDUP(sorted_ch)

    counts_ch = COUNTS_VECTOR(dedup_ch)

    emit:
    counts = counts_ch
}

workflow MERGE_LANES {
    take:
    samples
    
    main:
    merging_ch = samples
        .map { sample ->
            def stem = sample.simpleName
            def base = sample.name.split('_L00')[0]
            def type = sample.name.substring(stem.size())
            tuple(base, type, sample)
        }
        .groupTuple(by: [0, 1])

    merged_ch = MERGING(merging_ch)

    emit:
    merged = merged_ch
}

workflow {
    main:
    if (params.merge_lanes) {

        unmerged_ch = channel.fromPath(
            "${params.unmerged}/*.fastq*",
            checkIfExists: true)

        MERGE_LANES(unmerged_ch)

        samples_ch = MERGE_LANES.out.merged

    } else {

        samples_ch = channel
            .fromPath("${params.merged}/*.fastq*", 
            checkIfExists: true)
            .map{ fastq -> tuple(fastq.simpleName, fastq)}
    }
    samples_ch.view {"Pipeline input: $samples_ch"}
    MIRNA_PIPELINE(samples_ch)

    publish:
    counts = MIRNA_PIPELINE.out.counts
}

output {
    counts { 
        mode 'copy'
    }
}