// Parameters
params.sample = "${projectDir}/../Samples/merged/*.fastq*" //Directory with relevant files (preferable only .fastq)
params.bowtie_index = "${projectDir}/../Reference/hsa_mature_index"

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

process COUNTS_VECTOR {
    tag "${sample}"

    input:
    tuple val(sample),path(sorted_dedup_bam),path(indexed_dedup_bam)

    output:
    path("${sample}.dedup.idxstats.csv")
    
    script:
    """
    (echo "Reference,Length,Mapped,Unmapped"; samtools idxstats ${sorted_dedup_bam} | tr '\t' ',') > "${sample}.dedup.idxstats.csv"
    """
}
workflow {
    main:
    reads_ch = channel
        .fromPath(params.sample, checkIfExists: true)
        .map { fastq ->
            tuple(fastq.simpleName, fastq)
        }
    // reads_ch.view { sample, fastq ->
    //     "INPUT -> ${sample} | ${fastq}"
    // }

    extracted_ch = UMI_EXTRACTION(reads_ch)

    // extracted_ch.view { sample, fastq ->
    //     "OUTPUT -> ${sample} | ${fastq}"
    // }

    aligned_ch = BOWTIE_ALIGNMENT(extracted_ch)

    sorted_ch = SORT_BAM(aligned_ch)

    dedup_ch = UMI_DEDUP(sorted_ch)

    counts_ch = COUNTS_VECTOR(dedup_ch)

    publish:
    counts = counts_ch
}

output {
    counts {mode 'copy'}
}