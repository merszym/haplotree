process SAMTOOLS_SORT{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "Filter and Sort: $meta.id"
    label "local"

    input:
    tuple val(meta), path(extracted_bam)

    output:
    tuple val(meta), path("sorted_${extracted_bam}"), emit: bam
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    samtools view -b -u -q 25 -o filtered_${extracted_bam} ${extracted_bam}
    samtools sort $args -o sorted_${extracted_bam}  filtered_${extracted_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}