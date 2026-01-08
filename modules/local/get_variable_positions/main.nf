process GET_VARIABLE_POSITIONS{
    container (workflow.containerEngine ? "pypy:3" : null)
    label 'local'
    tag "${meta.id}"

    input:
    tuple val(meta), path(pileup)

    output:
    tuple val(meta), path("positions.tsv"), emit: tsv
    path "versions.yml"                   , emit: versions

    script:
    """
    get_positions.py ${pileup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypy: \$(pypy3 --version | tail -1 | cut -d ' ' -f2)
    END_VERSIONS
    """
}