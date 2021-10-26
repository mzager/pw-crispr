mageckflute_container = "quay.io/biocontainers/bioconductor-mageckflute:1.12.0--r41hdfd78af_0"

// Process MAGeCKFlute RRA
process mageckflute_rra {

    container "${mageckflute_container}"
    label "mem_medium"
    
    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*"
    
    input:
        file(gene_summary)
        file(sgrna_summary)
        val(scale_cutoff)
        val(organism)
        val(prefix)
        
    output:
        path "*"
        
    script:
        """/bin/bash
        set -Eeuo pipefail

        mageckflute_rra.R \
            "${gene_summary}" \
            "${sgrna_summary}" \
            "${organism}" \
            "${scale_cutoff}"

        ls -lahtr
        """
}

// Process : MAGeCKFlute Mle
process mageckflute_mle {

    container "${mageckflute_container}"
    label "mem_medium"

    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*"
    

    input:
        file(gene_summary)
        file(depmap_effect)
        file(depmap_samples)
        val(prefix)

    output:
        path "*"

    script:
        """/bin/bash
        set -Eeuo pipefail

        mageckflute_mle.R \
            "${gene_summary}" \
            "${params.organism}" \
            "treatment_control" \
            "depmap" \
            "${depmap_effect}" \
            "${depmap_samples}"
            
        ls -lahtr
        """
}