mageck_container = "quay.io/biocontainers/mageck:0.5.9.4--py38h8c62d01_1"
suffix_list = "trimmed gz fq fastq fna fasta" // List of extensions to be remove to determine sample name

// Process : MAGeCK Count
process mageck_count {

    container "${mageck_container}"
    label "mem_medium"

    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*.txt"
    publishDir "${params.output}/${prefix}/log/", mode: "copy", overwrite: "true", pattern: "*.log"

    input:
        tuple file(fastq), file(library)
        val(prefix)

    output:
        tuple file("*.R"), file("*"), emit: r
        path '*.count.txt', emit: counts
        
    script:
        """/bin/bash
        set -Eeuo pipefail

        # Parse the name of the sample from the name of the FASTQ file
        SAMPLE_NAME="${fastq.name}"
        for suffix in ${suffix_list}; do
            SAMPLE_NAME=\$(echo \$SAMPLE_NAME | sed "s/.\$suffix\$//")
        done

        echo FASTQ file is ${fastq.name}
        echo Sample name is \$SAMPLE_NAME

        mageck count \
            -l ${library} \
            -n \$SAMPLE_NAME \
            --sample-label \$SAMPLE_NAME \
            --fastq ${fastq} \
            --trim-5 0 \
            --pdf-report

        ls -lahtr
        """
}

// Process : Mageck Merge
process mageck_merge {

    container "quay.io/fhcrc-microbiome/python-pandas:v1.2.1_latest"
    label "io_limited"

    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*.txt"

    input:
        file "treatment/treatment_*.txt"
        file "control/control_*.txt"
        val(prefix)

    output:
        tuple file("counts.txt"), file("treatment_sample_names.txt"), file("control_sample_names.txt"), emit: merged

    script:
        """/bin/bash
        set -Eeuo pipefail

        join_counts.py ""

        ls -lahtr
        """
}

// Process : MAGeCK RRA
process mageck_rra {

    container "${mageck_container}"
    label "mem_medium"
    
    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*.txt"
    publishDir "${params.output}/${prefix}/log/", mode: "copy", overwrite: "true", pattern: "*.log"

    input:
        tuple file(counts_tsv), file(treatment_samples), file(control_samples)
        val(output_prefix)
        val(prefix)

    output:
        tuple file("*.R"), file("*"), emit: r
        path "${output_prefix}.gene_summary.txt", emit: geneSummary
        path "${output_prefix}.sgrna_summary.txt", emit: sgrnaSummary
        
    script:
        """/bin/bash
        set -Eeuo pipefail

        echo COUNTS file is ${counts_tsv.name}
        echo TREATMENT file is \$(cat ${treatment_samples})
        echo CONTROLS file is \$(cat ${control_samples})
        echo OUTPUT prefix is ${output_prefix}

        mageck test \
            -k ${counts_tsv} \
            -t "\$(cat ${treatment_samples})" \
            -c "\$(cat ${control_samples})" \
            -n "${output_prefix}" \
            --normcounts-to-file \
            --keep-tmp \
            --pdf-report

        ls -lahtr
        """
}

// Process : MAGeCK MLE
process mageck_mle {

    container "${mageck_container}"
    label "mem_medium" 

    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*.txt"
    publishDir "${params.output}/${prefix}/log", mode: "copy", overwrite: "true", pattern: "*.log"

    input:
        tuple file(counts_tsv), file(sample_names), file(control_names)
        val(treatment)
        val(control)
        val(output_prefix)
        val(prefix)

    output:
        path "${output_prefix}.gene_summary.txt", emit: geneSummary
        path "${output_prefix}.sgrna_summary.txt", emit: sgrnaSummary

    script:
        """/bin/bash
        set -Eeuo pipefail

        echo 'Samples	baseline	treatment_control' > design.mtx
        for i in \$(echo ${treatment} | tr ',' '\n')
        do
            echo \$(basename \$i | cut -d. -f1)'  1    1' >> design.mtx
        done
        for i in \$(echo ${control} | tr ',' '\n')
        do
            echo \$(basename \$i | cut -d. -f1)'  1    0' >> design.mtx
        done
        
        mageck mle \
            -k ${counts_tsv} \
            -d 'design.mtx' \
            -n ${output_prefix}
            
        ls -lahtr
        """
}

// Process : MAGeCK Pathway
process mageck_pathway { 
    
    container "${mageck_container}"
    label "mem_medium"
    
    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*"
    
    output:
        path "gene.pathway_summary.txt", emit: gene_summary

    input:
        tuple file(counts_tsv), file(gmt)
        val(output_prefix)
        val(prefix)

    script:
        """/bin/bash
        set -Eeuo pipefail

        mageck pathway \
            --gene-ranking ${counts_tsv} \
            --gmt-file ${gmt} \
            -n "${output_prefix}"

        ls -lahtr
        """
}