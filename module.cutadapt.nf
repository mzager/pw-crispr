cutadapt_container = "quay.io/biocontainers/cutadapt:3.4--py37h73a75cf_1"
suffix_list = "trimmed gz fq fastq fna fasta" // List of extensions to be remove to determin sample name

// Process : Cutadapt Trim
process cutadapt_trim {
    
    container "${cutadapt_container}"
    label "mem_medium"


    publishDir "${params.output}/${prefix}/fastq/", mode: "copy", overwrite: "true"

    input:
        file(fastq)
        val(trim_5_prime)
        val(trim_3_prime)
        val(prefix)
        
    output: 
        file "${fastq.name.replaceAll(".gz","")}"

    script:
        """/bin/bash
        set -Eeuo pipefail
        
        echo FASTQ file is ${fastq.name}
        echo Prime5 is ${trim_5_prime}
        echo Prime3 is ${trim_3_prime}

        # Parse the name of the sample from the name of the FASTQ file
        SAMPLE_NAME="${fastq.name}"
        for suffix in ${suffix_list}; do
            SAMPLE_NAME=\$(echo \$SAMPLE_NAME | sed "s/.\$suffix\$//")
        done

        cutadapt \
            -u "${trim_5_prime}" \
            -u "${trim_3_prime}" \
            -o \$SAMPLE_NAME".fastq" "${fastq.name}" \
            
        ls -lahtr
        """
}