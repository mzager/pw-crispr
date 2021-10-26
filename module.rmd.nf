rmd_container = "rocker/r-rmd:latest"

// Process : Build PDFs Using RMD Files
process rmd_pdf { 

    container "${rmd_container}"
    label "mem_medium"

    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*.pdf"

    input:
        tuple file("*.R"), file("*")
        val(prefix)

    output:
        path "*"

    script:
        """/bin/bash
            set -Eeuo pipefail

            for rmd in *.R; do
                if [[ -s "\${rmd}" ]]; then
                    R < "\${rmd}\" --no-save
                fi
            done

            ls -lahtr
        """
}