pandas_container = "quay.io/fhcrc-microbiome/python-pandas:v1.2.1_latest"

// Process : Populate duplicated NTCs in the count table
// Context : Sometimes the input library will have duplicated guides used to generate
// synthetic NTC genes. Those duplicated guides will be lost by mageck count (only
// keeping one of the unique values). Using this process, we will add back in the counts
// for any duplicated guides
process populate_ntc {

    container "${pandas_container}"
    label "io_limited"

    input:
        tuple file("input.counts.txt"), file("treatment_sample_names.txt"), file("control_sample_names.txt"), file(library)

    output:
        tuple file("counts.txt"), file("treatment_sample_names.txt"), file("control_sample_names.txt")

    script:
    template 'populate_ntc.py'
}
