#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.scale_cutoff = 1
params.trim_3_prime = -8
params.trim_5_prime = 32

// Import the modules
include {
    cutadapt_trim as Process_Cutadapt_Trim_Treatment;
    cutadapt_trim as Process_Cutadapt_Trim_Control;
} from './module.cutadapt'

include {
    mageck_count as Process_Mageck_Count_Treatment;
    mageck_count as Process_Mageck_Count_Control;
    mageck_rra as Process_Mageck_Rra;
    mageck_mle as Process_Mageck_Mle;
    mageck_pathway as Process_Mageck_Pathway;
    mageck_merge as Process_Mageck_Merge;
} from './module.mageck'

include {
    mageckflute_rra as Process_MageckFlute_Rra;
    mageckflute_mle as Process_MageckFlute_Mle
} from './module.mageckflute'

include {
    rmd_pdf as Process_Rmd_Pdf_Treatment;
    rmd_pdf as Process_Rmd_Pdf_Control;
} from './module.rmd'

include {
    populate_ntc as Process_Populate_NTC
} from './module.general'

// Validate Input Parameters / Print Help
def validate(params) {

    if (params.treatment_fastq 
        && params.control_fastq 
        && params.library
        && params.trim_5_prime
        && params.trim_3_prime 
        && params.organism 
        && params.scale_cutoff 
        && params.gmt 
        && params.output)
        { return; }
    
    log.info"""
        Usage:
        nextflow run FredHutch/crispr-screen-nf

        Required Arguments:
            --treatment_fastq   Path to FASTQ data for treatment samples
            --control_fastq     Path to FASTQ data for control samples
                                Multiple files can be specified with wildcards and commas, e.g.
                                    /path/to/inputs/A/*.fastq.gz,/path/to/inputs/B/*.fq.gz
            --library           Text file describing sgRNA library
                                    As described at https://sourceforge.net/p/mageck/wiki/input/
            --trim_5_prime      Amount to trim from 5 prime end (default: 32)
            --trim_3_prime      Amount to trim from 3 prime end (default: -8)
            --organism          Organism string provided for MAGeCK-Flute (default: hsa)
            --scale_cutoff      Parameter 'scale_cutoff' for MAGeCK-Flute (default: 1)
            --gmt               Pathway GMT File
            --deepmap_effect    "https://ndownloader.figshare.com/files/20234073" - Effect File Can't Download From Different Region
            --deepmap_samples   "https://ndownloader.figshare.com/files/20274744" - Sample File Can't Download From Different Region
            --output            Path to output directory
            
        Open Questions:
            - ntc_list          Path to file describing negative controls
                                    As described in https://sourceforge.net/p/mageck/wiki/input/#negative-control-sgrna-list

        Design Decisions
            - treatname & controlname should be standardize for first round of portal development
            - mle design matrix will be automatically generated and include treatment v control + paired sample where appropriate
            
        Todo:
            - Ingest MSigDB info Into HutchCloud
            - Dropdown of GMT options
            - Autogenerate the mle design matrix https://sourceforge.net/p/mageck/wiki/demo/#the-fourth-tutorial-using-mageck-mle-module
            - If species is not (hsa) or (mle) skip flute

        """
    
    exit 1
}

workflow {

    main:

    // Validate Input + Print Help On Fail 
    validate(params)
    
    // Channel : Treatment FastQ
    Channel.fromPath(params.treatment_fastq.split(',').toList()).set{Channel_Fastq_Treatment}

    // Chanel : Control FastQ 
    Channel.fromPath(params.control_fastq.split(',').toList()).set{Channel_Fastq_Control}

    // Chanel : SGRNA Library
    Channel.fromPath(params.library).set{Channel_Library}

    // Process : Cutadapt Trim
    // Remove Extra Adaptor Sequences From Reads
    Process_Cutadapt_Trim_Treatment(Channel_Fastq_Treatment, params.trim_5_prime, params.trim_3_prime, 'cutadapt/trim/treatment')
    Process_Cutadapt_Trim_Control(Channel_Fastq_Control, params.trim_5_prime, params.trim_3_prime, 'cutadapt/trim/control')

    // Process : Mageck Count
    // Map the raw FASTQ data to reference library file and count the reads for each sgRNA
    Process_Mageck_Count_Treatment(Process_Cutadapt_Trim_Treatment.out.combine(Channel_Library), 'mageck/count/treatment')
    Process_Mageck_Count_Control(Process_Cutadapt_Trim_Control.out.combine(Channel_Library), 'mageck/count/control')

    // Process : Mageck Merge
    // Concat All Count Data
    Process_Mageck_Merge(Process_Mageck_Count_Treatment.out.counts.toSortedList(), Process_Mageck_Count_Control.out.counts.toSortedList(), 'mageck/count/combined')

    // Add back counts for duplicated guides (if any)
    Process_Populate_NTC(Process_Mageck_Merge.out.merged.combine(Channel_Library))
    
    // Process : Mageck Rra
    // MAGeCK RRA (identifying CRISPR screen hits by calculating the RRA enrichment score to indicate the essentiality of a gene)
    Process_Mageck_Rra(Process_Populate_NTC.out, 'rra', 'mageck/rra')

    // Process : Mageck MLE
    // MAGeCK MLE (identifying CRISPR screen hits by calculating a ‘beta score’ for each targeted gene to measure the degree of selection after the target is perturbed)
    Process_Mageck_Mle(Process_Populate_NTC.out, params.treatment_fastq, params.control_fastq, 'mle', 'mageck/mle')

    // Process : Mageck Pathway
    Process_Mageck_Pathway(Process_Mageck_Rra.out.geneSummary.combine(Channel.fromPath(params.gmt)), 'gene', 'mageck/rra/pathway')
    
    // Run Mageck Flute RRA
    Process_MageckFlute_Rra(Process_Mageck_Rra.out.geneSummary, Process_Mageck_Rra.out.sgrnaSummary, params.scale_cutoff, params.organism, 'mageckflute/rra')

    // Process Mageck Flute MLE
    Process_MageckFlute_Mle(Process_Mageck_Mle.out.geneSummary, Channel.fromPath(params.depmap_effect), Channel.fromPath(params.depmap_samples), 'mageckflute/mle')
    
    // Process : Rmd To Pdf
    // Process_Rmd_Pdf_Treatment(Process_Mageck_Count_Treatment.out.r, 'mageck/count/treatment/pdf')
    // Process_Rmd_Pdf_Control(Process_Mageck_Count_Control.out.r, 'mageck/count/control/pdf')
   
}
