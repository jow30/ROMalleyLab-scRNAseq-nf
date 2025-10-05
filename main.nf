#!/usr/bin/env nextflow

nextflow.enable.dsl=2 // enables the modern DSL2 syntax

// Define a workflow
workflow {
    SAMPLESHEET_VERIFICATION() | view()
}

process SAMPLESHEET_VERIFICATION {
    output:
    stdout

    script:
    """
    echo "Samplesheet looks good!"
    """
}
