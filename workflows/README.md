# README for Snakemake Subworkflows

This pipeline contains several Snakemake subworkflows designed to process DNA methylation sequencing data. Each subworkflow is focused on a specific stage of the analysis pipeline, from preprocessing and alignment to analysis and report generation. Below is an overview of the main subworkflows and the key rules they include.
Subworkflows Overview
## preprocessing.smk

This subworkflow handles the initial preprocessing of raw sequencing data. It includes quality control and trimming of reads.

    fastp_pe: Processes paired-end sequencing data with Fastp for quality control, trimming, and generation of QC reports.
    fastp_se: Similar to fastp_pe, but for single-end sequencing data.

## alignment.smk

This subworkflow is responsible for aligning the preprocessed reads to a reference genome using Bismark, followed by deduplication.

    bismark_genome_preparation: Prepares the reference genome for bisulfite conversion.
    bismark_bowtie2: Aligns the trimmed reads to the prepared genome using Bowtie2 within Bismark.
    bismark_deduplicate: Removes duplicate reads from the aligned BAM files to ensure accurate methylation analysis.

## analysis.smk

This subworkflow focuses on analyzing the aligned data to extract methylation information and summarize coverage.

    bismark_methylation_extractor: Extracts methylation calls from the deduplicated BAM files, generating various reports.
    coverage2cytosine: Generates a coverage report for each cytosine in the genome, showing how many times each was covered by sequencing reads.
    coverage_table: Summarizes the coverage reports per chromosome.
    coverage_summary: Compiles a summary report of all coverage tables.

## reports.smk

This subworkflow generates various summary reports and a comprehensive MultiQC report that combines results from all samples.

    bismark2report: Creates an HTML report for each sample, summarizing alignment and deduplication metrics.
    bismark2summary: Generates an overall summary of all alignments.
    coverage_multiqc: Prepares coverage data for integration into a MultiQC report.
    multiqc: Generates a comprehensive MultiQC report combining results from Fastp, Bismark, and coverage analyses.

## trackhub.smk

This subworkflow prepares data for visualization on the UCSC Genome Browser by generating track hubs.

    filter_and_sort_bedgraph: Filters out unwanted data and sorts bedGraph files for conversion.
    bedGraph_to_bigwig: Converts bedGraph files to BigWig format, suitable for UCSC visualization.
    generate_trackhub: Creates a track hub containing all the necessary files for UCSC Genome Browser visualization.