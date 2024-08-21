# DNA Methylation Sequencing Analysis Pipeline

This pipeline automates the analysis of DNA methylation sequencing data, from preprocessing raw data to generating comprehensive reports and visualization files. It is ideal for researchers conducting epigenetic studies who require a streamlined and reproducible workflow.

## Usage

Clone the repository: Start by cloning the pipeline repository to your local machine.
```bash
cd <installation directory> # Git clone makes a subdirectory with the name of the repository
git clone git@github.com:srinzema/dna_methylation.git
cd dna_methylation
```

## Configure the pipeline:
Edit the `config.yaml` file to specify the necessary paths and settings for your analysis, and the `samples.tsv` to contain the names of your samples excluding file and sequencing extensions.

### Config.yaml

```
fastq_dir: Directory containing the raw FASTQ files.
genome_dir: Path to the reference genome you are using for alignment.
results: Directory where all processed data, reports, and visualizations will be saved.
samplesheet: A tab-separated file that lists your samples and any relevant metadata.
```

Example config.yaml:

```
fastq_dir: "/vol/sequencing/samples"
genome_dir: "/home/srinzema/genomes/hg38"
results: "results_test" # Relative to the pipeline
samplesheet: "samples.tsv" # Relative to the pipeline
```

### Samples.tsv

```
samples
sample_wildtype
sequenced_sample_1
sequenced_sample_7
```

Each row in the samples.tsv file represents a sample to be processed by the pipeline. The samples column contains the identifiers for each sample, which should correspond to the filenames of your FASTQ files in the `fastq_dir` specified in the `config.yaml`.

## Run the pipeline

Execute the pipeline: Run the entire pipeline using Snakemake. This will execute all necessary subworkflows in the correct order.

```bash
snakemake --cores <number_of_cores>
```

Check results: After the pipeline completes, find all processed data, reports, and visualization files in the output directories specified in the configuration.

This pipeline is flexible, allowing individual subworkflows to be run separately if needed, but the global usage will handle the entire analysis from start to finish. For more information on the subworkflows and their rules see the [workflow readme.md](workflows/README.md)