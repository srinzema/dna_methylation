from scripts import utils
from pathlib import Path
import pandas as pd


configfile: "config.yaml"
include: "workflows/preprocessing.smk"
include: "workflows/alignment.smk"
include: "workflows/analysis.smk"
include: "workflows/reports.smk"
include: "workflows/trackhub.smk"

GENOME_DIR = Path(config["genome_dir"])
RESULTS = Path(config["results"])
samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])
utils.samples = samples
utils.config = config


def check_input(wildcards):
    input = [
        f"{RESULTS}/qc/multiqc/multiqc_report.html",
        f"{RESULTS}/coverage/summary_report.tsv",
    ]

    if config["trackhub"]:
        input.append(f"{RESULTS}/trackhub/hub.txt")
    print(input)
    return input


rule all:
    input:
        expand(f"{RESULTS}/reports/{{sample}}.html", sample=samples.alias),
        check_input
