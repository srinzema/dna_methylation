from scripts import utils
from pathlib import Path
import pandas as pd


configfile: "config.yaml"
GENOME_DIR = Path(config["genome_dir"])
RESULTS = Path(config["results"])

samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])

rule all:
    input: expand(f"{RESULTS}/methylation/{{sample}}_R1.trimmed_bismark_bt2_pe.deduplicated_splitting_report.txt", sample=samples.alias)


rule bismark_methylation_extractor:
    input: f"{RESULTS}/deduplicated/{{sample}}_R1.trimmed_bismark_bt2_pe.deduplicated.bam"
    output: f"{RESULTS}/methylation/{{sample}}_R1.trimmed_bismark_bt2_pe.deduplicated_splitting_report.txt"
    params: f"{RESULTS}/methylation/"
    log: f"{RESULTS}/logs/bismark_methylation_extractor/{{sample}}.log"
    threads: 1
    shell: "bismark_methylation_extractor {input} --no_overlap --output_dir {params} > {log} 2>&1"


rule deduplicate_bismark:
    input: f"{RESULTS}/alignment/{{sample}}_R1.trimmed_bismark_bt2_pe.bam"
    output: f"{RESULTS}/deduplicated/{{sample}}_R1.trimmed_bismark_bt2_pe.deduplicated.bam"
    params: f"{RESULTS}/deduplicated/"
    log: f"{RESULTS}/logs/deduplicate_bismark/{{sample}}.log"
    threads: 1
    shell: "deduplicate_bismark -p --output_dir {params} {input} > {log} 2>&1"


def get_trimmed_reads(wildcards):
    sample_info = samples.loc[samples["alias"] == wildcards.sample].iloc[0]
    if pd.isna(sample_info["read2"]): # SE sample
        return [f"{RESULTS}/trimmed/{wildcards.sample}.trimmed.fastq.gz"]
    return [f"{RESULTS}/trimmed/{wildcards.sample}_R1.trimmed.fastq.gz", f"{RESULTS}/trimmed/{wildcards.sample}_R2.trimmed.fastq.gz"]


rule bismark_bowtie2:
    input:
        genome = GENOME_DIR,
        CT = expand(f"{GENOME_DIR}/Bisulfite_Genome/CT_conversion/BS_CT.{{ext}}.bt2", ext=[1, 2, 3, 4]),
        GA = expand(f"{GENOME_DIR}/Bisulfite_Genome/GA_conversion/BS_GA.{{ext}}.bt2", ext=[1, 2, 3, 4]),
        reads = get_trimmed_reads
    output:
        bam = f"{RESULTS}/alignment/{{sample}}_R1.trimmed_bismark_bt2_pe.bam"
    params:
        working_dir = f"{RESULTS}/trimmed",
        out = f"{RESULTS}/alignment/"
    log: f"{RESULTS}/logs/bismark_bowtie2/{{sample}}.log"
    threads: 16
    shell: "bismark -N 1 -p 8 --output_dir {params.out} --bowtie2 {input.genome} -1 {input.reads[0]} -2 {input.reads[1]} > {log} 2>&1"

rule bismark_genome_preparation:
    input: GENOME_DIR
    output: 
        CT = expand(f"{GENOME_DIR}/Bisulfite_Genome/CT_conversion/BS_CT.{{ext}}.bt2", ext=[1, 2, 3, 4]),
        GA = expand(f"{GENOME_DIR}/Bisulfite_Genome/GA_conversion/BS_GA.{{ext}}.bt2", ext=[1, 2, 3, 4]),
    log: f"{RESULTS}/logs/bismark/genome_preparation.log"
    threads: 16
    shell: "bismark_genome_preparation --bowtie2 --parallel {threads} {input} > {log} 2>&1"

# FASTP RULES

def original_read_1(wildcards):
    alias = wildcards.sample
    sample_info = samples.loc[samples["alias"] == alias]
    return sample_info.read1.iloc[0]


def original_read_2(wildcards):
    alias = wildcards.sample
    sample_info = samples.loc[samples["alias"] == alias]
    return sample_info.read2.iloc[0]


rule fastp_pe:
    input:
        read1 = original_read_1,
        read2 = original_read_2,
    output:
        json=f"{RESULTS}/qc/fastp/{{sample}}.json",
        html=f"{RESULTS}/qc/fastp/{{sample}}.html",
        out1=f"{RESULTS}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        out2=f"{RESULTS}/trimmed/{{sample}}_R2.trimmed.fastq.gz",
    log: f"{RESULTS}/logs/fastp/{{sample}}.log"
    threads: 2
    shell:
        """
        fastp -w {threads} --in1 {input.read1} --in2 {input.read2} \
        --out1 {output.out1} --out2 {output.out2} \
        -h {output.html} -j {output.json} > {log} 2>&1
        """


rule fastp_se:
    input: original_read_1
    output:
        json = f"{RESULTS}/qc/fastp/{{sample}}.json",
        html = f"{RESULTS}/qc/fastp/{{sample}}.html",
        out = f"{RESULTS}/trimmed/{{sample}}.trimmed.fastq.gz"
    log: f"{RESULTS}/logs/fastp/{{sample}}.log"
    wildcard_constraints: sample = "(?!.*_R\d)"
    threads: 2
    shell:
        """
        fastp -w {threads} --in1 {input} \
        --out1 {output.out} \
        -h {output.html} -j {output.json} > {log} 2>&1
        """