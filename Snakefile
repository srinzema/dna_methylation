from scripts import utils
from pathlib import Path
import pandas as pd


configfile: "config.yaml"
GENOME_DIR = Path(config["genome_dir"])
RESULTS = Path(config["results"])


def get_trimmed_reads(wildcards):
    sample_info = samples.loc[samples["alias"] == wildcards.sample].iloc[0]
    if pd.isna(sample_info["read2"]): # SE sample
        return [f"{RESULTS}/trimmed/{wildcards.sample}.trimmed.fastq.gz"]
    return [f"{RESULTS}/trimmed/{wildcards.sample}_R1.trimmed.fastq.gz", f"{RESULTS}/trimmed/{wildcards.sample}_R2.trimmed.fastq.gz"]


rule bismark_bowtie2:
    input:
        genome = GENOME_DIR,
        CT = expand(f"{GENOME_DIR}/Bisulfite_Genome/CT_conversion/{file}", file = ["BS_CT.1.bt2", "BS_CT.2.bt2", "BS_CT.3.bt2", "BS_CT.4.bt2", "BS_CT.rev.1.bt2", "BS_CT.rev.2.bt2", "genome_mfa.CT_conversion.fa"]),
        GA = expand(f"{GENOME_DIR}/Bisulfite_Genome/GA_conversion/{file}", file = ["BS_GA.1.bt2", "BS_GA.2.bt2", "BS_GA.3.bt2", "BS_GA.4.bt2", "BS_GA.rev.1.bt2", "BS_GA.rev.2.bt2", "genome_mfa.GA_conversion.fa"]),
        reads = get_trimmed_reads
    output:
        bam = f"{RESULTS}/alignment/{sample}.fastq_bismark.bam"
    log: f"{RESULTS}/logs/bismark_bowtie2/{sample}.log"
    threads: 8
    shell:
        "bismark -N 1 -p {threads} --bowtie2 {input.genome} -q {input.reads} > {log} 2>&1"

rule bismark_genome_preparation:
    input: GENOME_DIR
    output: 
        CT = expand(f"{GENOME_DIR}/Bisulfite_Genome/CT_conversion/{file}", file = ["BS_CT.1.bt2", "BS_CT.2.bt2", "BS_CT.3.bt2", "BS_CT.4.bt2", "BS_CT.rev.1.bt2", "BS_CT.rev.2.bt2", "genome_mfa.CT_conversion.fa"]),
        GA = expand(f"{GENOME_DIR}/Bisulfite_Genome/GA_conversion/{file}", file = ["BS_GA.1.bt2", "BS_GA.2.bt2", "BS_GA.3.bt2", "BS_GA.4.bt2", "BS_GA.rev.1.bt2", "BS_GA.rev.2.bt2", "genome_mfa.GA_conversion.fa"])
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
        fastp -w {threads} --in1 {input.read1} --in2 {input.read2} 
        --out1 {output.out1} --out2 {output.out2} 
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
        fastp -w {threads} --in1 {input}
        --out1 {output.out} 
        -h {output.html} -j {output.json} > {log} 2>&1
        """