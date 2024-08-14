from scripts import utils
from pathlib import Path
import pandas as pd


configfile: "config.yaml"
GENOME_DIR = Path(config["genome_dir"])
SPIKE_IN_1 = Path(config["spikein_1"])
SPIKE_IN_2 = Path(config["spikein_2"])
RESULTS = Path(config["results"])
samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])


rule all:
    input: 
        f"{RESULTS}/qc/multiqc/multiqc_report.html",
        f"{RESULTS}/reports/bismark_summary_report.html",
        expand(f"{RESULTS}/reports/{{sample}}.html", sample=samples.alias),
        expand(f"{RESULTS}/methylation/{{sample}}_pe.deduplicated_splitting_report.txt", sample=samples.alias),
        f"{RESULTS}/coverage/summary_report.tsv",
        f"{RESULTS}/trackhub/hub.txt",

    # input: expand(f"{RESULTS}/methylation/{{sample}}_R1.trimmed_bismark_bt2_pe.deduplicated_splitting_report.txt", sample=samples.alias)


rule generate_trackhub:
    input:
        bedgraphs = expand(f"{RESULTS}/methylation/{{sample}}_pe.deduplicated.bedGraph.gz", sample=samples.alias)
    output:
        folder = directory(f"{RESULTS}/trackhub"),
        hub_file = f"{RESULTS}/trackhub/hub.txt",
        track_file = f"{RESULTS}/trackhub/tracks/trackDb.txt",
        genome_file = f"{RESULTS}/trackhub/genomes.txt",
    log: f"{RESULTS}/logs/generate_trackhub/trackhub.log"
    threads: 1
    shell: "scripts/generate_trackhub.py --output-dir {output.folder} --assembly hg38 --email slrinzema@science.ru.nl --trackfiles {input.bedgraphs} > {log} 2>&1"


rule bismark_processing_report:
    input:
        alignment = f"{RESULTS}/alignment/{{sample}}_PE_report.txt",
        dedup = f"{RESULTS}/deduplicated/{{sample}}_pe.deduplication_report.txt"
    output:
        f"{RESULTS}/reports/{{sample}}.html"
    params: 
        dir = f"{RESULTS}/reports",
        filename = lambda w: f"{w.sample}.html"
    log: f"{RESULTS}/logs/processing_report/{{sample}}.log"
    threads: 1
    shell: "bismark2report --dir {params.dir} --output {params.filename} --alignment_report {input.alignment} --dedup_report {input.dedup} > {log} 2>&1"


rule bismark_summary_report:
    input: expand(f"{RESULTS}/alignment/{{sample}}_pe.bam", sample=samples.alias)
    output: 
        f"{RESULTS}/reports/bismark_summary_report.html",
        f"{RESULTS}/reports/bismark_summary_report.txt"
    params: f"{RESULTS}/reports/bismark_summary_report"
    log: f"{RESULTS}/logs/bismark_summary_report.log"
    threads: 1
    shell: "bismark2summary --basename {params} {input} > {log} 2>&1"


rule coverage_multiqc:
    input: expand(f"{RESULTS}/coverage/{{sample}}.coverage.summary.txt", sample=samples.alias)
    output: f"{RESULTS}/qc/coverage_mqc.tsv"
    threads: 1
    shell: "scripts/coverage_multiqc.sh {output} {input}"


rule coverage_summary:
    input: expand(f"{RESULTS}/coverage/{{sample}}.coverage.summary.txt", sample=samples.alias)
    output: f"{RESULTS}/coverage/summary_report.tsv"
    log: f"{RESULTS}/logs/coverage_table/coverage_table.log"
    threads: 1
    shell: "scripts/coverage_summary.py {output} {input} > {log} 2>&1"


rule coverage_table:
    input: f"{RESULTS}/coverage/{{sample}}.coverage.CpG_report.txt"
    output: f"{RESULTS}/coverage/{{sample}}.coverage.summary.txt"
    log: f"{RESULTS}/logs/coverage_table/{{sample}}.log"
    threads: 1
    shell: "scripts/coverage_table.sh {input} {output} > {log} 2>&1"


rule coverage2cytosine:
    input:
        genome = GENOME_DIR,
        coverage = f"{RESULTS}/methylation/{{sample}}_pe.bismark.cov.gz"
    output: f"{RESULTS}/coverage/{{sample}}.coverage.CpG_report.txt"
    params: lambda w: f"{RESULTS}/coverage/{w.sample}"
    log: f"{RESULTS}/logs/coverage2cytosine/{{sample}}.log"
    threads: 1
    shell: "coverage2cytosine --genome_folder {input.genome} -o {params} {input.coverage}"


rule bismark_methylation_extractor:
    input:
        genome = GENOME_DIR,
        bam = f"{RESULTS}/deduplicated/{{sample}}_pe.deduplicated.bam"
    output: 
        report = f"{RESULTS}/methylation/{{sample}}_pe.deduplicated_splitting_report.txt",
        cov = f"{RESULTS}/methylation/{{sample}}_pe.deduplicated.bismark.cov.gz",
        bedgraph = f"{RESULTS}/methylation/{{sample}}_pe.deduplicated.bedGraph.gz",
    params: f"{RESULTS}/methylation/"
    log: f"{RESULTS}/logs/bismark_methylation_extractor/{{sample}}.log"
    threads: 12 # This is divided by three in the command, because bismark uses three times as much cores than the specified amount. 
    shell: #"bismark_methylation_extractor {input} --no_overlap --output_dir {params} > {log} 2>&1"
        """
        bismark_methylation_extractor {input.bam} --genome_folder {input.genome} \
        --no_overlap --comprehensive --gzip --CX --cytosine_report \
        --parallel 4 --output_dir {params} > {log} 2>&1
        """


rule deduplicate_bismark:
    input: f"{RESULTS}/alignment/{{sample}}_pe.bam"
    output: 
        bam = f"{RESULTS}/deduplicated/{{sample}}_pe.deduplicated.bam",
        report = f"{RESULTS}/deduplicated/{{sample}}_pe.deduplication_report.txt"
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
        CT = expand(f"{GENOME_DIR}/Bisulfite_Genome/CT_conversion/BS_CT.{{n}}.bt2", n=[1, 2, 3, 4]),
        GA = expand(f"{GENOME_DIR}/Bisulfite_Genome/GA_conversion/BS_GA.{{n}}.bt2", n=[1, 2, 3, 4]),
        genome = GENOME_DIR,
        reads = get_trimmed_reads
    output:
        bam = f"{RESULTS}/alignment/{{sample}}_pe.bam",
        report = f"{RESULTS}/alignment/{{sample}}_PE_report.txt"
    params:
        output_dir = f"{RESULTS}/alignment/",
        sample = lambda w: w.sample
    log: f"{RESULTS}/logs/bismark_bowtie2/{{sample}}.log"
    threads: 16
    shell: "bismark -N 1 -p 8 --bowtie2 --basename {params.sample} --output_dir {params.output_dir} --genome_folder {input.genome} -1 {input.reads[0]} -2 {input.reads[1]} > {log} 2>&1"


rule bismark_genome_preparation:
    input: GENOME_DIR
    output:
        CT = expand(f"{GENOME_DIR}/Bisulfite_Genome/CT_conversion/BS_CT.{{n}}.bt2", n=[1, 2, 3, 4]),
        GA = expand(f"{GENOME_DIR}/Bisulfite_Genome/GA_conversion/BS_GA.{{n}}.bt2", n=[1, 2, 3, 4]),
    log: f"{RESULTS}/logs/bismark/genome_preparation.log"
    threads: 16
    shell: "bismark_genome_preparation --bowtie2 --parallel {threads} {input} > {log} 2>&1"


rule multiqc:
    input: 
        jsons = expand(f"{RESULTS}/qc/fastp/{{sample}}.json", sample=samples.alias),
        alignment_reports = expand(f"{RESULTS}/alignment/{{sample}}_PE_report.txt", sample=samples.alias),
        deduplication_reports = expand(f"{RESULTS}/deduplicated/{{sample}}_pe.deduplication_report.txt", sample=samples.alias),
        summary_report = f"{RESULTS}/reports/bismark_summary_report.txt",
        coverage_reports = expand(f"{RESULTS}/coverage/{{sample}}.coverage.CpG_report.txt", sample=samples.alias),
        summary_coverage = f"{RESULTS}/qc/coverage_mqc.tsv"
    output:
        f"{RESULTS}/qc/multiqc/multiqc_report.html"
    params: 
        indir = f"{RESULTS}/qc/",
        outdir = f"{RESULTS}/qc/multiqc"
    log: f"{RESULTS}/logs/multiqc.log"
    threads: 1
    shell: "printf '%s\n' {input} > {params.outdir}/files.txt; multiqc --file-list {params.outdir}/files.txt --outdir {params.outdir} -v -f > {log} 2>&1"

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