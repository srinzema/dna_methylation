from scripts import utils
from pathlib import Path
import pandas as pd


configfile: "config.yaml"
GENOME_DIR = Path(config["genome_dir"])
RESULTS = Path(config["results"])
samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])

# TODO; get this from config, maybe extra file?
filter_genes = ["chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr14_GL000009v2_random", "chr14_GL000194v1_random", "chr14_KI270722v1_random", "chr14_KI270726v1_random", "chr15", "chr15_KI270727v1_random", "chr16", "chr16_KI270728v1_random", "chr17", "chr17_GL000205v2_random", "chr18", "chr19", "chr1_KI270706v1_random", "chr1_KI270708v1_random", "chr1_KI270711v1_random", "chr1_KI270712v1_random", "chr1_KI270713v1_random", "chr2", "chr20", "chr21", "chr22", "chr22_KI270731v1_random", "chr22_KI270733v1_random", "chr3", "chr3_GL000221v1_random", "chr4", "chr4_GL000008v2_random", "chr5", "chr6", "chr7", "chr8", "chr9", "chr9_KI270718v1_random", "chr9_KI270719v1_random", "chr9_KI270720v1_random", "chrM", "chrUn_GL000195v1", "chrUn_GL000213v1", "chrUn_GL000214v1", "chrUn_GL000218v1", "chrUn_GL000219v1", "chrUn_GL000220v1", "chrUn_GL000224v1", "chrUn_KI270442v1", "chrUn_KI270741v1", "chrUn_KI270742v1", "chrUn_KI270743v1", "chrUn_KI270744v1", "chrUn_KI270745v1", "chrUn_KI270746v1", "chrUn_KI270748v1", "chrUn_KI270750v1", "chrUn_KI270751v1", "chrUn_KI270754v1", "chrUn_KI270755v1", "chrX", "chrY"]


rule all:
    input: 
        f"{RESULTS}/qc/multiqc/multiqc_report.html",
        f"{RESULTS}/reports/bismark_summary_report.html",
        expand(f"{RESULTS}/reports/{{sample}}.html", sample=samples.alias),
        expand(f"{RESULTS}/methylation/{{sample}}_pe.deduplicated_splitting_report.txt", sample=samples.alias),
        expand(f"{RESULTS}/methylation/{{sample}}.filtered.bedGraph.gz", sample=samples.alias),
        f"{RESULTS}/coverage/summary_report.tsv",
        expand(f"{RESULTS}/trackhub/hg38/{{sample}}.filtered.bigwig", sample=samples.alias),
        f"{RESULTS}/trackhub/hub.txt"


rule generate_trackhub:
    input: expand(f"{RESULTS}/trackhub/hg38/{{sample}}.filtered.bigwig", sample=samples.alias),
    output:
        hub_file = f"{RESULTS}/trackhub/hub.txt",
        track_file = f"{RESULTS}/trackhub/hg38/trackDb.txt",
        genome_file = f"{RESULTS}/trackhub/genomes.txt",
    params: trackhub = f"{RESULTS}/trackhub"
    log: f"{RESULTS}/logs/generate_trackhub/trackhub.log"
    threads: 1
    shell: "scripts/generate_trackhub.py --output-dir {params.trackhub} --assembly hg38 --email slrinzema@science.ru.nl --trackfiles {input} > {log} 2>&1"


rule bedGraph_to_bigwig:
    input:
        bedgraph = f"{RESULTS}/{{sample}}_sorted.bg",
        chrom_sizes = f"{GENOME_DIR}/hg38.fa.sizes"
    output: bigwig = f"{RESULTS}/trackhub/hg38/{{sample}}.filtered.bigwig",
    params: trackhub = f"{RESULTS}/trackhub/hg38/"
    log: f"{RESULTS}/logs/bedgraph_to_bigwig/{{sample}}.log"
    shell: "bedGraphToBigWig {input.bedgraph} {input.chrom_sizes} {output.bigwig} > {log} 2>&1"


rule filter_and_sort_bedgraph:
    input: f"{RESULTS}/methylation/{{sample}}_pe.deduplicated.bedGraph.gz" #TODO; make .filtered a wildcard
    output: f"{RESULTS}/{{sample}}_sorted.bg"
    params: filter = "\|".join(filter_genes)
    log: f"{RESULTS}/logs/sort_bedgraph/{{sample}}.log"
    shell: "zcat {input} | grep '{params.filter}' | sort -k1,1 -k2,2n > {output}"


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
