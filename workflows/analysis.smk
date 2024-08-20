from scripts import utils

configfile: "config.yaml"
GENOME_DIR = Path(config["genome_dir"])
RESULTS = Path(config["results"])
samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])
utils.samples = samples
utils.config = config
 

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
        coverage = f"{RESULTS}/methylation/{{sample}}_pe.deduplicated.bismark.cov.gz"
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
