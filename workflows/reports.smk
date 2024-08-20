from scripts import utils

configfile: "config.yaml"
GENOME_DIR = Path(config["genome_dir"])
RESULTS = Path(config["results"])
samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])
utils.samples = samples
utils.config = config


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


rule coverage_multiqc:
    input: expand(f"{RESULTS}/coverage/{{sample}}.coverage.summary.txt", sample=samples.alias)
    output: f"{RESULTS}/qc/coverage_mqc.tsv"
    threads: 1
    shell: "scripts/coverage_multiqc.sh {output} {input}"


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
