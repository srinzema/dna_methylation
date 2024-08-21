from scripts import utils

configfile: "config.yaml"
GENOME_DIR = Path(config["genome_dir"])
RESULTS = Path(config["results"])
samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])
utils.samples = samples
utils.config = config
 

rule coverage_summary:
    "This rule makes a summary of all reports generated by the rule coverage_table."

    input: expand(f"{RESULTS}/coverage/{{sample}}.coverage.summary.txt", sample=samples.alias)
    output: f"{RESULTS}/coverage/summary_report.tsv"
    log: f"{RESULTS}/logs/coverage_table/coverage_table.log"
    threads: 1
    shell: "scripts/coverage_summary.py {output} {input} > {log} 2>&1"


rule coverage_table:
    "Here I summarize the output from coverage2cytosine per chromosome."

    input: f"{RESULTS}/coverage/{{sample}}.coverage.CpG_report.txt"
    output: f"{RESULTS}/coverage/{{sample}}.coverage.summary.txt"
    log: f"{RESULTS}/logs/coverage_table/{{sample}}.log"
    threads: 1
    shell: "scripts/coverage_table.sh {input} {output} > {log} 2>&1"


rule coverage2cytosine:
    """
    Counts how many times each cytosine was covered by your DNA sequencing reads.
    Results in a complete map, including areas that might not be well-covered by sequencing.
    """

    input:
        genome = GENOME_DIR,
        coverage = f"{RESULTS}/methylation/{{sample}}_pe.deduplicated.bismark.cov.gz"
    output: f"{RESULTS}/coverage/{{sample}}.coverage.CpG_report.txt"
    params: lambda w: f"{RESULTS}/coverage/{w.sample}"
    log: f"{RESULTS}/logs/coverage2cytosine/{{sample}}.log"
    threads: 1
    shell: "coverage2cytosine --genome_folder {input.genome} -o {params} {input.coverage} > {log} 2>&1"


rule bismark_methylation_extractor:
    """
    This rule processes the output files generated by bismark_bowtie2 and bismark_deduplicate.
    Extracts the methylation calls for individual cytosines and outputs by context (Cpg, CHG, CHH) 
    and methylation state (methylated '+' or non-methylated '-').

    Parameters:
    --no_overlap: In case of overlap between R1 and R2, only R1 is used for cytosine methylation calls.
    --comprehensive: Pools methylation information from all strands into context-dependent files.
    --cytosine_report: Generates a comprehensive report on the methylation status of all cytosines in the genome.
    --CX: Allows analysis of all cytosine contexts beyond just CpG.
    """
    
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
    shell:
        """
        bismark_methylation_extractor {input.bam} \
        --genome_folder {input.genome} --output_dir {params} \
        --no_overlap --comprehensive --gzip --CX --cytosine_report \
        --parallel 4 > {log} 2>&1
        """
