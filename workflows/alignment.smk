from scripts import utils

configfile: "config.yaml"
GENOME_DIR = Path(config["genome_dir"])
RESULTS = Path(config["results"])
samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])
utils.samples = samples
utils.config = config


rule deduplicate_bismark:
    input: f"{RESULTS}/alignment/{{sample}}_pe.bam"
    output: 
        bam = f"{RESULTS}/deduplicated/{{sample}}_pe.deduplicated.bam",
        report = f"{RESULTS}/deduplicated/{{sample}}_pe.deduplication_report.txt"
    params: f"{RESULTS}/deduplicated/"
    log: f"{RESULTS}/logs/deduplicate_bismark/{{sample}}.log"
    threads: 1
    shell: "deduplicate_bismark -p --output_dir {params} {input} > {log} 2>&1"


rule bismark_bowtie2:
    input:
        CT = expand(f"{GENOME_DIR}/Bisulfite_Genome/CT_conversion/BS_CT.{{n}}.bt2", n=[1, 2, 3, 4]),
        GA = expand(f"{GENOME_DIR}/Bisulfite_Genome/GA_conversion/BS_GA.{{n}}.bt2", n=[1, 2, 3, 4]),
        genome = GENOME_DIR,
        reads = utils.get_trimmed_reads
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
