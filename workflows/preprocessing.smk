from scripts import utils

configfile: "config.yaml"
RESULTS = Path(config["results"])
samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])
utils.samples = samples
utils.config = config


rule fastp_pe:
    "Runs fastp on a paired-end sample."

    input:
        read1 = utils.original_read_1,
        read2 = utils.original_read_2,
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
    "Runs fastp on a single-end sample."

    input: utils.original_read_1
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