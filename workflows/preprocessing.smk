configfile: "config.yaml"
WORKING_ROOT = Path(config["results"])
WORKING_DIR = WORKING_ROOT / "preprocessing"
samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])


rule all_preprocessing:
    input: expand(WORKING_DIR / "{sample}.json", sample=samples.alias)


rule fastp_pe:
    "Runs fastp on a paired-end sample."

    input:
        read1 = utils.original_read_1,
        read2 = utils.original_read_2,
    output:
        json = WORKING_DIR / "{sample}.json",  # TODO: make this temp?
        out1 = WORKING_DIR / "{sample}_R1.trimmed.fastq.gz",
        out2 = WORKING_DIR / "{sample}_R2.trimmed.fastq.gz",
    log: 
        WORKING_ROOT / "logs/fastp/{sample}.log"
    threads: 2
    shell:
        """
        fastp -w {threads} \
              --in1 {input.read1} \
              --in2 {input.read2} \
              --out1 {output.out1} \
              --out2 {output.out2} \
              -h /dev/null \
              -j {output.json} > {log} 2>&1
        """


rule fastp_se:
    "Runs fastp on a single-end sample."

    input: utils.original_read_1
    output:
        json = WORKING_DIR / "{sample}.json",  # TODO: make this temp?
        out = WORKING_DIR / "{sample}.trimmed.fastq.gz",
    log: 
        WORKING_ROOT / "logs/fastp/{sample}.log"
    wildcard_constraints: 
        sample = "(?!.*_R\d)"  # Ensures the sample does not contain _R1 or _R2
    threads: 2
    shell:
        """
        fastp -w {threads} \
              --in1 {input} \
              --out1 {output.out} \
              -h /dev/null \
              -j {output.json} > {log} 2>&1
        """
