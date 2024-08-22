configfile: "config.yaml"
GENOME_DIR = Path(config["genome_dir"])
WORKING_ROOT = Path(config["results"])
WORKING_DIR = WORKING_ROOT / "alignment"

RESULTS = Path(config["results"])
samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])
utils.samples = samples
utils.config = config


rule all_alignment:
    input: expand(WORKING_DIR / "{sample}_pe.deduplicated.bam", sample=samples.alias)


rule bismark_deduplicate:
    """
    Deduplicates the BAM file from the bismark_bowtie2 rule. 
    Removes all but one read that align to the exact same position & orientation.
    """

    input: 
        WORKING_DIR / "{sample}_pe.bam"
    output: 
        bam = WORKING_DIR / "{sample}_pe.deduplicated.bam",
        report = WORKING_DIR / "{sample}_pe.deduplication_report.txt"
    params: 
        WORKING_DIR
    log: 
        WORKING_ROOT / "logs/bismark_deduplicate/{sample}.log"
    threads: 1
    shell: 
        """
        deduplicate_bismark -p \
            --output_dir {params} \
            {input} > {log} 2>&1
        """


rule bismark_bowtie2:
    "Runs bismark using bowtie2 on the prepared genome."

    input:
        CT = expand(GENOME_DIR / "Bisulfite_Genome/CT_conversion/BS_CT.{n}.bt2", n=[1, 2, 3, 4]),
        GA = expand(GENOME_DIR / "Bisulfite_Genome/GA_conversion/BS_GA.{n}.bt2", n=[1, 2, 3, 4]),
        genome = GENOME_DIR,
        reads = utils.get_trimmed_reads
    output:
        bam = WORKING_DIR / "{sample}_pe.bam",
        report = WORKING_DIR / "{sample}_PE_report.txt"
    params:
        output_dir = WORKING_DIR,
        sample = lambda w: w.sample
    log: 
        WORKING_ROOT / "logs/bismark_bowtie2/{sample}.log"
    threads: 16
    shell: 
        """
        bismark -N 1 -p 8 --bowtie2 \
                --basename {params.sample} \
                --output_dir {params.output_dir} \
                --genome_folder {input.genome} \
                -1 {input.reads[0]} \
                -2 {input.reads[1]} > {log} 2>&1
        """


rule bismark_genome_preparation:
    "Prepares the genome to be used by Bismark."

    input: 
        GENOME_DIR
    output:
        CT = expand(f"{GENOME_DIR}/Bisulfite_Genome/CT_conversion/BS_CT.{{n}}.bt2", n=[1, 2, 3, 4]),
        GA = expand(f"{GENOME_DIR}/Bisulfite_Genome/GA_conversion/BS_GA.{{n}}.bt2", n=[1, 2, 3, 4]),
    log: 
        WORKING_ROOT / "logs/bismark_genome_preparation/genome_preparation.log"
    threads: 16
    shell: 
        """
        bismark_genome_preparation \
            --bowtie2 \
            --parallel {threads} \
            {input} > {log} 2>&1
        """
