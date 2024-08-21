from scripts import utils

configfile: "config.yaml"
GENOME_DIR = Path(config["genome_dir"])
RESULTS = Path(config["results"])
samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])
utils.samples = samples
utils.config = config
 
# TODO; get this from config, maybe extra file?
filter_genes = ["chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr14_GL000009v2_random", "chr14_GL000194v1_random", "chr14_KI270722v1_random", "chr14_KI270726v1_random", "chr15", "chr15_KI270727v1_random", "chr16", "chr16_KI270728v1_random", "chr17", "chr17_GL000205v2_random", "chr18", "chr19", "chr1_KI270706v1_random", "chr1_KI270708v1_random", "chr1_KI270711v1_random", "chr1_KI270712v1_random", "chr1_KI270713v1_random", "chr2", "chr20", "chr21", "chr22", "chr22_KI270731v1_random", "chr22_KI270733v1_random", "chr3", "chr3_GL000221v1_random", "chr4", "chr4_GL000008v2_random", "chr5", "chr6", "chr7", "chr8", "chr9", "chr9_KI270718v1_random", "chr9_KI270719v1_random", "chr9_KI270720v1_random", "chrM", "chrUn_GL000195v1", "chrUn_GL000213v1", "chrUn_GL000214v1", "chrUn_GL000218v1", "chrUn_GL000219v1", "chrUn_GL000220v1", "chrUn_GL000224v1", "chrUn_KI270442v1", "chrUn_KI270741v1", "chrUn_KI270742v1", "chrUn_KI270743v1", "chrUn_KI270744v1", "chrUn_KI270745v1", "chrUn_KI270746v1", "chrUn_KI270748v1", "chrUn_KI270750v1", "chrUn_KI270751v1", "chrUn_KI270754v1", "chrUn_KI270755v1", "chrX", "chrY"]


rule generate_trackhub:
    "Generates a trackhub for UCSC."

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
    "Converts the filtered and sorted bedGraph into a bigwig, for UCSC to display."

    input:
        bedgraph = f"{RESULTS}/{{sample}}_sorted.bg",
        chrom_sizes = f"{GENOME_DIR}/hg38.fa.sizes"
    output: bigwig = f"{RESULTS}/trackhub/hg38/{{sample}}.filtered.bigwig",
    params: trackhub = f"{RESULTS}/trackhub/hg38/"
    log: f"{RESULTS}/logs/bedgraph_to_bigwig/{{sample}}.log"
    shell: "bedGraphToBigWig {input.bedgraph} {input.chrom_sizes} {output.bigwig} > {log} 2>&1"


rule filter_and_sort_bedgraph:
    "Removes data from the spike in genomes and sorts + unpacks the resulting bedGraph."

    input: f"{RESULTS}/methylation/{{sample}}_pe.deduplicated.bedGraph.gz" #TODO; make .filtered a wildcard
    output: f"{RESULTS}/{{sample}}_sorted.bg"
    params: filter = "\|".join(filter_genes)
    log: f"{RESULTS}/logs/sort_bedgraph/{{sample}}.log"