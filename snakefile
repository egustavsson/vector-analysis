from snakemake.utils import min_version, validate

min_version("7.18")

# ----------------------------------------------------------------

configfile: "config.yml"
validate(config, schema="schema/config_schema.yaml")
WORKDIR = config["workdir"]

sample = config["sample_name"]

target_list = [
    "results/" + sample + "_minimap.sam"
]

rule all:
    input: 
        target_list

# ----------------------------------------------------------------

rule bam_to_fasta:
    input:
        bam = config["CCS_bam"]

    output:
        fa = "results/" + sample + ".fasta"
    
    threads: config["threads"]

    log: "logs/sample_bam_to_fasta.log"

    shell: """
        bam2fasta -u -o {output.fa} {input.bam} 2> {log}
    """

# ----------------------------------------------------------------

rule mapping:
    input:
        fa = rules.bam_to_fasta.output.fa,
        genome = config["genome"]

    output:
        sam = "results/" + sample + "_minimap.sam"
    
    params:
        opts=config["minimap2_opts"]

    log: "logs/sample_minimap.log"

    threads: config["threads"]

    shell: """
        minimap2 {params.opts} -t {threads} {input.genome} {input.fa} |
        samtools sort -n -@ {threads} -O SAM > {output.sam}
    """

# ----------------------------------------------------------------


# Sort this way to be able to index
# samtools sort -O SAM -o ./results/sorted.sam -@ 100 ./results/test_minimap.sam