from os import path
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("7.18")

# ----------------------------------------------------------------

configfile: "config.yml"
validate(config, schema="schema/config_schema.yaml")
workdir: path.join(config["workdir"], config["pipeline"])

WORKDIR = path.join(config["workdir"], config["pipeline"])
SNAKEDIR = path.dirname(workflow.snakefile)

sample = config["sample_name"]

target_list = [
    "results/" + sample + ".bam"
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
        fa = "results/" + sample + ".fasta" + ".fasta",
        genome = config["genome"]

    output:
        sam = "results/" + sample + "_minimap.sam"
    
    log: "logs/sample_minimap.log"

    threads: config["threads"]

    shell: """
        minimap2 --eqx -a --secondary=no -t {threads} {input.genome} {input.fa} -o {output.sam} 2> {log}
    """

# ----------------------------------------------------------------

rule sam_to_bam:
    input: 
        sam = "results/" + sample + "_minimap.sam"
    
    output:
        bam = "results/" + sample + ".bam"
    
    log: "logs/sample_samtools.log"

    threads: config["threads"]

    shell: """
        samtools sort -O BAM -o {output.bam} -@ {threads} {input.sam};
        samtools index {output.bam}

    """