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
    "results/" + sample + ".bam",
    "results/" + sample + ".vcf"
]

rule all:
    input: 
        target_list

# ----------------------------------------------------------------

rule mapping:
    input:
        fa = config["CCS_fasta"],
        genome = config["genome"]

    output:
        sam = "results/" + sample + "_minimap.sam"
    
    params:
        mm_opts = config["minimap_opts"]
    
    log: "logs/" + sample + "_minimap.log"

    threads: config["threads"]

    shell: """
        minimap2 --eqx -a {params.mm_opts} --secondary=no -t {threads} {input.genome} {input.fa} -o {output.sam} 2> {log}
    """

# ----------------------------------------------------------------

rule sam_to_bam:
    input: 
        sam = "results/" + sample + "_minimap.sam"
    
    output:
        bam = "results/" + sample + ".bam"
    
    log: "logs/" + sample + "_samtools.log"

    threads: config["threads"]

    shell: """
        samtools sort -O BAM -o {output.bam} -@ {threads} {input.sam};
        samtools index {output.bam}
    """

# ----------------------------------------------------------------

rule sniffles:
    input:
       bam = "results/" + sample + ".bam"
    
    output:
        vcf = "results/" + sample + ".vcf"
    
    params: 
        sn_opts = config["sniffles_opts"]
    
    log: "logs/" + sample + "_sniffles.log"

    threads: config["threads"]

    shell: """
        sniffles -i {input.bam}  -v {output.vcf} {params.sn_opts} --threads {threads} 2> {log}
        """