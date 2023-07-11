# Vector stability analysis

This reposotory provides the code and describes the analysis steps for assesing vector stability and integration using long-read sequencing. The first step is to create the appropriate reference for mapping that incudes the vector sequence. The rest of the steps are performed by a `snakemake` pipeline. Description of the steps are found [here](#analysis-steps). This pipeline borrows from the [AAV](https://github.com/Magdoll/AAV) analysis by Elizabeth Tseng ([Magdoll](https://github.com/Magdoll)).

## Pipeline 

### Installation

#### Depedencies

- [miniconda](https://conda.io/miniconda.html)
- The rest of the dependencies (including `snakemake`) are installed via conda through the `environment.yml` file

#### Installation process

Clone the directory:

```bash
git clone --recursive https://github.com/egustavsson/vector-analysis.git
```

Create conda environment for the pipeline which will install all the dependencies:

```bash
cd vector-analysis
conda env create -f environment.yml
```

### Input

- PacBio CCS reads as unmapped BAM
- Reference genome assembly in FASTA format. Described in the analysis [tutorial](#1-preparing-the-genome-and-annotation-file).

### How to use

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, make sure to first activate the conda environment using the command `conda activate vector-analysis`.

```bash
cd vector-analysis
conda activate vector-analysis
snakemake --use-conda -j <num_cores> all
```
It is a good idea to do a dry run (using -n parameter) to view what would be done by the pipeline before executing the pipeline.

```bash
snakemake --use-conda -n all
```

To exit a running `snakemake` pipeline, hit `ctrl+c` on the terminal. If the pipeline is running in the background, you can send a `TERM` signal which will stop the scheduling of new jobs and wait for all running jobs to be finished.

```bash
killall -TERM snakemake
```

To deactivate the conda environment:
```bash
conda deactivate
```

## Analysis steps

### 1. Preparing the genome and annotation file
The following fasta files (if available) should be combined into a single "genome" fasta file:
- host genome (ex: hg38). Can be downloaded from [Ensembl](https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/), [GENCODE](https://www.gencodegenes.org/human/) or [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/)
- vector (including the vector + plasmid backbone as a single sequence)

*NOTE* the sequence IDs should be free of blank spaces and symbols. Stick with numbers, alphabet letters, and _ and -. If necessary, rename the sequence IDs in the combined fasta file.

Create a annotation.txt file according to the following format:
```
NAME=<sequence id>;TYPE={vector|helper|repcap|host|lambda};REGION=<start>-<end>;
```
Only the `vector` annotation is required and must be marked with `REGION=`. All other types are optional. For example:

```
NAME=my_plasmid;TYPE=vector;REGION=100-2000;
NAME=chr1;TYPE=host;
NAME=chr2;TYPE=host;
NAME=chr3;TYPE=host;
```

**IMPORTANT!!!** you must have exactly the same number of chromosomes in the reference fasta file as annotations file. This is especially common if you are including human genome (hg38) which has a lot of alternative chromosomes. It is recommended that you use a version of hg38 that only lists the major chromosomes.

Combinig the host genome and the vector into a single fasta can be done by: 

```bash
paste host.fasta vector.fasta > combined.fasta
```

### 2. Running the snakemake
After generated required genomes used for mapping are generated, make sure the `config.yml` is edited. These are the parameters:

| Parameter | Description |
| --- | --- |
| workdir | This is the working directory and should be set to `./vector-analysis/` |
| sample_name | sample name which will be the prefix of output. Default is `Sample` |
| genome | genome fasta that will be used to map against |
| CCS_bam | The unmapped HiFi CCS bam file |
| minimap2_opts | options passed to minimap2 for mapping. Default is `--eqx -a --secondary=no` |
| threads | Number of threads to use. Default is 10 |