# Vector stability analysis

This reposotory provides the code and describes the analysis steps for assesing vector stability and integration using long-read sequencing.

## Installation

### Depedencies

- [miniconda](https://conda.io/miniconda.html)
- The rest of the dependencies (including `snakemake`) are installed via conda through the `environment.yml` file

### Installation process

Clone the directory:

```bash
git clone --recursive https://github.com/egustavsson/vector-stability.git
```

Create conda environment for the pipeline which will install all the dependencies:

```bash
cd vector-stability
conda env create -f environment.yml
```

## Input

- PacBio CCS reads as unmapped BAM
- Reference genome assembly in FASTA format. Described in the tutorial.

## How to use

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, make sure to first activate the conda environment using the command `conda activate vector-stability`.

```bash
cd vector-stability
conda activate vector-stability
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






## Long-read sequencing

## Analysis

1. Preparing the genome and annotation file
This section follows what has previously been described for the [AAV](https://github.com/Magdoll/AAV) analysis by Elizabeth Tseng ([Magdoll](https://github.com/Magdoll)).
The following fasta files (if available) should be combined into a single "genome" fasta file:

(required) vector (including the AAV vector + plasmid backbone as a single sequence)
host genome (ex: hg38)
NOTE the sequence IDs should be free of blank spaces and symbols. Stick with numbers, alphabet letters, and _ and -. If necessary, rename the sequence IDs in the combined fasta file.

Create a annotation.txt file according to the following format:

NAME=<sequence id>;TYPE={vector|helper|repcap|host|lambda};REGION=<start>-<end>;
Only the vector annotation is required and must be marked with REGION= (the position from ITR to ITR) as well. All other types are optional.

For example:

NAME=myAAV_plasmid;TYPE=vector;REGION=100-2000;
NAME=myRepCap_plasmid;TYPE=repcap;REGION=500-1500;
NAME=myHelper_plasmid;TYPE=helper;
NAME=chr1;TYPE=host;
NAME=chr2;TYPE=host;
NAME=chr3;TYPE=host;
IMPORTANT!!! you must have exactly the same number of chromosomes in the reference fasta file as annotations file. This is especially common if you are including human genome (hg38) which has a lot of alternative chromosomes. It is recommended that you use a version of hg38 that only lists the major chromosomes.

- paste file1.fasta file2.fasta > combined.fasta
