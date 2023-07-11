# Vector stability analysis

This reposotory provides the code and describes the analysis steps for assesing vector stability and integration using long-read sequencing.

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
