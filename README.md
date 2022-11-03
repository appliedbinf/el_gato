# el_gato
Epidemiology of *Legionella*: Genome-bAsed Typing

* [Installation](#installation)
   * [Method 1: using conda ](#method-1:-using-conda)
   * [Method 2: using pip](#method-2:-using-pip)
* [Quickstart guide](#Quickstart-guide)
   * [A typical run](#A-typical-run)      
* [All available arguments](#Usage)
* [Dependencies](#Dependencies)

Currently in development
Codebase stage: development   
Developers and maintainers, Testers: [Andrew Conley](https://github.com/abconley), [Lavanya Rishishwar](https://github.com/lavanyarishishwar), [Will Overholt](https://github.com/), [Dev Mashruwala](https://github.com/), [Alan Collins](https://github.com/Alan-Collins)

# Installation 

## Method 1: using conda

```
# Create environment named elgato and install el_gato.py plus all dependencies
conda create -n elgato -c bioconda -c appliedbinf elgato

# Activate the environment to use el_gato.py
conda activate elgato
```

## Method 2: using pip

**N.B.** Using this method requires you to manually install all [dependencies](#Dependencies).

```
# Download el_gato by cloning the git repository
git clone https://github.com/appliedbinf/el_gato.git

# Move into the el_gato directory and install with pip
cd el_gato/
python3 -m pip install .
```

# Quickstart Guide

## A typical run

Here is an example of a basic run using paired end reads, assemblies, or both as input.

```
# Paired-end:
el_gato.py --read1 read1.fastq.gz --read2 read2.fastq.gz --out output_folder/

# Assembly:
el_gato.py --assembly assembly_file.fna --out output_folder/

# Combined:
el_gato.py --read1 read1.fastq.gz --read2 read2.fastq.gz --assembly assembly_file.fna --out output_folder/

```
# Usage

```
Required arguments:  
--read1 (paired end read1)  
--read2 (paired end read2)  
--assembly (assembly file)  

Optional arguments:   
--out (output folder name)  
--spades, -g (Runs SPAdes when paired-end samples are given (default: False))  
--help, -h (help)  
--threads, -t (threads (default: 1))  
--sample, -n (sample name)    
--overwrite, -w (overwrites output folder name)   
--sbt, -s (database containing SBT allele and mapping files))   
--suffix, -x (suffix of SBT allele files (default: _alleles.tfa))  
--profile, -p (name of allele profile in ST mapping file)   
--verbose -v (print what the script is doing (default: False))    
```

# Dependencies

[bwa](https://github.com/lh3/bwa)
[sambamba](https://github.com/biod/sambamba)
[freebayes](https://github.com/ekg/freebayes)
[SAMTools](https://github.com/samtools/samtools)
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
[isPcr](https://users.soe.ucsc.edu/~kent/)
[SPAdes](http://cab.spbu.ru/software/spades/)
[stringMLST](https://github.com/jordanlab/stringMLST)
