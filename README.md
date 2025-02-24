# el_gato
**E**pidemiology of ***L**egionella* : **G**enome-b**A**sed **T**yping:  

El_gato is a bioinformatics tool that utilizes either Illumina paired-end reads (.fastq) or a genome assembly (.fasta) as input to replicate *Legionella pneumophila* Sanger-based Sequence Based Typing (SBT). From the input, seven loci (*flaA*, *pilE*, *asd*, *mip*, *mompS*, *proA*, and *neuA/neuAh*) are identified and compared to a [database of sequence types](https://github.com/appliedbinf/el_gato/blob/main/el_gato/db/lpneumophila.txt). The unique combination of the allelic identities of the seven target loci determines the sequence type for each input sample. 

* [Installation](#installation)
   * [Method 1: Using Conda](#method-1-using-conda)
   * [Method 2: Using pip](#method-2-using-pip)
     * [Dependencies](#dependencies)
* [Usage](#usage)
   * [Quickstart guide](#quickstart-guide)
   * [All available arguments](#all-available-arguments)
* [Input and Output](docs/input_output.md)
  * [Input files](docs/input_output.md/#input-files)
  * [Output files](docs/input_output.md/#output-files)
* [How does el_gato work?](docs/approach.md)
* [Reporting Module](docs/reporting_module.md)
* [Notices](docs/notices.md)

**Codebase stage:** XXX  
**Developers, maintainers, and testers:** [Alan Collins](https://github.com/Alan-Collins), [Will Overholt](https://github.com/waoverholt/), [Jenna Hamlin](https://github.com/jennahamlin)  

**Previous developers and maintainers:** [Dev Mashruwala](https://github.com/dmashruwala), [Andrew Conley](https://github.com/abconley), [Lavanya Rishishwar](https://github.com/lavanyarishishwar), [Emily T. Norris](https://github.com/norriset), [Anna Gaines](https://github.com/annagaines)

# Installation 

## Method 1: Using Conda
```
# Create an environment, here named elgato, and install el_gato.py plus all dependencies
conda create -n elgato -c bioconda -c conda-forge el_gato

[XX JH confirm final time that this will make el_gato and report module work XX]
[XX maybe but mamba in parenthesis after conda or vice versus b/c we use mamba XX]

# Activate the environment to use el_gato.py
conda activate elgato
```

## Method 2: Using pip
**Note** Using this method requires you to install all [Dependencies](#dependencies)
```
# Download el_gato by cloning the git repository
git clone https://github.com/appliedbinf/el_gato.git

# Move into the el_gato directory and install with pip
cd el_gato/
python3 -m pip install .
```
### Dependencies
* [python3](https://www.python.org/downloads/) (>=3.8,<3.12)
* [minimap2](https://github.com/lh3/minimap2) (2.28-r1209)
* [SAMTools](https://github.com/samtools/samtools) (1.20)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (2.16.0+)
* [isPcr](https://users.soe.ucsc.edu/~kent/) (33x2)
* [fpdf2]() [what version]
  
# Usage

## Quickstart Guide
An example of a basic run using paired-end reads or assemblies as input. **We recommend using reads whenever available**, as read-based sequence typing is more reliable ([see input for more information](docs/input_output.md)).
```
# Paired-end:
el_gato.py --read1 read1.fastq.gz --read2 read2.fastq.gz --out output_folder/

# Assembly:
el_gato.py --assembly assembly_file.fna --out output_folder/
```

## All available arguments
[XX add in new arguments like the kmer option once included XX]
Usage information printed when running el_gato.py with `-h` or `--help`.
```
usage: el_gato.py [--read1 Read 1 file] [--read2 Read 2 file] [--assembly Assembly file] [--help]  
[--threads THREADS] [--depth DEPTH]  [--out OUT] [--sample SAMPLE] [--overwrite] [--sbt SBT]  
[--suffix SUFFIX] [--profile PROFILE] [--verbose] [--header] [--length LENGTH] [--sequence SEQUENCE]
[--samfile SAMFILE]

Legionella in silico sequence-based typing (SBT) script.
    Requires paired-end reads files or a genome assembly.

    Notes on arguments:
    (1) If only reads are provided, SBT is called using a mapping/alignment approach and BLAST.
    (2) If only an assembly is provided, a BLAST and *in silico* PCR-based approach is adopted.

Input files:
  Please specify either reads files and/or a genome assembly file

  --read1 Read 1 file      -1 Read 1 file
                            Input Read 1 (forward) file
  --read2 Read 2 file      -2 Read 2 file
                            Input Read 2 (reverse) file
  --assembly Assembly file -a Assembly file
                            Input assembly fasta file

Optional arguments:
  --help,                  -h  
                            Show this help message and exit
  --threads THREADS        -t THREADS
                            Number of threads to run the programs (default: 1)
  --depth DEPTH            -d DEPTH
                            Variant read depth cutoff (default: 10)
  --out OUT                -o OUT  
                            Output folder name (default: out)
  --sample SAMPLE          -n SAMPLE
                            Sample name (default: <Inferred from input file>)
  --overwrite              -w  
                            Overwrite output directory(default: False)
  --sbt SBT                -s SBT  
                            Database containing SBT allele and ST mapping files  
                            (default: .../el_gato/el_gato/db)
  --suffix SUFFIX          -x SUFFIX
                            Suffix of SBT allele files (default: _alleles.tfa)
  --profile PROFILE        -p PROFILE  
                            Name of allele profile to ST mapping file  
                            (default: ../el_gato/el_gato/db/lpneumophila.txt) 
  --verbose                -v  
                            Print what the script is doing (default: False)
  --header                 -e  
                            Include column headers in the output table (default: False)
--length LENGTH            -l LENGTH  
                            Specify the BLAST hit length threshold for  
                            identifying multiple loci in an assembly  
                            (default: 0.3)
--sequence SEQUENCE        -q SEQUENCE  
                            Specify BLAST hit percent identity threshold for  
                            identifying multiple loci in an assembly  
                            (default: 95.0)
--samfile SAMFILE          -m SAMFILE
                            Allows the user to include the reads_vs_all_ref_filt.sam
                            file to be included in the output directory
                            (default: False)
```
