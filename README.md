# el_gato
**E**pidemiology of ***L**egionella* : **G**enome-b**A**sed **T**yping:  

El_gato is a bioinformatics tool that utilizes either Illumina paired-end reads (.fastq) or a genome assembly (.fasta) as input to derive *Legionella pneumophila* Sequence Type (ST) from a database in contrast to the original method which relied on Sanger seqeunces. 

ST is used to describe relatedness of *L. pneumophila* isolates. The sequence of a portion of seven *L. pneumophila* genes (*flaA*, *pilE*, *asd*, *mip*, *mompS*, *proA*, and *neuA/neuAh*) is compared to a curated database of alleles and STs maintained by the European Society of Clinical Microbiology and Infectious Diseases Study Group for Legionella Infections (ESGLI) in which each unique allele is denoted with an allele number. The combination of allele numbers for all seven genes reported in order, corresponds to an allelic profile. The allelic profile, in turn, denotes a unique ST. 

* [Installation](#installation)
   * [Method 1: Using Conda](#method-1-using-conda)
   * [Method 2: Using pip](#method-2-using-pip)
     * [Dependencies](#dependencies)
* [Usage](#usage)
   * [Quickstart guide](#quickstart-guide)
   * [All available arguments](#all-available-arguments)
* [Acknowledgements](#acknowledgements)
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
# Create an environment, here named elgato, and install el_gato.py
# along with all dependencies
conda create -n elgato -c bioconda -c conda-forge el_gato

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
* [minimap2](https://github.com/lh3/minimap2) (>=2.24)
* [SAMTools](https://github.com/samtools/samtools) (>=1.15.1)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (>=2.13.0)
* [isPcr](https://users.soe.ucsc.edu/~kent/) (>=33)
* [fpdf2](https://pyfpdf.github.io/fpdf2/) (>=2.7.8)

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

```
Legionella in silico SBT script. 
    Requires paired-end reads files (preferred) or a genome assembly.

    Notes on arguments:
    (1) If only reads are provided, SBT is called using a mapping/alignment approach.
    (2) If only an assembly is provided, a BLAST and in silico PCR based approach is adopted. 

Input files:
  Please specify either reads files and/or a genome assembly file

  --read1 Read 1 file, -1 Read 1 file
                        Input Read 1 (forward) file
  --read2 Read 2 file, -2 Read 2 file
                        Input Read 2 (reverse) file
  --assembly Assembly file, -a Assembly file
                        Input assembly fasta file

Optional arguments:
  --help, -h            Show this help message and exit
  --version, -v         Print the version
  --threads THREADS, -t THREADS
                        Number of threads to run the programs (default: 1)
  --depth DEPTH, -d DEPTH
                        Specify the minimum depth used to identify loci in paired-end reads (default: 10)
  --kmer-size KMER_SIZE, -k KMER_SIZE
                        Specify the kmer sized used for mapping by minimap2. Max acceptable: 28. (default: 21)
  --out OUT, -o OUT     Output folder name (default: out)
  --sample SAMPLE, -n SAMPLE
                        Sample name (default: <Inferred from input file>)
  --overwrite, -w       Overwrite output directory (default: False)
  --sbt SBT, -s SBT     Database containing SBT allele and ST mapping files (default: /scicomp/home-pure/ptx4/el_gato/el_gato/db)
  --profile PROFILE, -p PROFILE
                        Name of allele profile to ST mapping file (default: /scicomp/home-pure/ptx4/el_gato/el_gato/db/lpneumophila.txt)
  --verbose             Print what the script is doing (default: False)
  --header, -e          Include column headers in the output table (default: False)
  --length LENGTH, -l LENGTH
                        Specify the BLAST hit length threshold for identifying multiple loci in assembly (default: 0.3)
  --sequence SEQUENCE, -q SEQUENCE
                        Specify the BLAST hit percent identity threshold for identifying multiple loci in assembly (default: 95.0)
  --samfile, -m         Specify whether or not the SAM file is included in the output directory (default: False)
```
# Acknowledgements

We greatly appreciate the United Kingdom Health Security Agency (UKHSA) for curating and sharing the *L. pneumophila* database. You can learn more about UKHSA here: https://www.gov.uk/government/organisations/uk-health-security-agency. Please contact legionella-sbt@ukhsa.gov.uk for enquiries about the database. 

 