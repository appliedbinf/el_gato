# el_gato
Epidemiology of *Legionella*: Genome-bAsed Typing

Currently in development

# Installation 
  1. Download or clone el_gato into your home directory
  2. Create a new conda environment to install dependencies:  
      conda create --name el_gato
  3. Activate your newly created conda environment:  
      conda activate el_gato
  4. Install the following dependencies via bioconda and conda-forge:  
      bwa, sambamba, freebayes, SAMTools, BLAST, isPcr, SPAdes, stringMLST

# Usage

Required arguments:  
--read1 (paired end read1)  
--read2 (paired end read2)  
--assembly (assembly file)  

Optional arguments:   
--out <<output folder name)  
--help, -h help  
--threads, -t threads (default: 1)  
--sample, -n sample name  
--overwrite, -w overwrites output folder name  
--sbt, -s database containing SBT allele and mapping files  
--suffix, -x suffix of SBT allele files (default: _alleles.tfa)  
--profile, -p name of allele profile in ST mapping file  
--verbose -v print what the script is doing (default: False)  

Paired-end:  
   python3 el_gato.py --read1 read1 --read2 read2 --out output folder  

Assembly:  
   python3 el_gato.py --assembly assembly file --out output folder name  

Combined:  
   python3 el_gato.py --read1 read1 --read2 read2 --assembly assembly file --out output folder
