# el_gato
Epidemiology of *Legionella*: Genome-bAsed Typing

* [Installation](#installation)
   * [Method 1. using conda yml script](#method-1-using-conda-yml-script)
   * [Method 2. Creating conda environment and loading the dependencies](#method-2-Creating-conda-environment-and-loading-the-dependencies)
   * [Quickstart guide](#Quickstart-guide)
      * [A typical run](#A-typical-run)
      * [General program structure overview](#General-program-structure-overview)
      * [Input file format](#Input-file-format)      
   * [All available arguments](#All-available-arguments)
   * [Literature Citations](#Literature-Citations)
   * [Software](#software)

Currently in development
Codebase stage: development   
Developers and maintainers, Testers: [Andrew Conley](https://github.com/abconley), [Lavanya Rishishwar](https://github.com/lavanyarishishwar), [Will Overholt](https://github.com/), [Dev Mashruwala](https://github.com/)

# Installation 

  1. Download or clone el_gato into your home directory
  2. Create a new conda environment to install dependencies:  
  ```
  # You can install the dependencies while creating the el_gato base environment.
    conda install mamba
    conda create -n el_gato 

  # Activate the conda environement using
    conda activa el_gato
    mamba install bwa, sambamba, freebayes, SAMTools, BLAST, isPcr, SPAdes, stringMLST
  ```

# Usage

Required arguments:  
--read1 (paired end read1)  
--read2 (paired end read2)  
--assembly (assembly file)  

Optional arguments:   
--out (output folder name)  
--help, -h (help)  
--threads, -t (threads (default: 1))  
--sample, -n (sample name)    
--overwrite, -w (overwrites output folder name)   
--sbt, -s (database containing SBT allele and mapping files))   
--suffix, -x (suffix of SBT allele files (default: _alleles.tfa))  
--profile, -p (name of allele profile in ST mapping file)   
--verbose -v (print what the script is doing (default: False))    

Paired-end:  
   python3 el_gato.py --read1 read1 --read2 read2 --out output folder  

Assembly:  
   python3 el_gato.py --assembly assembly file --out output folder name  

Combined:  
   python3 el_gato.py --read1 read1 --read2 read2 --assembly assembly file --out output folder
