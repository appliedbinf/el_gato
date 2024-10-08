# el_gato
**E**pidemiology of ***L**egionella* : **G**enome-b**A**sed **T**yping:  

El_gato is a bioinformatics tool that utilizes either a genome assembly (.fasta) or Illumina paired-end reads (.fastq) to replicate *Legionella pneumophila* Sequence Based Typing (SBT). From the input, 7 loci (*flaA*, *pilE*, *asd*, *mip*, *mompS*, *proA*, *neuA/neuAh*) are identified and compared to a database of sequence types. The sequence type provided for each input sample is based on the unique combination of the allelic identities of the 7 target loci. 

* [Installation](#installation)
   * [Method 1: using conda](#method-1-using-conda)
   * [Method 2: using pip](#method-2-using-pip)
     * [Dependencies](#dependencies)
* [Usage](#usage)
   * [Quickstart guide](#quickstart-guide)
   * [All available arguments](#all-available-arguments)
* [Input and Output](#input-and-output)
  * [Input files](#input-files)
     * [Paired-end reads](#pair-end-reads)
     * [Genome assemblies](#genome-assemblies)      
  * [Output files](#output-files)
     * [standard out](#standard-out)
     * [possible_mlsts.txt](#possible_mlststxt)
     * [intermediate_outputs.txt](#intermediate_outputstxt)
     * [identified_alleles.fna](#identified_allelesfna)
     * [run.log](#runlog)
     * [reads_vs_all_ref_filt_sorted.bam](#reads_vs_all_ref_filt_sortedbam-reads-only)
     * [report.json](#reportjson)
* [How does el_gato work?](#approach)
* [Using Nextflow](#using-nextflow)
* [Reporting Module](#reporting-module)

Codebase stage: development   
Developers and maintainers, Testers: [Alan Collins](https://github.com/Alan-Collins), [Dev Mashruwala](https://github.com/dmashruwala), [Andrew Conley](https://github.com/abconley), [Lavanya Rishishwar](https://github.com/lavanyarishishwar), [Emily T. Norris](https://github.com/norriset), [Anna Gaines](https://github.com/annagaines), [Will Overholt](https://github.com/waoverholt/)

# Installation 

## Method 1: using conda
```
# Create an environment named elgato and install el_gato.py plus all dependencies
conda create -n elgato -c bioconda -c conda-forge el_gato

# Activate the environment to use el_gato.py
conda activate elgato
```

## Method 2: using pip
**Note** Using this method requires you to install all [dependencies](#dependencies) manually.
```
# Download el_gato by cloning the git repository
git clone https://github.com/appliedbinf/el_gato.git

# Move into the el_gato directory and install with pip
cd el_gato/
python3 -m pip install .
```
### Dependencies
* [minimap2](https://github.com/lh3/minimap2)
* [SAMTools](https://github.com/samtools/samtools)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [isPcr](https://users.soe.ucsc.edu/~kent/)
  
# Usage

## Quickstart Guide
Here is an example of a basic run using paired-end reads or assemblies as input. We recommend using reads whenever available as read-based sequence typing is more reliable.
```
# Paired-end:
el_gato.py --read1 read1.fastq.gz --read2 read2.fastq.gz --out output_folder/

# Assembly:
el_gato.py --assembly assembly_file.fna --out output_folder/
```

## All available arguments
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
--sequence SEQUENCE        -s SEQUENCE  
                            Specify BLAST hit percent identity threshold for  
                            identifying multiple loci in an assembly  
                            (default: 95.0)
--samfile SAMFILE          -m SAMFILE
                            Allows the user to include the reads_vs_all_ref_filt.sam
                            file to be included in the output directory
                            (default: False)
```

# Input and Output

## Input files

If available, we recommend using raw or trimmed reads instead of assemblies as the extra data contained in reads is valuable for the process used by el_gato to identify sample ST. When run with reads, el_gato is able to use read quality and coverage information to apply quality control rules. When run using assemblies, el_gato is unable to identify errors that have been incorporated into the assembly and may therefore report incorrect results. For example, while many isolates encode two copies of *mompS*, in some cases only one copy of the locus is assembled. If only the secondary *mompS* locus is included in the assembly then el_gato will report that allele.

#### Pair-end reads
When running on a directory of reads, files are associated as pairs using the pattern `R{1,2}.fastq`. i.e., filenames should be identical except for containing either "R1" or "R2" and can be .fastq or .fastq.gz format. Any files for which a pair can not be identified using this pattern will not be processed.

#### Genome assemblies
When running on a directory of assemblies, all files in the target directory will be processed, and there is no filename restrictions.

## Output files
After a run, el_gato will print the identified ST of your sample to your terminal ([stdout](#standard-out)) and write several files to the specified output directory (default: out/).  A subdirectory is created for each sample processed, and each subdirectory is named with its sample name and contains output files specific to that sample (see below). 

### The files included in the output directory for a sample are: 

### standard out
ST profile is written as a tab-delimited table without the headings. Headings are included if el_gato.py is run with `-e` flag and are displayed like so:

`Sample  ST flaA  pilE  asd   mip   mompS proA  neuA_neuAH`    

 The sample column contains the user-provided or inferred sample name. The ST column contains the overall sequence type of the sample. 

The ST column can contain two kinds of values. If the identified ST corresponds to a profile found in the database, the corresponding number is given. If no matching ST profile is found or el_gato was unable to make a confident call, then this will be reflected in the value displayed in the ST column.

The corresponding allele number is reported for each gene if an exact allele match is found in the database. Alternatively, el_gato may also note the following symbols:

| Symbol | Meaning |
|:------:|:---------|
|Novel ST    | Novel Sequence Type: All 7 target genes were found, but not present in the profile - most likely a novel sequence type. |
|Novel ST*    | Novel Sequence Type due to novel allele: One or multiple target genes have a novel allele found. |
|MD-     | Missing Data: ST is  unidentifiable as a result of or more of the target genes that are unidentifiable.  |
|MA?     | Multiple Alleles: ST is ambiguous due to multiple alleles that could not be resolved. |
| NAT    | Novel Allele Type: BLAST cannot find an exact allele match - most likely a new allele. |
| -      | Missing Data: Both percent and length identities are too low to return a match or N's in sequence. |
| ?      | Multiple Alleles: More than one allele was found and could not be resolved. |

If symbols are present in the ST profile, the other output files produced by el_gato will provide additional information to understand what is being communicated.

### possible_mlsts.txt
This file would contain all possible ST profiles if el_gato identified multiple possible alleles for any ST loci. In addition, if multiple *mompS* alleles were found, the information used to determine the primary allele is reported in two columns: "mompS_reads_support" and "mompS_reads_against." mompS_reads_support indicates the number of reads associated with each allele that contains the reverse sequencing primer in the expected orientation, which suggests that this is the primary allele. mompS_reads_against indicates the number of reads containing the reverse sequencing primer in the wrong orientation and thus demonstrates that this is the secondary allele. These values are used to infer which allele is the primary *mompS* allele, and their values can be considered to represent the confidence of this characterization. [See Approach subsection for more details](#reads).

### intermediate_outputs.txt
el_gato calls other programs to perform intermediate analyses. The outputs of those programs are provided in this file. In addition, to help with troubleshooting issues, important log messages are also written in this file. The following information may be contained in this file, depending on if the input is reads or assembly:

* Reads-only - Samtools coverage command output. [See samtools coverage documentation for more information about headers](https://www.htslib.org/doc/samtools-coverage.html) or [below.](#samtools-coverage-headers)
* Reads-only - Information about the orientation of *mompS* sequencing primer in reads mapping to biallelic sites. [See Approach subsection for more details](#reads).
* BLAST output indicating the best match for identified alleles. [See BLAST output documentation for more information about headers](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/) or [below.](#blastn-output-headers)

Headers are included in outputs for the samtools coverage command and blast results. Header definitions are as follows:

#### samtools coverage headers

| Column header | Meaning 
|:-------------:|:---------------------------------------------------:|
| rname         | Locus name                                          |
| numreads      | Number reads aligned to the region (after filtering)|
| covbases      | Number of covered bases with depth >= 10             |
| coverage      | Percentage of covered bases [0..100]                |
| meandepth     | Mean depth of coverage                              |
| meanbaseq     | Mean baseQ in covered region                        |
| meanmapq      | Mean mapQ of selected reads                         | 

### BLASTn output headers

| Column header | Meaning                             |
|:-------------:|:-----------------------------------:|
| qseqid        | Query sequence id                   |
| sseqid        | Subject (matched allele) id         |
| pident        | Percentage of identical matches     |
| length        | Alignment length (sequence overlap) |
| mismatch      | Number of mismatches                |
| gapopen       | Mumber of gap openings              |
| qstart        | Start of alignment in query         |
| qend          | End of alignment in query           |
| sstart        | Start of alignment in subject       |
| send          | End of alignment in subject         |
| evalue        | Expect value                        |
| bitscore      | Bit score                           |
| qlen          | Query sequence length               |
| slen          | Subject sequence length             |


### identified_alleles.fna
The nucleotide sequence of all identified alleles is written in this file. If more than one allele is determined for the same locus, they are numbered arbitrarily. Fasta headers of sequences in this file correspond to the query IDs in the BLAST output reported in the intermediate_outputs.txt file.

### run.log
A detailed log of the steps taken during el_gato's running includes the outputs of any programs called by el_gato and any errors encountered. Some command outputs include headers (e.g., samtools coverage and BLAST).

### reads_vs_all_ref_filt_sorted.bam (reads only)
el_gato maps the provided reads to [a set of reference sequences in the el_gato db directory](https://github.com/appliedbinf/el_gato/blob/main/el_gato/db/ref_gene_regions.fna). The mapped reads are then used to extract the sequences present in the sample for identifying the alleles and, ultimately, the ST. reads_vs_all_ref_filt_sorted.bam and its associated file reads_vs_all_ref_filt_sorted.bai contains the mapping information that was used by el_gato. The BAM file can be viewed using software such as [IGV](https://software.broadinstitute.org/software/igv/) to understand better the data used by el_gato to make allele calls. Additionally, this file is a good starting point for investigating the cause of incorrectly resolved loci.

**Note:** A SAM file is also present, which has the same information as in the BAM file.

### report-json
Each sample outputs a json file that contains relevant information about the run that will be included in the report PDF.   

Summary page metadata: Complete MLST profile of the sample and the abbreviation key for the symbols.  

Run-specific data:  

Paired-end reads: Locus coverage information and *mompS* primer information.  

Assembly: BLAST hit length and sequence identity thresholds and locus location information.  

# Approach

At its core, el_gato uses BLAST to identify the closest match to each allele in your input data. For the loci *flaA*, *pilE*, *asd*, *mip*, and *proA*, this process is straight forward. Whereas loci *mompS* and *neuA/neuAh* require more involved processing, with neuA/neuAh being an issue only when processing reads—the specifics of these loci are discussed in the corresponding sections below. 


First for the simple loci (*flaA*, *pilE*, *asd*, *mip*, and *proA*), the following processes are used:

## Reads

When processing reads, identification of both *mompS* and *neuA*/*neuAh* requires additional analyses (described below). The other five loci are processed by mapping the provided reads to reference loci from *L. pneumophila* strain Paris and identifying the consensus sequence. Then, all alleles are determined using BLAST against the SBT allele database.

A couple of quality control steps are applied when processing the reads:

   1. **Base quality:** Any bases with quality scores below 20 are not included when calculating coverage at each position or identifying alternate base calls. The lowest number of bases with quality over 20 that map to a single position is reported in the log for each locus [add exact name of the identifier].
   2. **Coverage:** After excluding low-quality bases, if there is not at least one read covering 100% of the locus (<99% for *neuA*/*neuAh* - see below), then no attempt to identify the allele is made, and a "-" will be reported. A minimum depth of 10 is applied as a cutoff. 

### *neuA/neuAh*

[The sequence of *neuA*/*neuAh* loci can differ dramatically.](https://doi.org/10.1111/1469-0691.12459) The differences in sequence between *neuA*/*neuAh* alleles are sufficient that reads from some alleles will not map to others. Accordingly, we map reads to [X number of] reference sequences that cover the sequence variation currently represented in the SBT database. The [X number of] reference alleles used are the *neuA* allele from strain Paris (neuA_1), the *neuAh* allele from strain Dallas-1E (neuA_201), and [description of other reference sequences]. The reference sequence with the best mapping [what defines best? Is it the number of reads?] is identified using `samtools coverage` with the caveat that >99% of the *neuA*/*neuAh* locus must have coverage of at least one read (some alleles contain small indels, so 100% is too strict); otherwise a "-" will be reported. Once the reference sequence is selected, the processing using BLAST is the same as described above.

### *mompS*

[As described in the assembly section](#assembly), *mompS* is sometimes present in multiple copies in the genome of *L. pneumophila* isolates though it is typically two copies. Duplicate gene copies pose an obvious challenge for a read-mapping approach: if two similar sequence copies are present in a genome, reads from both copies may map to the same reference sequence.

el_gato resolves this issue by taking advantage of the proximity of the genome's two copies of *mompS* [(a schematic of the organization of the two *mompS* copies can be found in Fig. 1 in this paper).](https://doi.org/10.1016/j.cmi.2017.01.002). The sequence context of the two *mompS* copies is such that the correct copy is immediately upstream of the incorrect copy. Only the right copy is flanked on either side by sequences corresponding to primers (*mompS*-450F and *mompS*-1116R) used for conventional SBT. In contrast, while *mompS*-450F sequences are present upstream of the incorrect copy, the corresponding *mompS*-1116R sequences are not found downstream. For this reason, only the correct copy is amplified in conventional SBT [(See below schematic.)](#momps-read-mapping-schematic)

The sequence of the two copies of *mompS* and the identity of the correct allele is then resolved through the following process:

1. Reads from both *mompS* copies are mapped to a single *mompS* reference sequence flanked by the *mompS*-450F and *mompS*-1116R primer sequences. 

2. The nucleotide sequence of reads is recorded for each position within the *mompS* sequence. If the base at a particular position is heterogeneous in more than 30% of reads mapped to that position, the position is considered biallelic, and both bases are recorded. If the sequence contains only one copy of *mompS* or no sequence variation between the duplicate copies, then no biallelic sites will be found. Then, the sequence will be extracted, and an allele can be identified using BLAST. 

3. If multiple biallelic positions are identified, all sequences are recorded, and individual read pairs are identified that map to each of the biallelic position(s). 

4. Once the two alleles have been resolved, the correct allele for SBT is identified by analyzing the reads associated with each allele. Reads from each allele are searched for the *mompS*-1116R sequence. The orientation of the reads that contain the primer sequence is assessed. If the primer maps 3'-5' relative to the reference sequence (i.e., in the reverse direction), this is consistent with the read pair originating from the correct copy of *mompS*. However, if the read containing the primer maps 5'-3' (i.e., in the forward direction) relative to *mompS*, this is consistent with the read pair originating from the wrong copy of *mompS*. 

5. The number of reads associated with each allele that contains the primer in the correct orientation relative to *mompS* is counted and compared. The correct allele is then chosen using the following criteria:  

   a. Only one allele has associated reads with primer *mompS*-1116R correctly oriented.  

   b. One allele has more than three times as many reads with correctly oriented primer as the other.  

   c. One allele has no associated reads with the primer *mompS*-1116R in either orientation, but the other allele has associated reads with the primer in only the wrong orientation. In this case, the allele with no associated reads with the primer in either orientation is considered the primary locus by the process of elimination.

6. The allele number of all alleles is then determined using BLAST and the ST is generated using the correct allele. 

Note that as the above process depends upon read pairs mapping to biallelic sites and the primer region, sequence data characteristics such as read length and insert size can impact the ability of el_gato to determine *mompS* alleles. 

If the above process cannot identify the correct sequence, a ? will be returned as the *mompS* allele, and el_gato will report information about the steps in this process in the [output files](#output-files).

## *momps* Read Mapping Schematic
![mompS read mapping schematic](https://github.com/appliedbinf/el_gato/blob/images/images/mompS_allele_assignment.png)

## Assembly

Six of the seven loci (*flaA*, *pilE*, *asd*, *mip*,*proA*, and *neuA/neuAh*) are identified using BLAST. For each, the best BLAST result is returned as the allele. The closest match is returned with an \* if loci have no exact match. [add in updated code changes if necessary] When processing an assembly, only *mompS* requires extra processing. 

### *mompS*

[*mompS* is sometimes present in multiple copies in *Legionella pneumophila*, though typically two copies.](https://doi.org/10.1016/j.cmi.2017.01.002) When typing *L. pneumophila* using Sanger sequencing, primers amplify only the correct *mompS* locus. We, therefore, use *in silico* PCR to extract the correct *mompS* locus sequence from the assembly. The primers used for *in silico* PCR are *mompS*-450F (TTGACCATGAGTGGGATTGG) and *mompS*-1116R (TGGATAAATTATCCAGCCGGACTTC) [as described in this protocol](https://doi.org/10.1007/978-1-62703-161-5_6). The *mompS* allele is then identified using BLAST.

# Using NextFlow 
## Using NextFlow with Singularity Container

We provide a singularity container that can be run using the nextflow workflow for el_gato on a directory of either reads or assemblies. In both cases the target directory must contain only paired reads files (in .fastq or .fastq.gz format) or assembly files (in fasta format).

```
# Reads
nextflow run_el_gato.nf --reads_dir <path/to/reads/directory> --threads <threads> --out <path/to/output/directory> -profile singularity -c nextflow.config

# Assemblies
nextflow run_el_gato.nf --assembly_dir <path/to/assemblies/directory> --threads <threads> --out <path/to/output/directory> -profile singularity -c nextflow.config
```

**Note:** To run nextflow without the singularity container, uncomment conda environment installation on line 10 and line 47 of the run_el_gato.nf file and use the following commands:

```
# Reads
nextflow run_el_gato.nf --reads_dir <path/to/reads/directory> --threads <threads> --out <path/to/output/directory>

# Assemblies
nextflow run_el_gato.nf --assembly_dir <path/to/assemblies/directory> --threads <threads> --out <path/to/output/directory>
```
## Output files for Nextflow
At the completion of a run, the specified output directory (default: el_gato_out/) will contain a file named "all_mlst.txt" (the MLST profile of each sample) and one directory for each sample processed. Each sub-directory is named with a sample name and contains output files specific to that sample. These files include the el_gato log file and files providing more details about the sequences identified in the sample.  

Additionally, the specified output directory will contain a combined json file (report.json) that contains all of the data from the individual sample-level json files and the report PDF (report.pdf).

# Reporting Module  
## Dependencies
  * [fpdf2](https://github.com/py-pdf/fpdf2)

We provide a script that generates a PDF report of each el_gato run using the report.json file generated in the output folder for each sample.  
This report generated by Nextflow by default, but must be run manually if running el_gato for individual samples. 
A report can be generated for one or more samples and may include both assembly and read report jsons.
```
elgato_report.py -i <path/to/report1.json> [<path/to/report2.json> ...] -o <path/to/output/report.pdf>
```
