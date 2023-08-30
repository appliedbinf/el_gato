# el_gato
**E**pidemiology of ***L**egionella* : **G**enome-b**A**sed **T**yping:  

[Add two sentence summary of what el gato is]

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
* [How does el_gato work?](#approach)
* [Using Nextflow](#using-nextflow)

Codebase stage: development   
Developers and maintainers, Testers: [Andrew Conley](https://github.com/abconley), [Lavanya Rishishwar](https://github.com/lavanyarishishwar), [Emily T. Norris](https://github.com/norriset), [Anna Gaines](https://github.com/annagaines), [Will Overholt](https://github.com/waoverholt/), [Dev Mashruwala](https://github.com/dmashruwala), [Alan Collins](https://github.com/Alan-Collins)

# Installation 

## Method 1: using conda
```
# Create an environment named elgato and install el_gato.py plus all dependencies
conda create -n elgato -c bioconda -c conda-forge -c appliedbinf elgato

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
Here is an example of a basic run using paired-end reads or assemblies as input.
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

Legionella in silico sequence-based typing (SBT) script.
    Requires paired-end reads files or a genome assembly.

    Notes on arguments:
    (1) If only reads are provided, SBT is called using a mapping/alignment approach and BLAST.
    (2) If only an assembly is provided, a BLAST and *in silico* PCR-based approach is adopted.

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
  --threads THREADS, -t THREADS
                        Number of threads to run the programs (default: 1)
  --depth DEPTH, -d DEPTH
                        Variant read depth cutoff (default: 3)
  --out OUT, -o OUT     Output folder name (default: out)
  --sample SAMPLE, -n SAMPLE
                        Sample name (default: <Inferred from input file>)
  --overwrite, -w       Overwrite output directory(default: False)
  --sbt SBT, -s SBT     Database containing SBT allele and ST mapping files  
                        (default: .../el_gato/el_gato/db)
  --suffix SUFFIX, -x SUFFIX
                        Suffix of SBT allele files (default: _alleles.tfa)
  --profile PROFILE, -p PROFILE  
                        Name of allele profile to ST mapping file  
                        (default: ../el_gato/el_gato/db/lpneumophila.txt) 
  --verbose, -v         Print what the script is doing (default: False)
  --header, -e          Include column headers in the output table (default: False)
--length LENGTH, -l LENGTH  
                        Specify the BLAST hit length threshold for  
                        identifying multiple loci in an assembly  
                        (default: 0.3)
--sequence SEQUENCE, -s SEQUENCE  
                        Specify the BLAST hit percent identity threshold  for identifying multiple loci in an assembly  
                        (default: 95.0)

  --sequence SEQUENCE, -q SEQUENCE  
                        Specify the BLAST hit percent identity threshold  for identifying multiple loci in assembly (default: 95.0)
```

# Input and Output

## Input files

#### Pair-end reads
When running on a directory of reads, files are associated as pairs using the pattern `R{1,2}.fastq`. i.e., filenames should be identical except for containing either "R1" or "R2" and can be .fastq or .fastq.gz format. Any files for which a pair can not be identified using this pattern will not be processed.

#### Genome assemblies
When running on a directory of assemblies, all files in the target directory will be processed, and there will be no filename restrictions.

## Output files
After a run, el_gato will print the identified ST of your sample to your terminal (stdout) and write several files to the specified output directory (default: out/). 
A file named "all_mlst.txt" (i.e., the ST profile for each tested sample) is written to the output directory and a subdirectory for each sample processed. Each directory is named with its sample name and contains output files specific to that sample (see below). 

### The files included in the output directory for a sample are: 

### standard out
ST profile is written as a tab-delimited table without the headings  
`Sample  ST flaA   pilE asd   mip mompS   proA  neuA_neuAH`   
Headings are included if el_gato.py is run with `-e` flag. The sample column contains the user-provided or inferred sample name. The ST column contains the overall sequence type of the sample. The remaining columns have the allele number of the corresponding gene and the number of reads that support the *mompS* call.

The ST column can contain two kinds of values. If the identified ST corresponds to a profile found in the database, the corresponding number is given. If no matching ST profile is found or el_gato was unable to make a confident call, then this will be reflected in the value displayed in the ST column.

The corresponding allele number is reported for each gene if an exact allele match is found in the database. Alternatively, el_gato may also note the following symbols:

| Symbol | Meaning |
|:------:|:---------|
| NAT    | Novel Allele Type: BLAST cannot find an exact allele match - most likely a new allele. |
| -      | Missing Data: Both percent and length identities are too low to return a match or N's in sequence. |
| ?      | Multiple alleles: More than one allele was found and could not be resolved. |

If symbols are present in the ST profile, the other output files produced by el_gato will provide additional information to understand what is being communicated.

### possible_mlsts.txt
This file would contain all possible ST profiles if el_gato identified multiple possible alleles for any ST loci. In addition, if multiple *mompS* alleles were found, the information used to determine the primary allele is reported in two columns: "mompS_reads_support" and "mompS_reads_against." mompS_reads_support indicates the number of reads associated with each allele that contains the reverse sequencing primer in the expected orientation, which suggests that this is the primary allele. mompS_reads_against indicates the number of reads containing the reverse sequencing primer in the wrong orientation and thus demonstrates that this is the secondary allele. These values are used to infer which allele is the primary *mompS* allele, and their values can be considered to represent the confidence of this characterization [confidence needs a subsection]. [See Approach subsection for more details](#reads).

### intermediate_outputs.txt
el_gato calls other programs to perform intermediate analyses. The outputs of those programs are provided here. In addition, to help with troubleshooting issues, important log messages are also written in this file. The following information may be contained in this file, depending on reads or assembly input:

* Reads-only - Mapping information showing coverage of ST loci by sequencing reads
* Reads-only - Information about the orientation of *mompS* sequencing primer in reads mapping to biallelic sites. [See Approach subsection for more details](#reads).
* Reads-only - Samtools coverage output. [See samtools coverage documentation for more information about headers](https://www.htslib.org/doc/samtools-coverage.html) or [below.](#samtools-coverage-headers)
* BLAST output indicating the best match for identified alleles. [See BLAST output documentation for more information about headers](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/) or [below.](#blastn-output-headers)

Headers are included in outputs for samtools coverage and blast results. Header definitions are as follows:

#### samtools coverage headers

| Column header | Meaning 
|:-------------:|:---------------------------------------------------:|
| rname         | Locus name                                          |
| numreads      | Number reads aligned to the region (after filtering)|
| covbases      | Number of covered bases with depth >= 1             |
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
When run on reads, el_gato maps the provided reads to [a set of reference sequences in the el_gato db directory](https://github.com/appliedbinf/el_gato/blob/main/el_gato/db/ref_gene_regions.fna). The mapped reads are then used to extract the sequences present in the sample to identify the ST. reads_vs_all_ref_filt_sorted.bam and its associated file reads_vs_all_ref_filt_sorted.bai contain the mapping information that was used by el_gato. The BAM file can be viewed using software such as [IGV](https://software.broadinstitute.org/software/igv/) to understand better the data used by el_gato to make allele calls. Additionally, this file is a good starting point for investigating the cause of incorrectly resolved loci.

**Note:** A SAM file is also present, which has the same information as in the BAM file.

# Approach

Most of the processing performed by el_gato is the identification of allele sequences. At its core, el_gato uses BLAST to identify the closest match to each allele in your input data. For the loci *flaA*, *pilE*, *asd*, *mip*, and *proA*, this process is straight forward. Whereas loci *mompS* and *neuA/neuAh* require more involved processing, with neuA/neuAh being an issue only when processing reads—the specifics of these loci are discussed in the corresponding sections below. 


For the simple loci (*flaA*, *pilE*, *asd*, *mip*, and *proA*), the following processes are used:

## Assembly

Six of the seven loci (*flaA*, *pilE*, *asd*, *mip*,*proA*, and *neuA/neuAh*) are identified using BLAST. For each, the best BLAST result is returned as the allele. The closest match is returned with an \* if loci have no exact match. When processing an assembly, only *mompS* requires extra processing. [add in updated code changes]

### *mompS*

[*mompS* is sometimes present in multiple copies in *Legionella pneumophila*, typically two.](https://doi.org/10.1016/j.cmi.2017.01.002) When typing *L. pneumophila* using Sanger sequencing, primers amplify only the correct *mompS* locus. We, therefore, use *in silico* PCR to extract the correct *mompS* locus sequence from the assembly. The primers used for *in silico* PCR are *mompS*-450F (TTGACCATGAGTGGGATTGG) and *mompS*-1116R (TGGATAAATTATCCAGCCGGACTTC) [as described in this protocol](https://doi.org/10.1007/978-1-62703-161-5_6). The *mompS* allele is then identified using BLAST.

## Reads

When processing reads, identification of both *mompS* and *neuA*/*neuAh* requires additional analyses (described below). The other five loci are processed by mapping the provided reads to reference loci from *L. pneumophila* strain Paris and identifying the consensus sequence. Then, all alleles are determined using BLAST against the SBT allele database.

A couple of quality control steps are applied when processing the reads that map to each locus:

   1. **Base quality:** Any bases with quality scores below 20 are not included when calculating coverage at each position or identifying alternate base calls. The lowest number of bases with quality over 20 that map to a single position is reported in the log for each locus [add exact name of the column header].
   2. **Coverage:** After excluding low-quality bases, if there is not at least one read covering 100% of the locus (<99% for *neuA*/*neuAh* - see below), then no attempt to identify the allele is made, and a "-" will be reported. No minimum depth cutoff is applied. [Is depth here true still?]

### *neuA/neuAh*

[The sequence of *neuA*/*neuAh* loci can differ dramatically.](https://doi.org/10.1111/1469-0691.12459) The differences in sequence between *neuA*/*neuAh* alleles are sufficient that reads from some alleles will not map to others. Accordingly, we map reads to [X number of] reference sequences that cover the sequence variation currently represented in the SBT database. The [X number of] reference alleles used are the *neuA* allele from strain Paris (neuA_1), the *neuAh* allele from strain Dallas-1E (neuA_201), and [description of other reference sequences]. Reads from samples are mapped well to the [X number of] reference sequences described above. The reference sequence with the best mapping [what defines best? Is it the number of reads?] is identified using `samtools coverage` with the caveat that > 99% of the locus must have coverage of at least one read (some alleles contain small indels, so 100% is too strict); otherwise a "-" will be reported. Once the reference sequence is selected, the processing using BLAST is the same as described above.

### *mompS*

[As described in the assembly section](#assembly), *mompS* is sometimes present in two copies in the genome of *L. pneumophila* isolates. Duplicate gene copies pose an obvious challenge for a read-mapping approach: if two similar sequence copies are present in a genome, reads from both copies may map to the same reference sequence.

el_gato resolves this issue by taking advantage of the proximity of the genome's two copies of *mompS*. [A schematic of the organization of the two *mompS* copies can be found in Fig. 1 in this paper](https://doi.org/10.1016/j.cmi.2017.01.002). As described in the assembly section, *mompS* is sometimes present in two copies in the genome of *L. pneumophila* isolates. The sequence context of the two mompS copies is such that the correct copy is immediately upstream of the incorrect copy. Only the right copy is flanked on either side by sequences corresponding to primers (*mompS*-450F and *mompS*-1116R) used for conventional SBT. In contrast, while *mompS*-450F sequences are present upstream of the incorrect copy, the corresponding *mompS*-1116R sequences are not found downstream. For this reason, only the correct copy is amplified in conventional SBT (See below schematic)[should be a link to section]. 

Duplicate gene copies pose an obvious challenge for a read mapping approach: reads from both *mompS* copies can map to the same reference sequence. el_gato resolves this issue by taking advantage of the presence of *mompS*-1116R downstream of only the correct copy of *mompS*.

The sequence of the two copies of *mompS* and the identity of the correct allele is then resolved through the following process:

1. Reads from both *mompS* copies are mapped to a single *mompS* reference sequence flanked by the *mompS*-450F and *mompS*-1116R primer sequences. 

2. The nucleotide sequence of reads is recorded for each position within the mompS sequence. If the base at a particular position is heterogeneous in more than 30% of reads mapped to that position, the position is considered biallelic, and both bases are recorded. If the sequence contains only one copy of mompS or no sequence variation between the duplicate copies, then no biallelic sites will be found. Then, the sequence will be extracted, and an allele can be identified using BLAST. 

3. If multiple biallelic positions are identified, all sequences are recorded, and individual read pairs are identified that map to each of the biallelic position(s). 

4. Once the two alleles have been resolved, the correct allele for SBT is identified by analyzing the reads associated with each allele. Reads from each allele are searched for the *mompS*-1116R sequence. The orientation of the reads that contain the primer sequence is assessed. If the primer maps 3'-5' relative to the reference sequence (i.e., in the reverse direction), this is consistent with the read pair originating from the correct copy of *mompS*. However, if the read containing the primer maps 5'-3' (i.e., in the forward direction) relative to *mompS*, this is consistent with the read pair originating from the wrong copy of *mompS*. 

5. The number of reads associated with each allele that contains the primer in the correct orientation relative to *mompS* is counted and compared. The correct allele is then chosen using the following criteria:  
   a. Only one allele has associated reads with primer *mompS*-1116R correctly oriented. 
   b. One allele has more than three times [still valid?] as many reads with correctly oriented primer as the other.   
   c. One allele has no associated reads with the primer *mompS*-1116R in either orientation, but the other allele has associated reads with the primer in only the wrong orientation. In this case, the allele with no associated reads with the primer in either orientation is considered the primary locus by the process of elimination. [still valid?] 
6. The allele number of all alleles is then determined using BLAST and the ST is generated using the correct allele. 

If the above process cannot identify the correct sequence, a ? will be returned as the *mompS* allele, and el_gato will report information about the steps in this process in the [output files](#output-files).

![mompS read mapping schematic](https://github.com/appliedbinf/el_gato/blob/images/images/mompS_allele_assignment.png)

# Using nextflow

[Need to work through these steps for better readme documentation]
We provide a simple nextflow workflow to run el_gato on a directory of either reads or assemblies. In both cases, the target directory must contain only paired reads files (in .fastq or .fastq.gz format) or assembly files (in fasta format).

Uncomment conda environment installation on line 10 and line 47 of the run_el_gato.nf file to run nextflow [is this still necessary?]

```
# Reads
nextflow run_el_gato.nf --reads_dir <path/to/reads/directory> --threads <threads> --out <path/to/output/directory>

# Assemblies
nextflow run_el_gato.nf --assembly_dir <path/to/assemblies/directory> --threads <threads> --out <path/to/output/directory>

```
# Using nextflow with Singularity Container

We provide a singularity container that can be run using the nextflow workflow for el_gato on a directory of either reads or assemblies. In both cases the target directory must contain only paired reads files (in .fastq or .fastq.gz format) or assembly files (in fasta format).

```
# Reads
nextflow run_el_gato.nf --reads_dir <path/to/reads/directory> --threads <threads> --out <path/to/output/directory> -profile singularity -c nextflow.config

# Assemblies
nextflow run_el_gato.nf --assembly_dir <path/to/assemblies/directory> --threads <threads> --out <path/to/output/directory> -profile singularity -c nextflow.config
```