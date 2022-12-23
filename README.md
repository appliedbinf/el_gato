# el_gato
Epidemiology of *Legionella*: Genome-bAsed Typing

* [Installation](#installation)
   * [Method 1: using conda ](#method-1-using-conda)
   * [Method 2: using pip](#method-2-using-pip)
* [Quickstart guide](#Quickstart-guide)
   * [A typical run](#A-typical-run)   
* [Output](#Output)
     * [stdout (MLST profile)](#stdout-MLST-profile)
     * [possible_mlsts.txt](#possible_mlststxt)
     * [intermediate_outputs.txt](#intermediate_outputstxt)
     * [identified_alleles.fna](#identified_allelesfna)
     * [run.log](#runlog)
     * [reads_vs_all_ref_filt_sorted.bam](#reads_vs_all_ref_filt_sortedbam-reads-only)
* [Using Nextflow](#Using-nextflow)
   * [Input files](#Input-files)
   * [Output files](#Output-files)
* [All available arguments](#Usage)
* [How does el_gato work?](#Approach)
* [Dependencies](#Dependencies)

Currently in development
Codebase stage: development   
Developers and maintainers, Testers: [Andrew Conley](https://github.com/abconley), [Lavanya Rishishwar](https://github.com/lavanyarishishwar), [Emily T. Norris](https://github.com/norriset), [Anna Gaines](https://github.com/annagaines), [Will Overholt](https://github.com/waoverholt/), [Dev Mashruwala](https://github.com/dmashruwala), [Alan Collins](https://github.com/Alan-Collins)

# Installation 

## Method 1: using conda

```
# Create environment named elgato and install el_gato.py plus all dependencies
conda create -n elgato -c bioconda -c conda-forge -c appliedbinf elgato

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

```

# Output

Upon the completion of a run, el_gato.py will print the identified MLST of your sample to your terminal (stdout) and will write a number of files to the spcified output directory.

## stdout (MLST profile)

MLST profile is written as a tab-delimited table with the headings `Sample  ST flaA   pilE asd   mip mompS   proA  neuA_neuAH`. The sample column contains the user-provided or inferred sample name. The ST column contains the overall sequence type of the sample. The remaining columns contain the allele number of the corresponding gene.

The ST column can contain two kinds of values. If the identified MLST corresponds to a profile found in the database, the corresponding number is given. If no matching MLST profile is found, "NF" is reported.

For each gene, if an exact allele match is found in the database, the corresponding allele number is reported. Alternatively, the following symbols may also be reported:

| Symbol | Meaning |
|:------:|:---------|
| *      | No exact allele match was found. The closest match allele is reported with an asterisk. |
| -      | Some or all of the locus is missing. No match can be found. |
| ?      | More than one allele was found. |

In the case that any of these symbols are present in the MLST profile, the other output files produced by el_gato will provide more information to understand what is being communicated.

## possible_mlsts.txt

In the case that multiple alleles were identified for any MLST loci, this file will contain all possible ST profiles. In addition, if multiple mompS alleles were found, the information that was used to try to identify the primary allele is reported ([see Approach section for more details](#Reads)).

## intermediate_outputs.txt

El_gato calls other programs to perform intermediate analyses. The outputs of those programs is provided here. In addition, to help with troubleshooting issues important log messages are also written to this file. The following information may be contained in this file, depending on reads or assembly input:

* (Reads-only) Mapping information showing coverage of MLST loci by sequencing reads
* (Reads-only) Information about the orientation of mompS sequencing primer in reads mapping to bialleleic sites ([see Approach section for more details](#Reads)).
* BLAST output indicating the best match for identified alleles.
* Important logging information.

## identified_alleles.fna

The sequence of all identified alleles are written to this file. If more than one allele is identified for the same locus, they are numbered in an arbitrary order. Fasta headers of sequences in this file correspond to the query IDs in the BLAST output reported in the intermediate_outputs.txt file.

## run.log

Detailed log of the steps taken during the running of el_gato including the outputs of any programs called by el_gato and any errors encountered.

## reads_vs_all_ref_filt_sorted.bam (reads only)

When run on reads, el_gato maps the provided reads to [a set of reference sequences in the el_gato db directory](https://github.com/appliedbinf/el_gato/blob/main/el_gato/db/ref_gene_regions.fna). The mapped reads are then used to extract the sequences present in the sample to identify the MLST. reads_vs_all_ref_filt_sorted.bam and its associated file reads_vs_all_ref_filt_sorted.bai contain the mapping information that was used by el_gato. The BAM file can be viewed using software such as [IGV](https://software.broadinstitute.org/software/igv/) to get a better understanding of the information used by el_gato to make allele calls. Additionally, if any loci were not properly resolved, this file is a good starting point for figuring out why.

N.B., a SAM file is also present. This is the same information as in the BAM file.

# Using nextflow

We provide a simple nextflow workflow to run el_gato on a directory of either reads or assemblies. In both cases the target directory must contain only paired reads files (in .fastq or .fastq.gz format) or assembly files (in fasta format).

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
nextflow run_el_gato.nf --reads_dir <path/to/reads/directory> --threads <threads> --out <path/to/output/directory> -profile sge,singularity -c nextflow.config

# Assemblies
nextflow run_el_gato.nf --assembly_dir <path/to/assemblies/directory> --threads <threads> --out <path/to/output/directory> -profile sge,singularity -c nextflow.config

```
## Input files

When running on a directory of reads, files are associated as pairs using the pattern `*R{1,2}*.fastq*`. I.e., filenames should be identical except for containing either "R1" or "R2" and can be .fastq or .fastq.gz format. Any files for which a pair can not be identified using this pattern will not be processed.

When running on a directory of assemblies, all files in the target directory will be processed and there are no filename restrictions.

While el gato performs Q20 processing for reads, it is highly recommended to perform preferred QA/QC on input files. 

````
Ex: fastp -i <input_R1.fastq.gz> -I <input_R2.fastq.gz> -o <trimmed_R1.fastq.gz> -O <trimmed_R2.fastq.gz>
````

## Output files

At the completion of a run, the specified output directory (default: el_gato_out/) will contain a file named "all_mlst.txt" (the MLST profile of each sample) and one directory for each sample processed. Each directory is named with a sample name and contains output files specific to tht sample. These files include the el_gato log file and files providing more details about the sequences identified in the sample. [See the Output section](#Output) for more details.


# Usage

```
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
```

# Approach

At its core, el_gato uses BLASTn to identify the closest match to each alleles present in your input sequence. Most of the processing performed by el_gato is the identification of allele sequences. For the loci *flaA*, *pilE*, *asd*, *mip*, and *proA*, this process is straight forward. However, there are two loci for which more involved processing is required: *mompS* and *neuA*/*neuAh* (N.B., *neuA*/*neuAh* processing is only more complex when processing reads.) The specifics of how those two loci are processed is discussed below.

For the simple loci, to following processes are used:

## Assembly

When processing an assembly, only *mompS* requires extra processing. The other 6 loci are identified using BLASTn. For each, the best BLAST result is returned as the allele. If any loci have no exact match, then the closest match is returned with an \*.

### mompS

[*mompS* is sometimes present in two copies in *Legionella pneumophila*](https://doi.org/10.1016/j.cmi.2017.01.002) and therefore requires additional processing. When typing *L. pneumophila* using Sanger sequencing, primers are used that amplify only the correct *mompS* locus. We therefore use *in silico* PCR to extract the correct *mompS* locus sequence from the assembly. The pimers used for *in silico* PCR are *mompS*-450F (TTGACCATGAGTGGGATTGG) and *mompS*-1116R (TGGATAAATTATCCAGCCGGACTTC) [as described in this protocol](https://doi.org/10.1007/978-1-62703-161-5_6). The *mompS* allele is then identified using BLASTn.

## Reads

When processing reads, both *mompS* and *neuA*/*neuAh* must be processed separately. The other 5 loci are processed by mapping the provided reads to reference loci from *L. pneumophila* strain Paris and identifying the consensus sequence. Alleles are then identified using BLASTn.

### neuA/neuAh

[The sequence of *neuA*/*neuAh* loci can differ dramatically.](https://doi.org/10.1111/1469-0691.12459) The differences in sequence between *neuA*/*neuAh* alleles is sufficient that reads from some alleles will not map to others. Accordingly, we map reads to three reference sequences that cover the sequence variation currently represented in the SBT. The three reference alleles used are the *neuA* allele from strain Paris (neuA_1), the *neuAh* allele from strain Dallas-1E (neuA_201) and a chimeric sequence composed of sequence 5' of the *neuA* locus from strain Paris, the sequence of allele neuA_206, and sequence 3' of the *neuA* locus from strain Dallas-1E.

Using the three reference sequences described above, reads from all samples we have tested map well to only one reference. The reference sequence with the best mapping is identified using `samtools coverage`; only those reference loci with coverage over 99% (some alleles contain small indels so 100% is too strict) are retained for future processing. Once the best reference sequence (or sequences if two map with sufficient coverage) is identified, processing is the same as described above.

### mompS

[As described in the assembly section](#Assembly), *mompS* is sometimes present in two copies in the genome of *L. pneumophila* isolates. This poses an obvious challenge for a read-mapping approach: if two similar copies of a sequence are present in a genome, both copies may map to the same reference sequence.

The approach taken by el_gato to resolve this issue takes advantage of the two copies of *mompS* being close together in the genome. [A schematic of the organization of the two *mompS* copies can be found in Fig. 1 in this paper](https://doi.org/10.1016/j.cmi.2017.01.002). The sequence context of the two *mompS* copies is such that only one copy has **both** of the primers *mompS*-450F and *mompS*-1116R flanking it in the correct orientation for PCR amplification, while the other copy is only flanked by primer *mompS*-450F in the correct orientation for amplification (See below schematic).

The sequence of the two copies of *mompS* and the identity of the correct allele is then resolved through the following process:
1. Reads from the two *mompS* copies are mapped to a reference sequence containing only one copy of *mompS* and with the *mompS*-450F and *mompS*-1116R primers present.
2. The sequence of reads mapped to each position within the *mompS* sequence is recorded. If more than 30% of reads mapped to a position contain a different base, then the position is considered to be bialleleic and both bases are recorded. If there was in fact only one copy of *mompS* in the sample, then no bialleleic sites will be found - the sequence has now been extracted and an allele can be identified using BLASTn.
3. If multiple bialleleic positions are identified, then the two sequences are resolved by identifying individual read pairs that map to all bialleleic positions.
4. Once the two sequences have been resolved, the correct sequence is identified is identified by analyzing the reads associated with each sequence (i.e., that map to one or more bialleleic position and contain a base from one of the sequences at that position). The read from each list is searched for the sequence of primer *mompS*-1116R and the orientation of reads that contain the primer is assessed. If a read contains the primer and is mapped 3'-5' relative to the reference sequence (i.e., in the reverse direction), this is consistent with the read pair originating from the correct copy of *mompS*. However, if the read containing the primer is mapped 5'-3' (i.e., in the forward direction), this is consistent with the read pair originating from the wrong copy of *mompS*.
5. The number of reads associated with each sequence that contain the primer *mompS*-1116R in the correct orientation is counted and compared. The correct allele is then decided using the following criteria:
   a. Only one sequence has associated reads with correctly oriented primer.
   b. One sequence has more than three times as many reads with correctly oriented primer as the other.
   c. One sequence has no associated reads with the primer in either orientation, but the other has associated reads with the primer in only the wrong orientation.
6. The allele of both identified sequences is then identified using BLASTn.

If the above process is unable to identify the correct sequence, a ? will be returned as the *mompS* allele and information about the steps in this process will be reported in the [possible_mlsts.txt](#possible_mlsts.txt), [intermediate_outputs.txt](#intermediate_outputs.txt), [identified_alleles.fna](#identified_alleles.fna), and [run.log](#run.log) files.

![mompS read mapping schematic](https://github.com/appliedbinf/el_gato/blob/images/images/mompS_allele_assignment.png)

# Dependencies

* [minimap2](https://github.com/lh3/minimap2)
* [SAMTools](https://github.com/samtools/samtools)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [isPcr](https://users.soe.ucsc.edu/~kent/)
