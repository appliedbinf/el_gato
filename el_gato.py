#!/usr/bin/env python3
import argparse
import logging
import os
import re
import subprocess
import sys

class Ref:
    file = "Ref_Paris_mompS_2.fasta"
    name = "Paris_mompS_R"
    source = "Contig = gi|54295983|ref|NC_006368.1| location on contig = 3453389_3455389"
    allele_start = 367
    allele_stop = 718
    flank_start = 15
    flank_stop  = 972
    seq = "GTTATCAATAAAATGGAAACTCAATAATAAACAAGTGGAGACAAGGCATGTTTAGTTTGAAAAAAACAGCAGTGGCAGTACTCGCCTTAGGAAGCGGTGCAGTGTTTGCTGGAACCATGGGACCAGTTTGCACCCCAGGTAATGTAACTGTTCCTTGCGAAAGAACTGCATGGGATATTGGTATCACCGCACTATATTTGCAACCAATCTATGATGCTGATTGGGGCTACAATGGTTTCACCCAAGTTGGTGGCTGGCAGCATTGGCATGATGTTGACCATGAGTGGGATTGGGGCTTCAAATTAGAAGGTTCTTATCACTTCAATACTGGTAATGACATCAATGTGAACTGGTATCATTTTGATAATGACAGTGATCACTGGGCTGATTTTGCTAACTGGCACAACTACAACAACAAGTGGGATGCTGTTAATGCTGAATTAGGTCAATTCGTAGATTTCAGCGCTAACAAGAAAATGCGTTTCCACGGCGGTGTTCAATACGCTCGCATTGAAGCTGATGTGAACCGTTATTTCAATAACTTTGCCTTTAACGGGTTCAACTCTAAGTTCAATGGCTTTGGTCCTCGCACTGGTTTAGACATGAACTATGTATTTGGCAATGGCTTTGGTGTTTATGCTAAAGGCGCTGCTGCTATTCTGGTTGGTACCAGCGATTTCTACGATGGAATCAACTTCATTACTGGTTCTAAAAATGCTATCGTTCCTGAGTTGGAAGCTAAGCTTGGTGCTGATTACACTTACGCAATGGCTCAAGGCGATTTGACTTTAGACGTTGGTTACATGTGGTTTAACTACTTCAACGCTATGCACAATACTGGCGTATTTAATGGATTTGAAACTGATTTCGCAGCTTCTGGTCCTTACATTGGCTTGAAGTATGTTGGTAATGTGTAATTTGTTAAGTTGATAAGAAATTTCAGCAATACTGTTGACTTTATAGAAGTCCGGCTGGATAATTTATCCA"
    analysis_path = ""
    locus_order = ["flaA", "pilE", "asd", "mip", "mompS", "proA", "neuA_neuAH"]
    ispcr_opt = "stdout -out=fa -minPerfect=5 -tileSize=6 -maxSize=1200 -stepSize=5"
    mompS_primer1="TTGACCATGAGTGGGATTGG\tTGGATAAATTATCCAGCCGGACTTC"
    mompS_primer1="TTGACCATGAGTGGGATTGG\tCAGAAGCTGCGAAATCAG"




""" Get commandline arguments """
parser = argparse.ArgumentParser(description="""Legionella in silico SBT script. Needs either reads file and/or genome assembly.
Notes on arguments:
(1) If only reads are provided, de novo assembly is performed and SBT is called using a combination of assembly and mapping. 
(2) If only an assembly is provided, a BLAST and in silico PCR based approach is adopted. 
(3) If both are provided, SBT is called using a combination of assembly and mapping.
""", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-r1", help="Read(s) file", type=str, required=False, metavar="Read 1 file")
parser.add_argument("-r2", help="Read 2 file", type=str, required=False, metavar="Read 2 file")
parser.add_argument("-a", help="Assembly file for isolate", type=str, required=False)
parser.add_argument("-t", help="Number of threads to run the programs", type=int, required=False, default=4)
parser.add_argument("-d", help="Variant read depth cutoff", type=int, required=False, default=3)
parser.add_argument("--prefix", help="Prefix for output files", type=str, required=False, default="run")
parser.add_argument("--out", help="Output folder name", type=str, required=False, default="out")
parser.add_argument("--log", help="Logging file prefix", type=str, required=False, default="run")
parser.add_argument("--overwrite", "-w", help="Overwrite output directory", action="store_true", required=False, default=False)
parser.add_argument("--sbt", "-s", help="Database containing SBT allele files", type=str, required=False, default="./")
parser.add_argument("--suffix", "-x", help="Suffix of SBT allele files", type=str, required=False, default="_alleles.tfa")
parser.add_argument("--profile", help="Allele profile to ST file", type=str, required=False, default="lpneumophila.txt")
args = parser.parse_args()

""" Configuring the logger """
logging.basicConfig(filename=args.log+".log", filemode="w", level=logging.DEBUG, format="[%(asctime)s] %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p")
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter("[%(asctime)s] %(message)s")
console.setFormatter(formatter)
logging.getLogger().addHandler(console)

""" Check for command line arguments """
error = "Not enough arguments! The script requires either both read files and/or a genome assembly file."
Ref.analysis_path=""
if not args.a:
    # assembly file is not supplied, make sure both reads file are supplied
    if not args.r1 or not args.r2:
        logging.critical(f"Error: {error}")
        sys.exit(1)
    Ref.analysis_path = "r"
elif args.r1:
    # assembly is supplied, do we have reads1?
    Ref.analysis_path="a"
    if args.r2:
        # we have reads and assembly
        Ref.analysis_path+="r"
    else:
        # reads2 is missing, which is incorrect
        logging.critical(f"Error: {error}")
        sys.exit(1)
else:
    # assembly is supplied, read1 is not. Is read2 supplied?
    Ref.analysis_path="a"
    if args.r2:
        # only reads2 is supplied, which is incorrect
        logging.critical(f"Error: {error}")
        sys.exit(1)



""" Check for files, folders and dependencies in environment """


#TODO: check if the dependency programs are installed
def check_program(command):
    # Check for bwa
    # Check for sambamba
    # Check for freebayes
    # Check for samtools
    # Check for BLAST (blastn and makeblastdb)
    # Check for isPcr -> packaged with this script
    pass


#TODO: check if the input files exist
def check_files():
    # Check for assembly files if Ref.analysis_path == a
    # Check for reads files if Ref.analysis_path == r
    # Check for assembly and reads files if Ref.analysis_path == ar
    # check if output directory exists, if it does, check if overwrite option is being used.
    # quit otherwise
    # check if the SBT directory exists
    # check if all the alleles within the SBT directory exists
    pass


#TODO: ensure number of threads supplied are not more than system capacity
def ensure_safe_threads():
    pass

""" Functions """
def run_command(command, tool, stdin=None):
    logging.info(f"Running command: {command}")
    result = ""
    if stdin is not None:
        result = subprocess.check_output(command.split(" "), stderr=subprocess.STDOUT, input=bytes(stdin,"utf-8")).decode("utf-8") 
    else: 
        result = subprocess.check_output(command.split(" "), stderr=subprocess.STDOUT).decode("utf-8") 
    logging.debug(f"Command log for {tool}:\n{result}")
    return result

#TODO: add a checkpoint
def add_checkpoint(step):
    pass


#TODO: check if checkpoint exists
def check_checkpoint(step):
    pass


#TODO: resume a run from checkpoint
def resume_checkpoint(step):
    pass


#TODO: Check for reference file and that it is indexed by BWA and FAIDX
def validate_ref():
    if not os.path.exist(Ref.file):
        with open(Ref.file, "w") as f:
            f.write(f">{Ref.name}\n{Ref.seq}\n")
            
    if not os.path.exist(Ref.file+".bwt"):
        run_command(f"bwa index {Ref.file}", "bwa index")
        run_command(f"samtools faidx {Ref.file}", "samtools faidx")\
    

def check_coverage(file):
    logging.info(f"Computing coverage for {file}")
    depth = subprocess.check_output(f"samtools depth -a -r {Ref.name}:{Ref.allele_start}-{Ref.allele_stop} {file}".split(" ")).decode("utf-8").split("\n")
    for d in depth:
        d = d.rstrip().split("\t")
        if len(d) < 3:
            break
        if int(d[2]) < args.d:
            logging.warn(f"Low depth base (depth={d[2]}) found in {file}, at pos {d[0]}:{d[1]}.")
            return False
    
    logging.info(f"File {file} passes depth check.")
    return True


def call_variants(prefix):
    # SAM -> BAM and sort
    sam2bam = f"sambamba view -f bam -S -t {args.t} -o {prefix}.bam {prefix}.sam"
    run_command(sam2bam, "sambamba SAM to BAM conversion")
        
    sort_bam = f"sambamba sort -t {args.t} {prefix}.bam"
    run_command(sort_bam, "sambamba sort BAM")
    
    # Call variants
    freebayes_call = f"freebayes -v {prefix}.vcf -f {Ref.file} {prefix}.sorted.bam"
    run_command(freebayes_call, "freebayes")
    
    # Check that there is coverage across the entire gene region
    return check_coverage(file=f"{prefix}.sorted.bam")
    

def filter_variant_calls(samfile, outfile):
    proper_pairs = f"samtools view -h -f 0x2 {samfile}"
    proper_pairs = run_command(proper_pairs, "samtools view").rstrip().split("\n")
    
    reads_of_interest = {}
    out_text = ""
    all_reads = {}
    for line in proper_pairs:
        if line.startswith("@"): 
            out_text += line + "\n"
            continue
        cols = line.rstrip().split("\t")
        read_id  = cols[0]
        read_start = int(cols[3])
        
        if read_id in all_reads:
            all_reads[read_id] += line
        else :
            all_reads[read_id] = line
        
        region_start = Ref.flank_start
        #TODO: Improve region_end calculation
        # region_end is a way to check that one of the pair is in the right flank.
        # a better way of checking this will be to look at the right-most mapping position
        region_end   = Ref.flank_stop - int(re.search("\d+M", cols[5]).group()[:-1])
        if (read_start < region_start or read_start > region_end):
            reads_of_interest[read_id] = 1
    
    with open(outfile, "w") as fh:
        fh.write(out_text)
        for row in all_reads:
            if row in reads_of_interest:
                fh.write(all_reads[row] + "\n")


def vcf_to_fasta(full_vcf, filtered_vcf):
    this_seq = ""
    start_anchor = 0
    
    filtered_call = {}
    with open(filtered_vcf, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            (chr, pos, id, ref, alt, qual, filter, info, format, gt) = line.rstrip().split("\t")
            pos = int(pos)
            if float(re.search("AB=[0-9.]+;", info).group()[3:-1]) == 0:
                # homozygous call
                filtered_call[pos] = alt
    
    with open(full_vcf, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            (chr, pos, id, ref, alt, qual, filter, info, format, gt) = line.rstrip().split("\t")
            pos = int(pos)
            # check if this position is within the typing allele
            if Ref.allele_start <= pos and pos <= Ref.allele_stop:
                # check for zygosity
                ab = re.search("AB=[0-9.,]+;", info).group()[3:-1]
                if re.search(",", ab) or float(ab) != 0:
                    # heterozygous call
                    # Two scenarios:
                    # 1. POS exist in filtered file          => add the alternative base from filtered_vcf
                    # 2. POS doesn't exist in filtered file  => add the reference base
                    if pos in filtered_call:
                        this_seq += Ref.seq[start_anchor:(pos-1)] + filtered_call[pos]
                        start_anchor = pos-1+len(ref)
                    else :
                        pos -= 1 # 0-based in list
                        this_seq += Ref.seq[start_anchor:pos] + ref
                        start_anchor = pos+len(ref)
                else :
                    # homozygous call
                    pos -= 1 # 0-based in list
                    this_seq += Ref.seq[start_anchor:pos] + alt
                    start_anchor = pos+len(ref)
                
    this_seq += Ref.seq[start_anchor:]
    return this_seq


def blast_mompS_allele(seq, db=os.path.join(args.sbt, "mompS"+args.suffix)):
    makeblastdb = f"makeblastdb -in {db} -dbtype nucl"
    run_command(makeblastdb, "makeblastdb/mompS")
    blastcmd = f"blastn -query - -db {db} -outfmt 6 -perc_identity 100"
    res = run_command(blastcmd, "blastn/mompS", seq).rstrip()
    if res == "":
        # TODO: run blast again with lower identity threshold and return allele*
        return "-"
    else :
        cols = res.split("\t")
        allele = cols[1].replace("mompS_", "")
        return allele
    

def call_mompS_mapping(r1, r2, outfile, filt_file = ""):
    
    if filt_file == "":
        filt_file = outfile+".filtered"
    
    # Map reads to mompS gene
    bwa_call = f"bwa mem -t {args.t} {Ref.file} {r1} {r2} -o {outfile}.sam"
    run_command(bwa_call, "bwa")
    
    # Create a separate file containing reads coming from the border regions
    filter_variant_calls(samfile=f"{outfile}.sam", outfile=f"{filt_file}.sam")
    
    # Call variants
    full_coverage = call_variants(prefix=outfile)
    filtered_coverage = call_variants(prefix=filt_file)
    
    allele_confidence = ""
    if not (full_coverage or filtered_coverage):
        allele_confidence = "?"
    
    allele_seq = vcf_to_fasta(full_vcf=outfile+".vcf", filtered_vcf=filt_file+".vcf")
    allele_seq = allele_seq[(Ref.allele_start-1):Ref.allele_stop]
    
    mompS_allele = blast_mompS_allele(allele_seq)
    if mompS_allele != "-":
        mompS_allele += allele_confidence
        
    return(mompS_allele)
    

def call_mompS_pcr(assembly_file = args.a, db=os.path.join(args.sbt, "mompS"+args.suffix)):
    makeblastdb = f"makeblastdb -in {db} -dbtype nucl"
    run_command(makeblastdb, "makeblastdb/mompS")
    blast_command = f"blastn -db {db} -outfmt 6 -query {assembly_file} -perc_identity 100"
    res = run_command(blast_command, "blastn/mompS").rstrip().split("\n")
    res = [ line.rstrip().split("\t")[1] for line in res]
    
    alleles = {}
    for allele in res:
        alleles[allele] = 1
    alleles = list(alleles.keys())
    
    if len(alleles) == 1:
        return alleles[0]
    else :
        primer1 = os.path.join(args.out, "mompS_primer1.tab")
        with open(primer1, "w") as f:
            f.write("mompS_1\t"+Ref.mompS_primer1)
        ispcr_command = f"isPcr {assembly_file} {primer1} {Ref.ispcr_opt}"
        primer1_res = run_command(ispcr_command, "mompS2 primer1")
        
        if len(primer1_res.rstrip().split("\n")) > 0 :
            # nested PCR
            primer2 = os.path.join(args.out, "mompS_primer2.tab")
            with open(primer2, "w") as f:
                f.write("mompS_2\t"+Ref.mompS_primer2)
            ispcr_command = f"isPcr {assembly_file} {primer2} {Ref.ispcr_opt}"
            primer2_res = run_command(ispcr_command, "mompS2 primer2", primer1_res).rstrip().split("\n")[1]
            logging.debug(f"Found the sequence: {primer2_res}")
            return blast_mompS_allele(primer2_res)
        else:
            # mompS can't be typed with assembly
            #TODO: recommend that the user supply reads file
            return "-"
    

def genome_assembly(r1, r2, out=os.path.join(args.out, "run_spades")):
    assembly_command = f"spades.py -1 {r1} -2 {r2} -o {out} --careful -t {args.t}"
    run_command(assembly_command)
    args.a = os.path.join(out, "scaffolds.fasta")


def blast_non_mompS(assembly_file = args.a):
    loci = ["flaA", "pilE", "asd", "mip", "proA", "neuA_neuAH"]
    calls = {}
    for locus in loci:
        db = os.path.join(args.sbt, locus+args.suffix)
        makeblastdb = f"makeblastdb -in {db} -dbtype nucl"
        run_command(makeblastdb, f"makeblastdb/{locus}")
        blastcmd = f"blastn -query {assembly_file} -db {db} -outfmt 6 -perc_identity 100"
        res = run_command(blastcmd, f"blastn/{locus}").rstrip()
        allele = ""
        if res == "":
            # TODO: run blast again with lower identity threshold and return allele*
            allele = "-"
        else :
            cols = res.split("\t")
            allele = cols[1].replace(locus+"_", "")
        calls[locus] = allele
    return calls


def get_st(allele_profile, profile_file = args.profile):
    with open(profile_file, "r") as f:
        f.readline()
        for line in f:
            line = line.rstrip()
            if line.endswith(allele_profile):
                st = line.split("\t")[0]
                return st
    # following system call craps out, will debug it someday
    # allele_profile = allele_profile.replace("\t", "\\t")
    # grep_command = f"grep -P \'{allele_profile}$\' {profile_file}"
    # st = run_command(grep_command, "Retreiving ST").rstrip().split("\t")[0]
    # return st


#TODO: make decision on the path to choose
# if only reads are provided, run genome_assembly and go to call_mompS_mapping
# if only assembly is provided, go to call_mompS_pcr
# if both are provided, go to call_mompS_mapping
def choose_path():
    alleles = {}
    if Ref.analysis_path == "ar":
        mompS_allele = call_mompS_mapping(r1=args.r1, r2=args.r2, outfile=os.path.join(args.out, args.prefix))
        alleles = blast_non_mompS()
        alleles["mompS"] = mompS_allele
    elif Ref.analysis_path == "a":
    elif Ref.analysis_path == "a":
    else:
        logging.critical("This path should not have been traversed. Is Ref.analysis_path being changed somewhere else?")


""" Main code """


allele_profile = ""
for locus in Ref.locus_order:
    allele_profile += alleles[locus] + "\t"
allele_profile = allele_profile.rstrip()
allele_profile = get_st(allele_profile) + "\t" + allele_profile

print("ST\t"+"\t".join(Ref.locus_order))
print(allele_profile)
