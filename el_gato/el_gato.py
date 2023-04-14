#!/usr/bin/env python3
import argparse
import inspect
import logging
import multiprocessing
import os
import re
import subprocess
import sys
import shutil
import shlex
import time
from collections import defaultdict, Counter
t0 = time.time()
script_filename = inspect.getframeinfo(inspect.currentframe()).filename
script_path = os.path.dirname(os.path.abspath(script_filename))

class Ref:
    file = "Ref_Paris_mompS_2.fasta"
    name = "Paris_mompS_R"
    source = "Contig = gi|54295983|ref|NC_006368.1| location on contig = 3453389_3455389"
    seq = "GTTATCAATAAAATGGAAACTCAATAATAAACAAGTGGAGACAAGGCATGTTTAGTTTGAAAAAAACAGCAGTGGCAGTACTCGCCTTAGGAAGCGGTGCAGTGTTTGCTGGAACCATGGGACCAGTTTGCACCCCAGGTAATGTAACTGTTCCTTGCGAAAGAACTGCATGGGATATTGGTATCACCGCACTATATTTGCAACCAATCTATGATGCTGATTGGGGCTACAATGGTTTCACCCAAGTTGGTGGCTGGCAGCATTGGCATGATGTTGACCATGAGTGGGATTGGGGCTTCAAATTAGAAGGTTCTTATCACTTCAATACTGGTAATGACATCAATGTGAACTGGTATCATTTTGATAATGACAGTGATCACTGGGCTGATTTTGCTAACTGGCACAACTACAACAACAAGTGGGATGCTGTTAATGCTGAATTAGGTCAATTCGTAGATTTCAGCGCTAACAAGAAAATGCGTTTCCACGGCGGTGTTCAATACGCTCGCATTGAAGCTGATGTGAACCGTTATTTCAATAACTTTGCCTTTAACGGGTTCAACTCTAAGTTCAATGGCTTTGGTCCTCGCACTGGTTTAGACATGAACTATGTATTTGGCAATGGCTTTGGTGTTTATGCTAAAGGCGCTGCTGCTATTCTGGTTGGTACCAGCGATTTCTACGATGGAATCAACTTCATTACTGGTTCTAAAAATGCTATCGTTCCTGAGTTGGAAGCTAAGCTTGGTGCTGATTACACTTACGCAATGGCTCAAGGCGATTTGACTTTAGACGTTGGTTACATGTGGTTTAACTACTTCAACGCTATGCACAATACTGGCGTATTTAATGGATTTGAAACTGATTTCGCAGCTTCTGGTCCTTACATTGGCTTGAAGTATGTTGGTAATGTGTAATTTGTTAAGTTGATAAGAAATTTCAGCAATACTGTTGACTTTATAGAAGTCCGGCTGGATAATTTATCCA"
    allele_start = 367
    allele_stop = 718
    flank_start = 15
    flank_stop = 972
    analysis_path = ""
    locus_order = ["flaA", "pilE", "asd", "mip", "mompS", "proA", "neuA_neuAH"]
    ispcr_opt = "stdout -out=fa -minPerfect=5 -tileSize=6 -maxSize=1200 -stepSize=5"
    mompS_primer1 = "TTGACCATGAGTGGGATTGG\tTGGATAAATTATCCAGCCGGACTTC"
    mompS_primer2 = "TTGACCATGAGTGGGATTGG\tCAGAAGCTGCGAAATCAG"
    prereq_programs = ["minimap2", "samtools", "makeblastdb", "blastn", "isPcr"]
    REF_POSITIONS = {
        "asd": {        
            'start_pos' : 351,
            'end_pos' : 823,
        },
        "flaA": {
            'start_pos' : 351,
            'end_pos' : 532,
        },
        "mip": {
            'start_pos' : 350,
            'end_pos' : 751,
        },
        "neuA": {
            'start_pos' : 350,
            'end_pos' : 703,
        },
        "neuAh": {
            'start_pos' : 350,
            'end_pos' : 703,
        },
        "neuA_207": {
            'start_pos' : 350,
            'end_pos' : 700,
        },
        "neuA_211": {
            'start_pos' : 350,
            'end_pos' : 700,
        },
        "pilE": {
            'start_pos' : 351,
            'end_pos' : 683,
        },
        "proA": {
            'start_pos' : 350,
            'end_pos' : 754,
        },
        "mompS": {
            'start_pos' : 367,
            'end_pos' : 718,
        },
    }

class SAM_data(object):
    """stores columns of SAM entry as attributes"""
    def __init__(self, object):
        self.qname = object.split('\t')[0]
        self.flag = int(object.split('\t')[1])
        self.rname = object.split('\t')[2]
        self.pos = int(object.split('\t')[3])
        self.mapq = int(object.split('\t')[4])
        self.cigar = object.split('\t')[5]
        self.rnext = object.split('\t')[6]
        self.pnext = object.split('\t')[7]
        self.seq = object.split('\t')[9]
        self.qual = object.split('\t')[10]
        self.ln = len(self.seq)
        self.end = self.pos + self.ln
        self.mod_seq = ''
        self.mod_qual = ''

        self.cigar_mod_read()
        self.refresh()


    def refresh(self):
        self.seq = self.mod_seq
        self.qual = self.mod_qual
        self.ln = len(self.seq)
        self.end = self.pos + self.ln

    def cigar_mod_read(self):
        """
        According to the M, I, and D components of a cigar string, modifies a seq and quality score strings so that they are in register with refernce sequence 
        returns: modified sequence, modified quality score string
        """
        cig_sects = re.findall('[0-9]+[A-Z]', self.cigar)

        count = 0

        for sect in cig_sects:
            letter = sect[-1]
            sect_len = int(sect[:-1])
            
            if letter == 'M':
                self.mod_seq += self.seq[count:count+sect_len] # Add corresponding portion of original seq
                self.mod_qual += self.qual[count:count+sect_len]
                count += sect_len

            elif letter == 'D':
                self.mod_seq += '*' * sect_len # Add asterisks for each deleted position relative to the reference sequence
                self.mod_qual += ' ' * sect_len


            elif letter == 'I' or letter == 'S':
                count += sect_len

    def get_base_calls(self, positions):
        """
        Args:
            positions: list
                list of positions of base calls desired

        returns:
            list
                base calls
        """

        bases = {k:'' for k in range(len(positions))}
 
        for n, p in enumerate(positions):
            p_mod = p+1 - self.pos
            if p_mod < 0:
                continue
            if p_mod < self.ln:
                bases[n] = self.seq[p_mod]
            else: 
                continue

        return bases

class Allele():
    def __init__(self, num: int = 0, basecalls: str = ''):
        self.num_locs = num
        self.basecalls = basecalls # base calls at multiallelic sites
        self.reads_at_locs = [[]]*num # [[reads, at, loc, 1], [reads, at, loc, 2], ...]
        self.confidence = {"for": 0, "against": 0} # How many reads contain primer sequence? i.e., evidence for native locus
        self.seq = ''
        self.location = ''
        self.fasta_header = ''
        self.allele_id = '-'

    def assess_conf(self):
        if self.confidence['for'] == 'NA':
            pass
        elif self.confidence['for'] > 0:
            self.location = f"_native_locus_{self.confidence['for']}_reads"
        elif self.confidence['against'] > 0:
            self.location = f"_non-native_locus"

def get_args() -> argparse.ArgumentParser:

    """ Get commandline arguments """
    parser = argparse.ArgumentParser(description="""Legionella in silico SBT script. 
    Requires paired-end reads files or a genome assembly.

    Notes on arguments:
    (1) If only reads are provided, SBT is called using a mapping/alignment approach.
    (2) If only an assembly is provided, a BLAST and in silico PCR based approach is adopted. 
    """, formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False)
    group1 = parser.add_argument_group(title='Input files',
                                       description="Please specify either reads files and/or a genome assembly file")
    group1.add_argument("--read1", "-1", help="Input Read 1 (forward) file", type=str, required=False,
                        metavar="Read 1 file")
    group1.add_argument("--read2", "-2", help="Input Read 2 (reverse) file", type=str, required=False,
                        metavar="Read 2 file")
    group1.add_argument("--assembly", "-a", help="Input assembly fasta file", type=str, required=False,
                        metavar="Assembly file")
    group2 = parser.add_argument_group(title='Optional arguments')
    group2.add_argument("--help", "-h", action="help", help="Show this help message and exit")
    group2.add_argument("--threads", "-t", help="Number of threads to run the programs (default: %(default)s)", type=int,
                        required=False, default=1)
    group2.add_argument("--depth", "-d", help="Variant read depth cutoff (default: %(default)s)", type=int, required=False,
                        default=3)
    group2.add_argument("--out", "-o", help="Output folder name (default: %(default)s)", type=str, required=False,
                        default="out")
    group2.add_argument("--sample", "-n", help="Sample name (default: %(default)s)", type=str, required=False,
                        default="<Inferred from input file>")
    group2.add_argument("--overwrite", "-w", help="Overwrite output directory (default: %(default)s)", action="store_true",
                        required=False, default=False)
    group2.add_argument("--sbt", "-s", help="Database containing SBT allele and ST mapping files (default: %(default)s)",
                        type=str, required=False, default=os.path.join(os.path.dirname(__file__), "db"))
    group2.add_argument("--suffix", "-x", help="Suffix of SBT allele files (default: %(default)s)", type=str,
                        required=False, default="_alleles.tfa")
    group2.add_argument("--profile", "-p", help="Name of allele profile to ST mapping file (default: %(default)s)",
                        type=str, required=False, default=os.path.join(os.path.dirname(__file__), "db", "lpneumophila.txt"))
    group2.add_argument("--verbose", "-v", help="Print what the script is doing (default: %(default)s)",
                        action="store_true", required=False, default=False)
    group2.add_argument("--header", "-e", help="Include column headers in the output table (default: %(default)s)", action="store_true", required=False, default=False),


    return parser


def fasta_to_dict(FASTA_file: str) -> dict: 
    """Read a fasta file into a dict    
    Dict has headers (minus the > symbol) as keys and the associated    
    sequence as values. 
        
    Args:   
      FASTA_file: str   
        path to fasta format file   
    Returns:    
      dict:     
        dict of format {fasta_header : sequence}    
    Raises: 
      TypeError: If FASTA_file is not a str 
      OSError: If FASTA_file is not the path to an existing file    
    """ 

    if type(FASTA_file) is not str: 
        raise TypeError(    
            "FASTA_file must be str, not {}.".format(type(FASTA_file).__name__))    

    if not os.path.exists(FASTA_file):  
        raise OSError(  
            "FASTA_file must be the path to an existing file.") 


    fasta_dict = {} 
    with open(FASTA_file, 'r') as f:    
        multifasta = f.read()   
    f.close()   
    fastas = multifasta.split(">")  
    trimmed_fastas = [] 
    for i in fastas:    
        if len(i) != 0: 
            trimmed_fastas.append(i)    

    fastas = trimmed_fastas 

    for i in fastas:    
        header = i.split()[0]   
        seq = "".join(i.split("\n")[1:])    
        fasta_dict[header] = seq    

    return fasta_dict


def rev_comp(string) -> str:
    """Reverse complement a string of nucleotide sequence
    
    Args:
      string: str
          Nucleotide sequence

    Returns:
      str:
        reverse complement of given nucleotide sequence
    
    Raises:
      TypeError: If string is not str.

    """
    if type(string) is not str:
        raise TypeError(
            "string must be str, not {}.".format(type(string).__name__))

    rev_str = ''
    rev_comp_lookup = {
    "A" : "T", 
    "T" : "A", 
    "C" : "G", 
    "G" : "C", 
    "a" : "t", 
    "t" : "a", 
    "c" : "g", 
    "g" : "c",
    }
    for i in reversed(string):
        if i in "ATCGatcg":
            rev_str += rev_comp_lookup[i]
        else:
            rev_str += i
    return rev_str


def check_input_supplied(
        args: argparse.ArgumentParser,
        parser: argparse.ArgumentParser,
        inputs: dict
        ) -> dict:
    """Checks if all the required arguments have been supplied

    Parameters
    ----------
    args: argparse.ArgumentParser
        The arguments with which program was executed, parsed using 
        argparse.

    parser: argparse.ArgumentParser
        Argparse command line parser object

    inputs: dict
        Run settings

    Returns
    -------
    dict
        modified inputs dict
        
        Exits the program if required arguments are not supplied or if they are mismatched
    """
    error = f"""Wrong arguments! The script requires either a genome assembly file or two paired-end read files.\n\nRun {parser.prog} -h to see usage."""
    if not args.assembly:
        # assembly file is not supplied, make sure both reads file are supplied
        if not args.read1 or not args.read2:
            print(f"Error: {error}")
            sys.exit(1)
        inputs["analysis_path"] = "r"
        inputs["logging_buffer_message"] += "User supplied both reads files, adopting the assembly/mapping/alignment route\n"
    else:
        # assembly is supplied, read1 is not. Is read2 supplied?
        inputs["analysis_path"] = "a"
        if args.read1 or args.read2:
            # assembly and reads file(s) is supplied, which is incorrect
            print(f"Error: {error}")
            sys.exit(1)
        inputs["logging_buffer_message"] += "User supplied the assembly file, adopting the alignment/in silico pcr route\n"

    return inputs


def set_inputs(
        args: argparse.ArgumentParser,
        inputs: dict
        ) -> dict:
    """Carries over the arguments supplied in argparse into a local class Inputs

    Parameters
    ----------
    args: argparse.ArgumentParser
        The arguments with which program was executed, parsed using 
        argparse.

    inputs: dict
        Run settings

    Returns
    -------
    dict
        modified inputs dict
    """
    if "r" in inputs["analysis_path"]:
        inputs["read1"] = args.read1
        inputs["read2"] = args.read2
    if "a" in inputs["analysis_path"]:
        inputs["assembly"] = args.assembly
    inputs["threads"] = args.threads
    inputs["out_prefix"] = args.out
    inputs["log"] = os.path.join(args.out, "run.log")
    inputs["sbt"] = args.sbt
    inputs["suffix"] = args.suffix
    inputs["profile"] = args.profile
    inputs["verbose"] = args.verbose
    inputs["overwrite"] = args.overwrite
    inputs["depth"] = args.depth
    inputs["header"] = args.header
    Ref.file = os.path.join(inputs["out_prefix"], Ref.file)
    if args.sample == inputs["sample_name"]:
        if inputs["read1"] is not None:
            inputs["sample_name"] = os.path.basename(os.path.splitext(inputs["read1"])[0])
            inputs["sample_name"] = re.sub("_R1.*", "", inputs["sample_name"])
        else:
            inputs["sample_name"] = os.path.basename(os.path.splitext(inputs["assembly"])[0])
    else:
        inputs["sample_name"] = args.sample

    return inputs


def make_output_directory(inputs: dict):
    """Makes the output directory
    
    Parameters
    ----------
    inputs: dict
        Run settings

    Returns
    -------
    None
        Output directory is created; if there are errors the program stops here 
    """
    if os.path.isdir(inputs["out_prefix"]):
        if inputs["overwrite"]:
            inputs["logging_buffer_message"] += "Output directory exists, removing the existing directory\n"
            try:
                shutil.rmtree(inputs["out_prefix"])
            except PermissionError:
                print("Failed to remove the existing directory. Do you have write permissions?")
                sys.exit(1)
            os.mkdir(inputs["out_prefix"])
            inputs["logging_buffer_message"] += f"New output directory created\n"
        else:
            print(f"Output directory '{inputs['out_prefix']}' exists and overwrite is turned off. Exiting")
            sys.exit(1)
    else:
        os.mkdir(inputs["out_prefix"])
        inputs["logging_buffer_message"] += f"New output directory created\n"


def configure_logger(inputs: dict):
    """Configures the logging for printing

    Parameters
    ----------
    inputs: dict
        Run settings

    Returns
    -------
    None
        Logger behavior is set based on the Inputs variable
    """
    try:
        logging.basicConfig(filename=inputs["log"], filemode="w", level=logging.DEBUG,
                            format=f"[%(asctime)s | {inputs['out_prefix']} ]  %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p")
    except FileNotFoundError:
        print(
            f"The supplied location for the log file '{inputs['log']}' doesn't exist. Please check if the location exists.")
        sys.exit(1)
    except IOError:
        print(
            f"I don't seem to have access to make the log file. Are the permissions correct or is there a directory with the same name?")
        sys.exit(1)

    if inputs["verbose"]:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter(f"[%(asctime)s | {inputs['out_prefix']} ] %(message)s")
        console.setFormatter(formatter)
        logging.getLogger().addHandler(console)


def get_inputs(inputs: dict):
    """Writes run settings to log

    Parameters
    ----------
    inputs: dict
        Run settings

    """
    
    logging.debug(f"""Parameters input:
                Read 1      {inputs["read1"]}
                Read 2      {inputs["read2"]}
                Assembly    {inputs["assembly"]}
                Threads     {inputs["threads"]}
                Out Prefix  {inputs["out_prefix"]}
                Log file    {inputs["log"]}
                SBT         {inputs["sbt"]}
                Suffix      {inputs["suffix"]}
                Profile     {inputs["profile"]}
                Verbose     {inputs["verbose"]}
                Overwrite   {inputs["overwrite"]}
                Depth       {inputs["depth"]}\n\n""")


def check_program(program_name: str, inputs: dict) -> None:
    """Checks if the supplied program_name exists

    Parameters
    ----------
    program_name : str
        Name of the program to check if it exists
    inputs: dict
        Run settings

    Returns
    -------
    None
        Exits the program if a dependency doesn't exist
    """
    logging.info(f"Checking for program {program_name}")
    path = shutil.which(program_name)
    if path is None:
        logging.critical(f"Program {program_name} not found! Cannot continue; dependency not fulfilled.")
        if not inputs["verbose"]:
            print(f"Program {program_name} not found! Cannot continue; dependency not fulfilled.")
        sys.exit(1)

    if "blast" in program_name:
        # only 1 '-'' for blast arguments
        command = f"{program_name} -version | head -1"
    elif program_name == "isPcr":
        # No version retrieval option. pull from usage head instead.
        command = 'isPcr 2>&1 | head -1 | grep -Po "(?<=v )\S+"'
    elif program_name == "samtools":
        # only want the first two lines
        command = "samtools --version | head -1"
    else:
        command = f"{program_name} --version"
    version = result = subprocess.check_output(command, shell=True, stderr=subprocess.DEVNULL, encoding='utf-8')
    logging.info(f"{program_name} version is {version.strip()}")

def check_files(inputs: dict) -> None:
    """Checks if all the input files exists; exits if file not found or if file is a directory
    
    Parameters
    ----------
    inputs: dict
        Run settings

    Returns
    -------
    None
        Exits the program if file doesn't exist
    """
    if inputs["read1"] and not os.path.isfile(inputs["read1"]):
        logging.critical(f"Read file 1: '{inputs['read1']}' doesn't exist. Exiting")
        if not inputs["verbose"]:
            print(f"Read file 1: '{inputs['read1']}' doesn't exist. Exiting")
        sys.exit(1)
    if inputs["read2"] and not os.path.isfile(inputs['read2']):
        logging.critical(f"Read file 2: '{inputs['read2']}' doesn't exist. Exiting")
        if not inputs["verbose"]:
            print(f"Read file 2: '{inputs['read2']}' doesn't exist. Exiting")
        sys.exit(1)
    if inputs["assembly"] and not os.path.isfile(inputs["assembly"]):
        logging.critical(f"Assembly file: '{inputs['assembly']}' doesn't exist. Exiting")
        if not inputs["verbose"]:
            print(f"Assembly file: '{inputs['assembly']}' doesn't exist. Exiting")
        sys.exit(1)
    if not os.path.isdir(inputs["sbt"]):
        logging.critical(f"SBT directory: '{inputs['sbt']}' doesn't exist. Exiting")
        if not inputs["verbose"]:
            print(f"SBT directory: '{inputs['sbt']}' doesn't exist. Exiting")
        sys.exit(1)
    for locus in Ref.locus_order:
        file = os.path.join(inputs["sbt"], locus + inputs["suffix"])
        if not os.path.isfile(file):
            logging.critical(f"Allele file: '{file}' for locus ({locus}) doesn't exist. Exiting")
            sys.exit(1)
    if not os.path.isfile(inputs["profile"]):
        logging.critical(f"Profile file: '{inputs['profile']}' doesn't exist. Exiting")
        if not inputs["verbose"]:
            print(f"Profile file: '{inputs['profile']}' doesn't exist. Exiting")
        sys.exit(1)


def ensure_safe_threads(inputs: dict, threads: int = 1) -> dict:
    """Ensures that the number of user supplied threads doesn't exceed system capacity

    Sets the number of threads to maximum available threads if threads exceed system capacity.

    Parameters
    ----------
    inputs: dict
        Run settings
    threads : int
        Number of threads supplied by the user.

    Returns
    -------
    dict
        modified inputs object with threads adjusted

    """
    if threads > multiprocessing.cpu_count():
        logging.critical("User has supplied more threads than processor capacity, resetting to max cores.")
        inputs["threads"] = multiprocessing.cpu_count()

    return inputs


def run_command(command: str, tool: str = None, stdin: str = None, shell: bool = False, desc_file: str = None, desc_header: str = None, column_headers: str = "") -> str:
    """Runs a command, and logs it nicely

    Wraps around logging and subprocess.check_output

    Parameters
    ----------
    command : str
        The command to be executed. Converted to a list internally.

    tool: str
        The name of the tool for logging purposes

    stdin: str, optional
        Optional text to be supplied as stdin to the command

    shell: bool, optional
        shell option passed to check_output

    desc_file: str, optional
        filename to write output of command 

    desc_header: str, optional
        String to write before command output (i.e., narrative description of command)

    column_headers: str, optional
        Line of text to prepend to command output in log and intermediate_outputs.txt

    Returns
    -------
    str
        The output generated by running the command

    """
    logging.debug(f"Running command: {command}")
    full_command = command
    if tool is not None:
        logging.info(f"Running {tool}")
    if not shell:
        command = shlex.split(command, posix=False)
    if stdin is not None:
        try:
            result = subprocess.check_output(command, stderr=subprocess.DEVNULL, shell=shell,
                                             input=stdin, encoding='utf-8')
        except subprocess.CalledProcessError:
            logging.critical(f"CRITICAL ERROR! The following command had an improper exit: \n{full_command}\n")
            sys.exit(1)
    else:
        try:
            result = subprocess.check_output(command, shell=shell, stderr=subprocess.DEVNULL, encoding='utf-8')
        except subprocess.CalledProcessError:
            logging.critical(f"CRITICAL ERROR! The following command had an improper exit: \n{full_command}\n")
            sys.exit(1)

    pretty_result = prettify("\n".join([column_headers, result]))
    if tool is not None:
        logging.debug(f"Command log for {tool}:\n{pretty_result}")
        logging.info(f"Finished running {tool}")
    else:
        logging.debug(f"Command log:\n{pretty_result}")

    # Write result to desc_file
    if desc_file:
        with open(desc_file, 'a') as f:
            if desc_header:
                f.write(desc_header + '\n\n')
            f.write(f"Command:\n\n{full_command}\n\n")
            f.write(f"Output:\n\n{pretty_result}\n\n")

    return result


def prettify(text, delim="\t"):
    lines = text.split("\n")
    # calculate padding for each column
    col_pad = defaultdict(int)
    for line in lines:
        for n, col in enumerate(line.split(delim)):
            col_pad[n] = max(col_pad[n], len(col))

    # Build new, padded output
    padded_text = ""
    for line in lines:
        line_new = []
        for n, col in enumerate(line.split(delim)): 
            line_new.append(col.ljust(col_pad[n]+2))
        padded_text += "".join(line_new) + "\n"

    return padded_text

def blast_momps_allele(seq: str, db: str) -> list:
    """BLAST the mompS allele in the isolate to find the allele number

    Parameters
    ----------
    seq : str
        mompS allele sequence found in the isolate

    db : str, optional
        location of the mompS alleles database file

    Returns
    -------
    str
        mompS allele number
    """
    logging.debug(f"Looking for \n{seq}")
    blastcmd = f"blastn -query - -db {db} -outfmt '6 std qlen slen sseqid' | awk -F'\\t' '{{OFS=FS}}{{gsub(/_.+/, \"\", $15)}}1' | sort -k15,15 -k12,12gr | sort --merge -u  -k15,15"
    column_headers = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tsseqid"
    res = run_command(blastcmd, "blastn/mompS", seq, shell=True, column_headers=column_headers).rstrip()
    if res == "":
        return [Allele()]
    else:
        bits = res.split()
        a = Allele()
        if float(bits[2]) == 100.00 and bits[3] == bits[13]:
            a.allele_id = bits[1].split("_")[-1]
        else:
            a.allele_id = bits[1].split("_")[-1]+"*"
        a.seq = "".join(seq.split("\n")[1:])[93:445].upper()

        return [a]


def call_momps_pcr(inputs: dict, assembly_file: str) -> list:
    """Find the mompS gene using an in silico PCR procedure

    Parameters
    ----------
    inputs: dict
        Run settings
    assembly_file : str, optional
        Read1 file name

    Returns
    -------
    str
        mompS allele number
    """
    primer1 = os.path.join(inputs["sbt"], "mompS_primer1.tab")
    ispcr_command = f"isPcr {assembly_file} {primer1} {Ref.ispcr_opt}"
    primer1_res = run_command(ispcr_command, "mompS2 primer1")

    if primer1_res != "":
        # nested PCR
        primer2 = os.path.join(inputs["sbt"], "mompS_primer2.tab")
        ispcr_command = f"isPcr stdin {primer2} {Ref.ispcr_opt}"
        primer2_res = run_command(ispcr_command, "mompS2 primer2", primer1_res).rstrip()
        logging.debug(f"Found the sequence: {primer2_res}")
        a_list = []
        for pcr_res in primer2_res.split(">")[1:]:
            a_list += blast_momps_allele(seq=">"+pcr_res, db=os.path.join(inputs["sbt"], "mompS" + inputs["suffix"]))
        identified_mompS_msg = f"mompS alleles identified: {', '.join([a.allele_id for a in a_list])}"
        with open(f"{inputs['out_prefix']}/intermediate_outputs.txt", 'a') as f:
            f.write(identified_mompS_msg + "\n\n")
        logging.info(identified_mompS_msg)
        return a_list 
    else:
        # try BLAST without isPCR
        error_msg = f"in silico PCR returned no results. The region around mompS may be missing from your assembly"
        logging.info(error_msg)
        with open(f"{inputs['out_prefix']}/intermediate_outputs.txt", 'a') as f:
            f.write(error_msg + '\n\n')
        blast_command = f"blastn -query {assembly_file} -db {inputs['sbt']}/mompS_alleles.tfa -outfmt '6 std qlen slen sseqid sseq' | awk -F'\\t' '{{OFS=FS}}{{gsub(/_.+/, \"\", $15)}}1' | sort -k15,15 -k12,12gr | sort --merge -u  -k15,15"
        desc_header = "Best match of mompS in provided assembly using BLASTN."
        column_headers = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tsseqid\tsseq"
        result = run_command(blast_command, tool='blast', shell=True, desc_file=f"{inputs['out_prefix']}/intermediate_outputs.txt", desc_header=desc_header, column_headers=column_headers)
        if result != "":
            error_msg = f"BLAST identified mompS"
            logging.info(error_msg)
            with open(f"{inputs['out_prefix']}/intermediate_outputs.txt", 'a') as f:
                f.write(error_msg + '\n\n')
            bits = result.strip().split()
            a = Allele()
            if float(bits[2]) == 100.00 and bits[3] == bits[13]:
                a.allele_id = bits[1].split("_")[-1]
            else:
                a.allele_id = bits[1].split("_")[-1]+"*"

            # Extract sequence of allele from assembly
            ass_contig = bits[0]
            ass_start = int(bits[6])
            ass_end = int(bits[7])
            db_start = int(bits[8])
            db_end = int(bits[9])

            a.seq = bits[15]
            
            return [a]

        else:
            error_msg = f"BLAST and in silico PCR returned no results. mompS may be missing from your assembly"
            logging.info(error_msg)
            with open(f"{inputs['out_prefix']}/intermediate_outputs.txt", 'a') as f:
                f.write(error_msg + '\n\n')
            return [Allele()]


def blast_non_momps(inputs: dict, assembly_file: str, ref: Ref) -> dict:
    """Find the rest of alleles (non-mompS) by BLAST search

    Parameters
    ----------
    inputs: dict
        Run settings
    assembly_file : str, optional
        Read1 file name
    ref: Ref
        Information about reference sequence

    Returns
    -------
    dict
        dictionary containing locus (key) to allele (value) mapping
    """
    loci = ["flaA", "pilE", "asd", "mip", "proA", "neuA_neuAH"]
    calls = dict.fromkeys(loci, '')

    assembly_dict = fasta_to_dict(assembly_file)

    blast_command = f"blastn -query {assembly_file} -db {inputs['sbt']}/all_loci.fasta -outfmt '6 std qlen slen sseqid' | awk -F'\\t' '{{OFS=FS}}{{gsub(/_.+/, \"\", $15)}}1' | sort -k15,15 -k12,12gr | sort --merge -u  -k15,15"
    desc_header = "Best match of each locus in provided assembly using BLASTN."
    column_headers = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tsseqid"
    result = run_command(blast_command, tool='blast', shell=True, desc_file=f"{inputs['out_prefix']}/intermediate_outputs.txt", desc_header=desc_header, column_headers=column_headers)

    for line in result.strip().split('\n'):
        bits = line.split()

        # Find best match in db
        locus = "_".join(bits[1].split("_")[:-1])
        a = Allele()
        if float(bits[2]) == 100.00 and bits[3] == bits[13]:
            a.allele_id = bits[1].split("_")[-1]
        else:
            a.allele_id = bits[1].split("_")[-1]+"*"

        # Extract sequence of allele from assembly
        ass_contig = bits[0]
        ass_start = int(bits[6])
        ass_end = int(bits[7])
        db_start = int(bits[8])
        db_end = int(bits[9])
        
        if db_start < db_end:
            a.seq = assembly_dict[ass_contig][ass_start-1:ass_end]
        else:
            a.seq = rev_comp(assembly_dict[ass_contig][ass_start-1:ass_end])

        calls[locus] = [a]

    not_found_loci = [k for k,v in calls.items() if v =='']

    if len(not_found_loci) != 0:
        error_msg = f"The following loci were not found in your assembly: {', '.join(not_found_loci)}\n"
        logging.info(error_msg)
        with open(f"{inputs['out_prefix']}/intermediate_outputs.txt", 'a') as f:
            f.write(error_msg + '\n\n')
        for locus in not_found_loci:
            calls[locus] = [Allele()]

    return calls

def get_st(allele_profile: str, Ref: Ref, profile_file: str) -> str:
    """Looks for the ST in the allele profile table (simple look-up)

    Parameters
    ----------
    allele_profile : str, optional
        allele profile as ordered tab-separated string
    Ref: Ref
        Reference sequence information
    profile_file : str, optional
        profile file containing ST as the first column, and allele profiles in the next columns

    Returns
    -------
    str
        ST number
    """
    novel_ST = "Novel ST" # all 7 target genes were found, but not present in the profile - probably a novel sequence type
    novel_allele = "Novel ST*" # one or multiple target genes have a novel allele found
    not_found = "NF-" # one or more of the target genes were unidentifiable - ST is also unidentifiable as a result
    multiple = "NF?" # one or more of the target genes have multiple alleles found - ST is ambiguous due to multiple alleles
    
    if "-" in allele_profile:
        return not_found
    elif "?" in allele_profile:
        return multiple
    elif "*" in allele_profile:
         return novel_allele
      
    with open(profile_file, "r") as f:
        f.readline()
        for line in f:
            line = line.rstrip()
            if line.endswith(allele_profile):
                st = line.split("\t")[0]
                return st
    return novel_ST


def read_sam_file(samfile: str):
    contig_dict = defaultdict(list) #{'contig': [SAM_data, SAM_data]}
    read_info_dict = defaultdict(list) #{'readname': [SAM_data, SAM_data]}

    with open(samfile, 'r') as sam:
        for line in sam.readlines():
            if line[0] == "@":
                continue
            
            entry = SAM_data(line)
            read_info_dict[entry.qname].append(entry)
            contig_dict[entry.rname].append(entry)

    return contig_dict, read_info_dict


def assess_allele_conf(bialleles, reads_at_locs, allele_idxs, read_info_dict, ref):
    """Checks mompS PCR primer presence absence to assess allele location at native locus

    Parameters
    ----------
    bialleles : list
        Bases calls for each alleles at each multiallelic site
    reads_at_locs: list
        List of reads that map to each multiallelic site
    allele_idxs: list
        Location of each allele within gene of interest
    read_info_dict: dict
        {read_name : [SAM_data, SAM_data]}
        information for all read pairs
    ref: Ref
        Information about reference sequence

    Returns
    -------
    list
        list of Allele class instances describing allele characteristics
    """

    # Make 2 empty Allele instances to store allele info
    alleles_info = [Allele(len(bialleles[0]), bialleles[0]), Allele(len(bialleles[0]), bialleles[1])]

    for n, reads in enumerate(reads_at_locs):
        idx = allele_idxs[n]
        for read in reads:
            mates = read_info_dict[read]
            allele_found = False
            for mate in mates:
                # Get base call at index.
                # Index result with 0 as dict returned with k,v per index given
                b = mate.get_base_calls([idx + ref.allele_start-1])[0]
                for allele in alleles_info:
                    if b == allele.basecalls[n]:
                        allele.reads_at_locs[n].append(read)
                        allele_found = True
                        break
                if allele_found:
                    break

    for allele in alleles_info:
        all_informative_reads = set(allele.reads_at_locs[0])
        if len(allele.reads_at_locs) > 1:
            for i in range(1, len(allele.reads_at_locs)):
                all_informative_reads = all_informative_reads.union(set(allele.reads_at_locs[i]))
        all_informative_reads = sorted([i for i in all_informative_reads])

        for read_pair in all_informative_reads:
            for mate in read_info_dict[read_pair]:
                if mate.flag > 2047:
                    continue
                
                if "GAAGTCCGGCTGGATAATTTATCCA" in mate.seq:
                    if 16 & mate.flag == 16: 
                        # If 16 in the flag
                        # i.e., read containing primer mapped in reverse
                        # meaning mate is 5' of primer
                        allele.confidence["for"] += 1
                    
                    else: # Primer indicates secondary allele
                        allele.confidence['against'] += 1
                    
                    break

    
    return alleles_info


def process_reads(contig_dict: dict, read_info_dict: dict, ref: Ref, outdir: str):

    # Check coverage of neuA/neuAh regions to infer which is present

    logging.info("Processing sam file for downstream analysis.")
    process_sam_command = f"samtools view -Sb {outdir}/reads_vs_all_ref_filt.sam | samtools sort - > {outdir}/reads_vs_all_ref_filt_sorted.bam; samtools index {outdir}/reads_vs_all_ref_filt_sorted.bam"
    run_command(process_sam_command, tool='samtools', shell=True)

    logging.info("Checking coverage of reference loci by mapped reads")
    coverage_command = f"samtools coverage -r flaA:351-531 {outdir}/reads_vs_all_ref_filt_sorted.bam | cut -f 1,4,5,6,7,8,9; samtools coverage -r pilE:351-682 {outdir}/reads_vs_all_ref_filt_sorted.bam | tail -1 | cut -f 1,4,5,6,7,8,9; samtools coverage -r asd:351-822 {outdir}/reads_vs_all_ref_filt_sorted.bam | tail -1 | cut -f 1,4,5,6,7,8,9; samtools coverage -r mip:350-750 {outdir}/reads_vs_all_ref_filt_sorted.bam | tail -1 | cut -f 1,4,5,6,7,8,9; samtools coverage -r mompS:367-717 {outdir}/reads_vs_all_ref_filt_sorted.bam | tail -1 | cut -f 1,4,5,6,7,8,9; samtools coverage -r proA:350-754 {outdir}/reads_vs_all_ref_filt_sorted.bam | tail -1 | cut -f 1,4,5,6,7,8,9; samtools coverage -r neuA:350-702 {outdir}/reads_vs_all_ref_filt_sorted.bam | tail -1 | cut -f 1,4,5,6,7,8,9; samtools coverage -r neuAh:350-702 {outdir}/reads_vs_all_ref_filt_sorted.bam | tail -1 | cut -f 1,4,5,6,7,8,9; samtools coverage -r neuA_207:350-702 {outdir}/reads_vs_all_ref_filt_sorted.bam | tail -1 | cut -f 1,4,5,6,7,8,9; samtools coverage -r neuA_211:350-702 {outdir}/reads_vs_all_ref_filt_sorted.bam | tail -1 | cut -f 1,4,5,6,7,8,9"
    desc_header = "Assessing coverage of MLST loci by provided sequencing reads."

    result = run_command(coverage_command, tool='samtools coverage', shell=True, desc_file=f"{outdir}/intermediate_outputs.txt", desc_header=desc_header)

    alleles = {}
    cov_results = {}
    for line in result.strip().split('\n')[1:]:
        gene, _, _, cov, depth, _, _ = line.split()
        cov = float(cov)
        depth = float(depth)
        cov_results[gene] = {'cov': cov, 'depth': depth}
        if cov != 100.:
            if 'neuA' in gene:
                if cov < 99:
                    del ref.REF_POSITIONS[gene]
            else:
                logging.info(f"Insufficient coverage of the {gene} locus to identify allele. This may indicate a gene deletion or a bad sequencing run.")
                a = Allele()
                a.allele_id = '-'
                alleles[gene] = [a]
                del ref.REF_POSITIONS[gene]

    cov_msg = ""

    if len([i for i in ref.REF_POSITIONS.keys() if 'neuA' in i]) == 0:
        # No neuA loci had sufficient coverage and were deleted
        logging.info(f"WARNING! Insufficient coverage of neuA to identify the allele. This may indicate a gene deletion or a bad sequencing run.")
        cov_msg += "minimum coverage of neuA locus is 0." + "\n"
        a = Allele()
        a.allele_id = '-'
        alleles['neuA_neuAH'] = [a]

    if len([i for i in ref.REF_POSITIONS.keys() if 'neuA' in i]) > 1:
        cov_sorted = sorted(
            [(k,v['depth']) for k,v in cov_results.items()],
            key=lambda x: x[1],
            reverse=True)
        # Try to use depth to determine the right one to use
        # if highest depth locus more than 3* the depth of next highest, keep
        if cov_sorted[0][1] > 3*cov_sorted[1][1]:
            logging.info(f"mean coverage of {cov_sorted[0][0]} is greater. Using {cov_sorted[0][0]}.")
            for gene in cov_sorted[1:]:
                del ref.REF_POSITIONS[gene[0]]
        else:
            logging.info("Analysis of read mapping to neuA locus variants was unsuccessful. Unclear which to use.")

    
    
    for locus in ref.REF_POSITIONS:
        locus_reads = contig_dict[locus]
        allele_start = ref.REF_POSITIONS[locus]['start_pos']
        allele_stop = ref.REF_POSITIONS[locus]['end_pos']

        # reads_dict {bases: {pos:base}, qualities: {pos:qual}, readnames: {pos:rname}}
        reads_dict ={
            'bases' : {k:[] for k in range(allele_start-1, allele_stop)},
            'qualities' : {k:[] for k in range(allele_start-1, allele_stop)},
            'readnames' : {k:[] for k in range(allele_start-1, allele_stop)}
            }

        for read in locus_reads:
            if ( 
                (read.pos <= allele_start and allele_start <= read.end) or
                (read.pos <= allele_stop and allele_stop <= read.end) or
                (allele_start <= read.pos and read.end <= allele_stop) or
                (read.pos <= allele_start and allele_stop <= read.end)
                ):
                if read.pos < allele_start:
                    pad = allele_start - read.pos
                else:
                    pad = 0
                

                for read_idx in range(pad, min([read.ln - 1, allele_stop+1-read.pos])):
                    ref_seq_idx = read.pos + read_idx - 1
                    reads_dict['bases'][ref_seq_idx].append(read.seq[read_idx])
                    reads_dict['qualities'][ref_seq_idx].append(read.qual[read_idx])
                    reads_dict['readnames'][ref_seq_idx].append(read.qname)

        seq = []
        cov = []
        for position in range(allele_start-1, allele_stop):
            good_basecalls = [
                base for base, qual in zip(
                    reads_dict['bases'][position],
                    reads_dict['qualities'][position]
                    )
                if ord(qual)-33 > 20 # Assuming Phred 33...
                ]
            count = Counter(good_basecalls)
            cov.append(sum([n for n in count.values()]))
            if len(count) == 1:
                seq.append([[b for b in count][0]])
            else:
                total = sum([i for i in count.values()])
                seq.append(
                    [base for base, num in count.items() if num > 0.3*total]
                    )
        min_cov = min(cov)

        if min_cov == 0:
            msg = f"WARNING: After applying a quality cutoff of 20 to basecalls, at least one position in {locus.split('_')[0]} has 0 coverage and can't be resolved"
            logging.info(msg)
            cov_msg += f"\n{msg}\n\n"
            a = Allele()
            a.allele_id = '?'
            alleles[locus] = [a]
            continue

        cov_msg += f"minimum coverage of {locus} locus is {min_cov}." + '\n'
        logging.info(f"minimum coverage of {locus} locus is {min_cov}.")

        num_alleles_per_site = [len(i) for i in seq]

        n_multiallelic = len([i for i in num_alleles_per_site if i > 1])

        if n_multiallelic == 0:
            alleles[locus] = [Allele()]
            alleles[locus][0].seq = "".join([b[0] for b in seq])
            if locus == 'mompS':
                # Don't assess confidence if only one allele
                alleles[locus][0].confidence = {'for': 'NA', 'against': 'NA'}

        elif n_multiallelic == 1:
            multi_allelic_idx = [n for n,i in enumerate(num_alleles_per_site) if i == 2]
            reads_at_a = reads_dict['readnames'][multi_allelic_idx[0]+allele_start-1]
            bialleles = seq[multi_allelic_idx[0]]
            
            if locus == 'mompS':
                alleles[locus] = assess_allele_conf(bialleles, [reads_at_a], multi_allelic_idx, read_info_dict, ref)
            else:
                alleles[locus] = [Allele(len(bialleles[0]), bialleles[i]) for i in range(len(bialleles))]

            for base in seq:
                for allele in alleles[locus]:
                    if len(base) == 1:
                        allele.seq += base[0]
                    else:
                        allele.seq += allele.basecalls[0]

        else:
            multi_allelic_idx = [n for n,i in enumerate(num_alleles_per_site) if i == 2]
            reads_at_locs = [reads_dict['readnames'][idx+allele_start-1] for idx in multi_allelic_idx]

            intersect = set(reads_at_locs[0]).intersection(set(reads_at_locs[1]))
            for i in range(2, len(reads_at_locs)):
                intersect = intersect.intersection(set(reads_at_locs[i]))

            intersect = sorted([i for i in intersect])

            conflicting_reads = []
            read_pair_base_calls = []
            for read_name in intersect:
                read_pair = read_info_dict[read_name]
                calls_list = []
                for mate in read_pair: # If read is supplementary alignment, ignore.
                    if mate.flag > 2047:
                        continue
                    calls_list.append(mate.get_base_calls(
                        [idx+allele_start-1 for idx in multi_allelic_idx]
                        ))

                # Assess base calls of mate reads that map to one or more positions
                agreeing_calls = []
                for k in range(len(multi_allelic_idx)):
                    calls = []
                    for mate_calls in calls_list:
                        if mate_calls[k] != '':
                            calls.append(mate_calls[k])
                    if len(set(calls)) > 1:
                        conflicting_reads.append(calls_list)
                        break
                    else:
                        agreeing_calls.append(calls[0])
                if len(agreeing_calls) != 0:
                    read_pair_base_calls.append("".join(agreeing_calls))


            if len(conflicting_reads) > 0.33 * len(reads_at_locs[0]):
                logging.info(f"more than 33% of reads disagree with which variant bases are in the same gene for {locus.split('_')[0]}")
                cov_message += f"\nmore than 33% of reads disagree with which variant bases are in the same gene for {locus.split('_')[0]}\n\n"
                a = Allele()
                a.allele_id = '?'
                alleles[locus] = [a]
                continue

            biallele_results_count = Counter(read_pair_base_calls)
            max_biallele_count = max([v for v in biallele_results_count.values()])

            bialleles = [k for k,v in biallele_results_count.items() if v > 0.66*max_biallele_count]

            # If more than 2 alleles found, can't resolve
            if len(bialleles) > 2:
                logging.info(f"ERROR: {len(bialleles)} well-supported {locus.split('_')[0]} alleles identified and can't be resolved. Aborting.")
                with open(f"{outdir}/intermediate_outputs.txt", 'a') as f:
                    f.write(f"\nERROR: {len(bialleles)} well-supported {locus.split('_')[0]} alleles identified and can't be resolved. Aborting.\n\n")
                sys.exit(1)
            if len(bialleles) > 1 and locus == 'mompS':
                alleles[locus] = assess_allele_conf(bialleles, reads_at_locs, multi_allelic_idx, read_info_dict, ref)

            else:
                alleles[locus] = [Allele(len(bialleles[0]), bialleles[i]) for i in range(len(bialleles))]

            biallic_count = 0
            for base in seq:
                for allele in alleles[locus]:
                    if len(base) == 1:
                        allele.seq += base[0]
                    else:
                        allele.seq += allele.basecalls[biallic_count]
                if len(base) > 1:
                    biallic_count+=1
    
    with open(f"{outdir}/intermediate_outputs.txt", 'a') as f:
        f.write(cov_msg + "\n\n")

    # Remove deletions from seq
    for allele_list in alleles.values():
        for allele in allele_list:
            allele.seq = allele.seq.replace("*", "")

    # Assess confidence in mompS alleles
    for allele in alleles['mompS']:
        allele.assess_conf()

    neuAs = [l for l in alleles if 'neuA' in l and l != 'neuA_neuAH']
    if len(neuAs) > 1:
        neuA_list = []
        for neu in neuAs:
            neuA_list += alleles[neu]
            del alleles[neu]
    # rename neuA and neuAH to neuA_neuAH
    elif len(neuAs) == 1:
        alleles['neuA_neuAH'] = alleles[neuAs[0]]
        del alleles[neuAs[0]]

    return alleles


def write_alleles_to_file(alleles: list, outdir: str):
    identified_allele_fasta_string = ""
    for locus, l in alleles.items():
        if len(l) > 1:
            for n, allele in enumerate(l):
                if locus == "mompS":
                    allele.fasta_header = f"{locus}_{n+1}{allele.location}"
                else:
                    allele.fasta_header = f"{locus}_{n+1}"
                identified_allele_fasta_string += f">{allele.fasta_header}\n{allele.seq}\n"
        else:
            allele = l[0]
            allele.fasta_header = f"{locus}"
            identified_allele_fasta_string += f">{allele.fasta_header}\n{allele.seq}\n"

    with open(f"{outdir}/identified_alleles.fna", "w") as fout:
        fout.write(identified_allele_fasta_string)

def map_alleles(inputs: dict, ref: Ref):
    """Map reads to reference loci sequences and identify one or more alleles

    Parameters
    ----------
    inputs: dict
        run settings
    ref: Ref
        Reference sequence information in Ref instance
    
    Returns
    -------
    list
        Identified alleles
    """
    r1 = inputs['read1']
    r2 = inputs['read2']
    threads = inputs['threads']
    outdir = inputs['out_prefix']
    db = inputs['sbt']
    sample_name = inputs['sample_name']
    profile = inputs['profile']

    # Run BWA mem
    logging.info("Mapping reads to reference sequence, then filtering unmapped reads from sam file")
    mapping_command = f"minimap2 -ax sr -t {threads} {db}/ref_gene_regions.fna {r1} {r2} | samtools view -h -F 0x4 -@ {threads} -o {outdir}/reads_vs_all_ref_filt.sam"
    run_command(mapping_command, tool='minimap2 -ax sr', shell=True)

    contig_dict, read_info_dict = read_sam_file(f"{outdir}/reads_vs_all_ref_filt.sam")

    alleles = process_reads(contig_dict, read_info_dict, ref, outdir)
    write_alleles_to_file(alleles, outdir)
    # BLAST alleles
    logging.info("BLASTing identified alleles against database")
    blast_command = f"blastn -query {outdir}/identified_alleles.fna -db {db}/all_loci.fasta -outfmt '6 std qlen slen' | sort -k1,1 -k12,12gr | sort --merge -u  -k1,1"
    desc_header = "Best match of each identified sequence determined using BLASTN"
    column_headers = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen"

    result = run_command(blast_command, tool='blast', shell=True, desc_file=f"{outdir}/intermediate_outputs.txt", desc_header=desc_header, column_headers=column_headers)
    if len(result) == 0:
        logging.info("WARNING: No allele matches found in the database. Can't resolve any alleles!")
        with open(f"{outdir}/intermediate_outputs.txt", "a") as fout:
            fout.write("WARNING: No allele matches found in the database. Can't resolve any alleles!")
        return alleles
        
    for line in result.strip().split('\n'):
        bits = line.split()
        if float(bits[2]) == 100.00 and bits[3] == bits[12] and bits[3] == bits[13]:
            for allele_list in alleles.values():
                for allele in allele_list:
                    if bits[0] == allele.fasta_header:
                        allele.allele_id = bits[1].split("_")[-1]
        else:
            for allele_list in alleles.values():
                for allele in allele_list:
                    if bits[0] == allele.fasta_header:
                        allele.allele_id = bits[1].split("_")[-1]+"*"

    if len(alleles['mompS']) == 1:
        logging.info("1 mompS allele identified.")
        message = "Identified allele information:\n\n1 mompS allele identified.\n\n"
    else:
        logging.info(f"Identified allele information:\n\n{len(alleles['mompS'])} mompS allele identified.")
        message = f"{len(alleles['mompS'])} mompS allele identified.\n"
        for a in alleles['mompS']:
            logging.info(f"mompS allele '{a.allele_id}' information")
            logging.info(f"lowest coverage of bialleleic site: {min([len(set(i)) for i in a.reads_at_locs])}")
            logging.info(f"number of reads from this allele containing outermost reverse primer sequence: {a.confidence['for']}")
            logging.info(f"number of reads from this allele containing outermost reverse primer sequence in the reverse orientation (indicating this is the secondary allele): {a.confidence['against']}")
            message += f"mompS allele '{a.allele_id}' information\n"
            message += f"lowest coverage of bialleleic site: {min([len(set(i)) for i in a.reads_at_locs])}\n"
            message += f"number of reads from this allele containing outermost reverse primer sequence: {a.confidence['for']}\n"
            message += f"number of reads from this allele containing outermost reverse primer sequence in the reverse orientation (indicating this is the secondary allele): {a.confidence['against']}\n\n"
    old_mompS = alleles['mompS']
    if len(alleles['mompS']) > 1:
        which_native = ["_native_locus" in a.location for a in alleles['mompS']]
        if not any(which_native):
            logging.info("Unable to determine which allele is present in native mompS locus")

        else:
            if len([i for i in which_native if i]) > 1:
                logging.info("Found evidence that multiple alleles may exist in a sequence context that is similar to the native locus. Trying to determine allele locations.")
 
            else:
                native_allele = [a for a in alleles['mompS'] if "_native_locus" in a.location][0]
                non_native_alleles = [a for a in alleles['mompS'] if "_native_locus" not in a.location]
                logging.info(f"Allele {native_allele.allele_id} was determined to be the primary mompS allele. {native_allele.confidence['for']} reads support this.")
                for a in non_native_alleles:
                    logging.info(f"Allele {a.allele_id} was determined not to be the secondary mompS allele.")

        # Pick primary and secondary allele
        # If ANY reads contain primer in correct orientation
        # That is evidence of primary allele
        alleles['mompS'] = [a for a in alleles['mompS'] if a.confidence['for'] > 0]
        if len(alleles['mompS']) > 1:
            if alleles['mompS'][0].confidence['for'] > 3*alleles['mompS'][1].confidence['for']:
                alleles['mompS'] = [alleles['mompS'][0]]
            elif alleles['mompS'][1].confidence['for'] > 3*alleles['mompS'][0].confidence['for']:
                alleles['mompS'] = [alleles['mompS'][1]]
            else:
                logging.info("Failed to determine primary mompS allele as both alleles appear to be flanked by primer in the expected orientation.")
                message += "Failed to determine primary mompS allele as both alleles appear to be flanked by primer in the expected orientation.\n\n"
                # save mompS info for detailed outfile
                alleles['mompS'] = [Allele()]
                alleles['mompS'][0].allele_id = '?'
        elif len(alleles['mompS']) == 0:
            # If only one allele has reads with primer in the wrong orientation, use the other one
            alleles['mompS'] = [a for a in alleles['mompS'] if a.confidence['against'] == 0]
            if len(alleles['mompS']) != 1:
                logging.info("Failed to determine primary mompS allele. Primary mompS allele is identified by finding read pairs that cover both biallelic positions and sequencing primer. In this sample, no such reads were found. Perhaps sequencing reads are too short.")
                message += "Failed to determine primary mompS allele. Primary mompS allele is identified by finding read pairs that cover both biallelic positions and sequencing primer. In this sample, no such reads were found. Perhaps sequencing reads are too short.\n\n"
                alleles['mompS'] = [Allele()]
                alleles['mompS'][0].allele_id = '?'

    temp_alleles = alleles.copy()
    temp_alleles['mompS'] = old_mompS
    write_possible_mlsts(inputs=inputs, alleles=temp_alleles, header=True, confidence=True)

    for locus in [i for i in alleles if i != 'mompS']:
        if len(alleles[locus]) > 1:
            logging.info(f"{len(alleles[locus])} {locus} alleles identified.")
            message += f"{len(alleles[locus])} {locus} alleles identified.\n"
            alleles_found = [a.allele_id for a in alleles[locus]]
            if len(set(alleles_found)) == 1:
                logging.info(f"both identified {locus} alleles are most similar to allele {alleles_found[0]}. That allele will be reported.")
                message += f"both identified {locus} alleles are most similar to allele {alleles_found[0]}. That allele will be reported.\n\n"
            else:
                logging.info(f"The following {locus} alleles were found: {', '.join(alleles_found)}. It is unclear which is the right allele. It may help to perform QC on your reads and rerun el_gato.")
                message += f"The following {locus} alleles were found: {', '.join(alleles_found)}. It is unclear which is the right allele. It may help to perform QC on your reads and rerun el_gato.\n\n"
                a = Allele()
                a.allele_id = '?'
                alleles[locus] = [a]


    with open(f"{outdir}/intermediate_outputs.txt", "a") as fout:
            fout.write(message)

    return alleles


def write_possible_mlsts(inputs: dict, alleles: dict, header: bool, confidence: bool):
    """ Write possible combinations of alleles and corresponding ST to file
    """
    possible_mlsts = ""

    if header:
        if confidence:
            possible_mlsts += "Sample\tST\t" + "\t".join(Ref.locus_order) + "\tmompS_reads_support\t" + "mompS_reads_against\n"
        else:
            possible_mlsts += "Sample\tST\t" + "\t".join(Ref.locus_order) + "\n"

    for flaA in alleles['flaA']:
        for pilE in alleles['pilE']:
            for asd in alleles['asd']:
                for mip in alleles['mip']:
                    for mompS in alleles['mompS']:
                        for proA in alleles['proA']:
                            for neuA_neuAH in alleles['neuA_neuAH']:
                                allele_profile = '\t'.join([a.allele_id for a in [flaA, pilE, asd, mip, mompS, proA, neuA_neuAH]])
                                if confidence:
                                    possible_mlsts += (inputs['sample_name'] + "\t" + get_st(allele_profile, Ref, profile_file=inputs["profile"]) + "\t" + allele_profile + f"\t{mompS.confidence['for']}\t{mompS.confidence['against']}\n")
                                else:
                                     possible_mlsts += (inputs['sample_name'] + "\t" + get_st(allele_profile, Ref, profile_file=inputs["profile"]) + "\t" + allele_profile + "\n")

    with open(f"{inputs['out_prefix']}/possible_mlsts.txt", 'w') as f:
        f.write(possible_mlsts)

def choose_analysis_path(inputs: dict, ref: Ref) -> str:
    """Pick the correct analysis path based on the program input supplied

    Parameters
    ----------
    inputs: dict
        Run settings

    ref: Ref
        Ref class instance

    Returns
    -------
    str
        formatted ST + allele profile (and optional header) of the isolate
    """
    alleles = {}
    if inputs["analysis_path"] == "a":
        alleles = blast_non_momps(inputs, assembly_file=inputs["assembly"], ref=ref)
        alleles["mompS"] = call_momps_pcr(inputs, assembly_file=inputs["assembly"])
        write_possible_mlsts(inputs=inputs, alleles=alleles, header=True, confidence=False)
        write_alleles_to_file(alleles, inputs['out_prefix'])
        for locus, a in alleles.items():
            if len(a) > 1:
                a = Allele()
                a.allele_id = '?'
                alleles[locus] = [a]
    elif inputs["analysis_path"] == "r":
        alleles = map_alleles(inputs=inputs, ref=ref)

    else:
        logging.critical(
            "This path should not have been traversed. Is inputs['analysis_path'] being changed somewhere else?")
        if not inputs["verbose"]:
            print(f"Something went wrong with genome assembly. Please check log.")

    return print_table(inputs,  Ref, alleles)


def print_table(inputs: dict, Ref: Ref, alleles: dict) -> str:
    """Formats the allele profile so it's ready for printing

    Parameters
    ----------
    inputs: dict
        Run settings
    Ref: Ref
        Reference sequence information
    alleles : dict
        The allele profile and the ST

    Returns
    -------
    str
        formatted ST + allele profile (and optional header) of the isolate
    """

    outlines = []
    if inputs['header']:
        header = "Sample\tST\t" + "\t".join(Ref.locus_order)
        if len(alleles['mompS']) > 1:
            header += f"\tmompS_min_coverage\tmompS_reads_with_primer"
        outlines.append(header)
    for n in range(len(alleles['mompS'])):
        allele_profile = ""
        for locus in Ref.locus_order:
            if locus != "mompS":
                allele_profile += alleles[locus][0].allele_id + "\t"
            else:
                allele_profile += alleles[locus][n].allele_id + "\t"
        allele_profile = allele_profile.rstrip()
        allele_profile = (inputs["sample_name"] 
            + "\t" + get_st(allele_profile, Ref,
                            profile_file=inputs["profile"])
            + "\t" + allele_profile)

        outlines.append(allele_profile)

    
    return "\n".join(outlines)


def pretty_time_delta(seconds: int):
    """Pretty print the time

    Parameters
    ----------
    seconds : int
        Time in seconds to convert to human readable

    Returns
    -------
    str
        human readable time
    """
    seconds = int(seconds)
    days, seconds = divmod(seconds, 86400)
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    if days > 0:
        return f"{days}d {hours}h {minutes}m {seconds}s"
    elif hours > 0:
        return f"{hours}h {minutes}m {seconds}s"
    elif minutes > 0:
        return f"{minutes}m {seconds}s"
    else:
        return f"{seconds}s"


def main():
    """ Main code """

    inputs = {
        'read1' : None,
        'read2' : None,
        'assembly' : None,
        'threads' : 1,
        'out_prefix' : "out",
        'sample_name' : "<Inferred from input file>",
        'log' : os.path.join("out", "run.log"),
        'sbt' : os.path.join(os.path.dirname(__file__), "db"),
        'suffix' : "_alleles.tfa",
        'profile' : os.path.join(os.path.dirname(__file__), "db", "lpneumophila.txt"),
        'verbose' : False,
        'overwrite' : False,
        'depth' : 3,
        'analysis_path' : "",
        'logging_buffer_message' : "",
        'header' : True
        }


    parser = get_args()
    args = parser.parse_args()
    inputs = check_input_supplied(args, parser, inputs)
    inputs = set_inputs(args, inputs)
    make_output_directory(inputs)
    configure_logger(inputs)
    logging.info("Starting preprocessing")
    for line in inputs["logging_buffer_message"].rstrip().split("\n"):
        logging.info(line)
    logging.info("Checking if all the prerequisite programs are installed")
    for program in Ref.prereq_programs:
        check_program(program, inputs)
    logging.info("All prerequisite programs are accessible")

    logging.info("Checking if all the required input files exist")
    check_files(inputs)
    logging.info("Input files are present")

    logging.info("Ensuring thread counts are correct")
    inputs = ensure_safe_threads(inputs)
    logging.info("Thread count has been validated")
    logging.info("All reference files have been discovered")
    get_inputs(inputs)
    logging.info("Starting analysis")
    output = choose_analysis_path(inputs, Ref)
    logging.info("Finished analysis")

    logging.debug(f"Output = \n{output}\n")
    print(output)

    total_time = pretty_time_delta(int(time.time() - t0))
    logging.info(f"The program took {total_time}")


if __name__ == '__main__':
    main()
