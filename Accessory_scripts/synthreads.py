#!/usr/bin/env python

import sys
import os
import argparse
import random
import pathlib


GENE_UPDOWN = {
	'flaA':
		{'up': 'TAAGAATTAACAGTGCTAAAGATGATGCAGCAGGACTCGCCATTTCTCAACGAATGACTGCACAAATTCGCGGGATGAATCAAGCAGTTAGAAATGCTAACGATGGTATTTCTCTGGCTCAGGTTGCTGAAGGAGCTATGCAGGAAACAACCAATATTTTGCAGCGTATGCGTGAGCTTTCCGTTCAAGCTGCTAACTCAACCAATAACTCATCTGACCGGGCTTCAATTCAAAGCGAAATTTCACAATTAAAAAGTGAATTGGAGCGTATTGCTCAAAATACTGAATTCAATGGTCAAAGAATTCTCGATGGCTCTTTCTCTGGTGCAAGCTTCCAGGTTGGAGCGAAT',
		'down': 'AGATGCTACGTCTGCCTATGCCAAAGCAGCGGCTATCAACGATGCCGGGATAGGAGGACTATCAGTTACAGCATCTACCAGCGGCACACAAGCAGTTGGCGCAATAGGTGGAACCGCAGGTGATACCTACAACTTAACGATTAATGGTGTTGCAATATATACTAACCTTAACGTAGCAACGGCATTAACCAACTCCGATCTTCGGGATGCCATCAATGGTGTCAGTAACCAAACTGGAGTTGTCGCGTCACTTAACGGAGGTAATATGACACTCACTGCTGCTGATGGCAGAAACATTACGGTAACAGAAAGCGGTACAGGTTTTACTGCAGGTACTGACGGTTTAACT'},
	'pilE':
		{'up': 'GGAGCATTTATTGCAACGCTAAACACTCCTCCACTCAATTCGTGCGGTTTTGGTTTCACATCCATGCTTCTGGAGCTTAATTTCTTAAATGGCGGTGCTTTTCAATATGCGCGATTTGATATTTCCGGCGATGGAGGCTTTGATATAGCCGATCAGTATAATGGGGCCTACCCCGTTGGAATTGGCTTGTCCAACAGTTATGCCAACTCACCAACTGTTCTGGGGCCTAACAAGGATAACAATTTGGTCATTCTGATCACACAATCGGATGGAACACAAACTACTGTTATCGCTCCTAACATTACGCCAAGAAAAATAGGCTGGTGGCAAATTCAATAGAGGATACCCAA',
		'down': 'GTTGGAACTCAAGCCAAAGATACCGAGTGTGCCAGCATGTCAATTAACCAGGCCAATGTAAAAACAGCAGTAGATTCCTCCGCTAATGCGCAACCAGAATGCTGGAATCCCGGTTAATACTGTAAAAAATTCTTTACACAAAGGATTCAATTCGAGTAGATTTGCCAGTAAGGGAAAGTCCAACTTAATTTGAAATGGAAATTTTTTGCTTTTTAATAGGCGTACTTTATGCCTATACCTATAATATTTATCTAACTATTGCGATTCCTTTACTCTTTTATATCACACCAAGATATCCCATTGTCCTGTTCTTCATTGCCGGGTGTGGCTTTGCCCTGATTCATCAGTG'},
	'asd':
		{'up': 'CCAATTAAGTCCCTTTACCCATTAGCGAGTTCTCGCTCAGTAGGTAAGACAGTGACTTTTAGAGATCAGGAATTGGATGTTTTGGATTTAGCAGAATTTGATTTTAGCAAAGTCGATCTTGCTTTATTCTCTGCAGGCGGCGCTGTTTCAAAAGAGTATGCGCCAAAAGCTGTGGCTGCTGGGTGCGTGGTGGTTGATAATACCTCCTGTTTTCGTTATGAGGATGATATTCCTCTTGTTGTCCCTGAAGTGAATCCTCATCGAATCGCTGATTATACTAAAAGAGGCATTATTGCTAACCCTAATTGCTCCACCATTCAGATGGTCGTTGCTCTAAAGCCAATCTATGA',
		'down': 'ATATTTCTCATCCTTGTGGACTAAATTTGTGGATAGTCGCGGATAACATTCGCAAAGGTGCAGCAACTAACGCAGTACAAATTGCTGAAATACTGCAAAAAGAGTTTCTGTGAAGCTCCCACCCCCGCAATAGCATAGGTCCAGTATCTTGATAGAAAACCTGTATGCCCCCTAAGGTCATGCAGGTTATATTAGAAAAAATGCCATACTTTTTATATCCTAAGTTATACTTCTATTATGACTGTTAAATTTTGAAAAAGATGGCGAAAGAATATTATTTATCCGATGAGATCGAAGAAATAGTGCTCGCGAGTACTGTATTGTCTTACTCTCGACCCAATTTATCATC'},
	'mip':
		{'up': 'GAAAATAAATCAAAAACACTAATGTTCATCGCCGTTAAAATCTCTTGTTCATTATTAGGGGCAAGTGTAGAAGGATATTACCTTTTTGTCCATTATATAATTAAATGATAGCTTATGACTGGTAATTTTACGCAAATTTAGGCAGAATTAGTGGGCGATTTGTTTTTGCTTTATTTTGCTCAATTTATTGTGCAGTATGAGAACTTAAGTGTAAGACTAAAAGGGGATTGTTTATGAAGATGAAATTGGTGACTGCAGCTGTTATGGGGCTTGCAATGTCAACAGCAATGGCTGCAACCGATGCCACATCATTAGCTACAGACAAGGATAAGTTGTCTTATAGCATTGG',
		'down': 'CCAGGTTTCACAAGTTATCCCTGGATGGACAGAAGCTTTGCAATTGATGCCAGCTGGATCAACTTGGGAAATTTATGTTCCCTCAGGTCTTGCATATGGCCCACGTAGCGTTGGCGGACCTATTGGCCCAAATGAAACTTTAATATTTAAAATTCACTTAATTTCAGTGAAAAAATCATCTTAAGTTTTTTTGAATTAAAGTCATACAAAACGCATCCTTCTCATTTAGAGAGGGATGCTCTCTTTGTAAAGGCTAATGATCTTCATAAAAGGTGTCAGCCTAACACCACTGAGGAAATTAAAATATGTCTGTTTTAAAGGCTATGACTCATACTGAGTGAACAAAGGAA'},
	'mompS':
		{'up': 'GTTATCAATAAAATGGAAACTCAATAATAAACAAGTGGAGACAAGGCATGTTTAGTTTGAAAAAAACAGCAGTGGCAGTACTCGCCTTAGGAAGCGGTGCAGTGTTTGCTGGAACCATGGGACCAGTTTGCACCCCAGGTAATGTAACTGTTCCTTGCGAAAGAACTGCATGGGATATTGGTATCACCGCACTATATTTGCAACCAATCTATGATGCTGATTGGGGCTACAATGGTTTCACCCAAGTTGGTGGCTGGCAGCATTGGCATGATGTTGACCATGAGTGGGATTGGGGCTTCAAATTAGAAGGTTCTTATCACTTCAATACTGGTAATGACATCAATGTGAACTGGTATCATTTTGATA',
		'down': 'TATCGTTCCTGAGTTGGAAGCTAAGCTTGGTGCTGATTACACTTACGCAATGGCTCAAGGCGATTTGACTTTAGACGTTGGTTACATGTGGTTTAACTACTTCAACGCTATGCACAATACTGGCGTATTTAATGGATTTGAAACTGATTTCGCAGCTTCTGGTCCTTACATTGGCTTGAAGTATGTTGGTAATGTGTAATTTGTTAAGTTGATAAGAAATTTCAGCAATACTGTTGACTTTATAGAAGTCCGGCTGGATAATTTATCCA'},
	'proA':
		{'up': 'GGACTTGCCATTATTAGAAATCACACGTGATTCAAGTGTTGAAATGTGCTTTATGGAAAATACCGATGTTAAAGTAGTGGACATGGGACATAAGTACTATTCAAATAACAAACCTATGCAATTTACTTGTAAAGAGACTCCCGATACTCAATCTACTAAAACTTATTATACCGGATATTCCGCTGATGGTTATGATAGAGATAATGGAGCCGCTTCTCCTACCAATGACGCATTGTATGCTGGTTATGTCATTAAACACATGTACCACGATTGGTATGGTGTAGAGGCTTTAACTAAATCAGATGGATCGCCAATGCAATTAGTTATGCGAGTGCATTATGGTCAAGGC',
		'down': 'TGGAATCTTCGTATGGCTTTTGATGTTATGGTGAAGGCCAATATGGATTATTGGACACCTTATTCAACATTTGATGAGGGTGGTTGCGGTATGTTGAGTGCTGCCAAAGATTTGGGTTACAATTTGGATGATATTAAGAAGTCATTGAGTGAAGTGACTATAAATTATCAATCTTGTTATGTCGATTAATTCATTTTATTGATTGATAAAATTCAGCCAGAAGTCGACTGATGACTTCTGGCTTTTCCATATGAATATGGTGTCCACCATCTAACTTCTCAATTTTAATATTTTTTACTGCTTGAATCCTTGCTTTCATCTTATCCGAATCAAATGAAAATCCTTTGCTG'},
	'neuA':
		{'up': 'GGCTATTAGCGCAGTAATAAAAGATGAAAATATGTACGGATGCCAATTTCATCCTGAAAAAAGCGGGGAAGTAGGTTTAAGCATCATTCAACAGTTTTTGCAGATTTAGGGTGAATTAAAAATGAGAATATTGGCAGTAATCCCGGCAAGGGCTGGCTCAAAGAGGCTGCCCGGTAAAAATACCAGATTACTTGCCGGAAAGCCATTAATTGCACATACTATTGTTGCTGCCTTGCAGTCGTCTTGTTGTGAAGAAATCGTTGTTTCAACCGATAGTAAACAAATAGCAGACGTCGCCGTTCAATATGGGGCTTCAGTACCCTGGCTAAGATCGGAAGATTTAGCCACG',
		'down': 'GAACCAACCAAACCTTTATTGTTAAATAGTATTAGTGAATCCATCGACATCGATACGCCAATCGATTGGGCTCTAACAGAAAAATTAATGGAATTAAACCAAGAGGCTCTAGTATGACTTGTTTTATTATTGCTGAAGCAGGAGTAAACCATAATGGTGATCTTCAACTGGCTAAAGAATTAGTTTATGCAGCCAAAGAGTCAGGAGCTGATGCAGTAAAGTTTCAGACTTTCAAAGCGGACACCTTGGTTAATAAAACAGTAGAAAAAGCGGAATACCAAAAAAATAATGCCCCGGAATCCTCTACTCAGTATGAGATGCTCAAGGCACTGGAACTTTCAGAAGAGGAC'},
	'neuAH':
		{'up': 'CAGCCTATTAGTGCTGTTATAAACGATGGCAATATTTATGGATGTCAATTCCATCCAGAAAAAAGTGGAGAAGTAGGATTAAAGATCATCCAGCAGTTTTTACAAATTTAGGACGCAAGAAGTGAGAGTTTTAGCGATAATACCAGCAAGAGCTGGATCCAAAAGGCTACCAGGAAAAAATATAAGATTACTTGCTGGTAAACCACTTATCGCTCATACTATAAATGCGGCTCTACAATCCAATTGCTGCGAACAAATCGTAGTTTCAACTGATTGTCGAGAAATAGCTGATATCGCAATTCAATATGGAGCCTCTGTGCCATGGCTAAGGCCACACTCCTTAGCTGAT',
		'down': 'AATCCAACAAAGGCGCTGCTTATGGATAGTCCAAGTGAATCTATAGATATAGATACTCCAATGGATTGGGCATTAGTAGAAAAATTAATTGAACTAAAAGAAGAGATACTATTATGAGTTGTTTTATAATTGCAGAGGCAGGGGTAAATCATAATGGTGATCTTCAATTGGCTAAAGAATTAATTTATGCAGCCAAAGAGTCAGGAGCTGATGCGGTAAAATTTCAGACTTTCAAAGCTGACTCCCTGGTCAATAAAACAGTAGAAAAAGCGGAATACCAAAAAAATAATGCTCCAGAATCGACTACTCAGTATGAGATGCTTAAAGCTCTGGAGCTTACAGAAGAGGAT'}
}

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

class St():
	def __init__(self, line):
		try:
			bits = line.split()
			self.st = bits[0]
			self.flaA = bits[1]
			self.pilE = bits[2]
			self.asd = bits[3]
			self.mip = bits[4]
			self.mompS = bits[5]
			self.proA = bits[6]
			self.neuA_neuAH = bits[7]
			self.flaA_seq = ""
			self.pilE_seq = ""
			self.asd_seq = ""
			self.mip_seq = ""
			self.mompS_seq = ""
			self.proA_seq = ""
			self.neuA_neuAH_seq = ""
		except:
			print(line)
			sys.exit()


class Read():
	def __init__(self, n):
		self.rname = f'read_{n}'
		self.seq = ''
		self.qual = ''


def cmdline_args():

	p = argparse.ArgumentParser(
		description=""
		)
	p.add_argument(
		"-s", "--num_sts", required = False, type=int, default=1,
		help="number of distinct STs for which reads should be generated (default: 1)"
		)
	p.add_argument(
		"-n", "--num_reads", required = False, type=int, default=10_000,
		help="number of reads to be generated (default: 10,000)"
		)
	p.add_argument(
		"--av_frag_len", required = False, type=int, default=500,
		help="Approximate average fragment length (default: 500)"
		)
	p.add_argument(
		"--av_read_len", required = False, type=int, default=240,
		help="Approximate average read length (default: 240)"
		)
	p.add_argument(
		"-p", "--db_path", required = False, type=str, default="db/",
		help="path to allele sequence database (default: db/)"
		)
	p.add_argument(
		"-o", "--outprefix", required = False, type=str, default="./",
		help="prefix/path to prepend to output files (default: ./)"
		)
	p.add_argument(
		"--seed", required = False, type=int, default=None,
		help="Set the seed to control randomness of read generation"
		)
	p.add_argument(
		"sts", nargs="*", 
		help="list of ST profiles desired"
		)

	return p.parse_args()


def fasta_to_dict(FASTA_file):
	"""Read a fasta file into a dict 

	Dict has headers (minus the > symbol) as keys and the associated 
	sequence as values.
	
	Args:
	  FASTA_file (str): 
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
		header = i.split("\n")[0]
		seq = "".join(i.split("\n")[1:])
		fasta_dict[header] = seq

	return fasta_dict


def rev_comp(string):
	"""Reverse complement a string of nucleotide sequence
	
	Args:
	  string (str):
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


def load_sts(args):
	sts = {}
	with open(f"{args.db_path}/lpneumophila.txt", 'r') as fin:
		for line in fin.readlines()[1:]:
			if len(line.strip()) == 0:
				continue
			s = St(line)
			sts[s.st] = s

	return sts


def choose_sts(args, sts):
	st_choices = []
	if len(args.sts) == 0:
		choices = []
		for i in range(args.num_sts):
			choice = random.choice([k for k in sts.keys()])
			choices.append(choice)
			st_choices.append(sts[choice])
		print(f"STs selected: {' '.join(choices)}")
		return st_choices


	for s in args.sts:
		st_choices.append(sts[s])
	return st_choices


def load_allele_seqs(args):
	allele_seqs = fasta_to_dict(f"{args.db_path}/all_loci.fasta")
	return allele_seqs


def fill_st_seq(st_choices, allele_seqs):
	for choice in st_choices:
		choice.flaA_seq = (
			GENE_UPDOWN['flaA']['up']
			+ allele_seqs[f'flaA_{choice.flaA}']
			+ GENE_UPDOWN['flaA']['down']
			)
		choice.pilE_seq = (
			GENE_UPDOWN['pilE']['up']
			+ allele_seqs[f'pilE_{choice.pilE}']
			+ GENE_UPDOWN['pilE']['down']
			)
		choice.asd_seq = (
			GENE_UPDOWN['asd']['up']
			+ allele_seqs[f'asd_{choice.asd}']
			+ GENE_UPDOWN['asd']['down']
			)
		choice.mip_seq = (
			GENE_UPDOWN['mip']['up']
			+ allele_seqs[f'mip_{choice.mip}']
			+ GENE_UPDOWN['mip']['down']
			)
		choice.mompS_seq = (
			GENE_UPDOWN['mompS']['up']
			+ allele_seqs[f'mompS_{choice.mompS}']
			+ GENE_UPDOWN['mompS']['down']
			)
		choice.proA_seq = (
			GENE_UPDOWN['proA']['up']
			+ allele_seqs[f'proA_{choice.proA}']
			+ GENE_UPDOWN['proA']['down']
			)
		if int(choice.neuA_neuAH) < 200:
			var = 'neuA'
		else:
			var = 'neuAH'
		choice.neuA_neuAH_seq = (
			GENE_UPDOWN[var]['up']
			+ allele_seqs[f'neuA_neuAH_{choice.neuA_neuAH}']
			+ GENE_UPDOWN[var]['down']
			)

	return st_choices


def rand_seq(n):
	seq = ''
	for i in range(n):
		seq += random.choice(['A', 'T', 'C', 'G'])

	return seq


def sequence_fragment(seq, av_frag_len, av_read_len):
	mid = random.randint(0,len(seq))
	fragment_len = random.randint(0.8*av_frag_len,1.2*av_frag_len)
	start = int(mid-fragment_len/2)
	end = int(mid+fragment_len/2)
	if start < -(av_read_len-50):
		r1 = rand_seq(av_read_len)
	else:
		r1 = seq[max([start, 0]): start+av_read_len]
	if end > len(seq)+(av_read_len-50):
		r2 = rand_seq(av_read_len-50)
	else:
		r2 = rev_comp(seq[end-av_read_len: min([len(seq), end])])

	# make sure reads are long enough or make them up to between 100 and av_read_len with random seq
	if len(r1) < 100:
		add_len = random.randint(100-len(r1), av_read_len-len(r1))
		r1 = "".join([rand_seq(add_len), r1])
	if len(r2) < 100:
		add_len = random.randint(100-len(r2), av_read_len-len(r2))
		r2 = "".join([rand_seq(add_len), r2])

	switch = random.choice([True, False])

	if switch:
		r1, r2 = r2, r1
	return r1, r2


def generate_reads(n, st_choices, av_frag_len, av_read_len):
	reads1 = []
	reads2 = []
	for i in range(1, n+1):
		st = random.choice(st_choices)
		read1 = Read(i)
		read2 = Read(i)
		seq = random.choice([st.flaA_seq, st.pilE_seq, st.asd_seq, st.mip_seq, st.mompS_seq, st.proA_seq, st.neuA_neuAH_seq])
		read1.seq, read2.seq = sequence_fragment(seq, av_frag_len, av_read_len)
		read1.qual = "".join([chr(random.randint(58, 69)) for _ in range(len(read1.seq))])
		read2.qual = "".join([chr(random.randint(58, 69)) for _ in range(len(read2.seq))])
		reads1.append(read1)
		reads2.append(read2)

	return reads1, reads2


def write_reads_file(reads, filename):
	contents = ""
	for read in reads:
		contents += f"@{read.rname}\n{read.seq}\n+\n{read.qual}\n"

	with open(filename, 'w') as fout:
		fout.write(contents)


def main(args):
	# make outdir if doesn't exist
	if args.outprefix[-1] == "/":
		if not os.path.exists(args.outprefix):
			pathlib.Path(args.outprefix).mkdir(parents=True, exist_ok=True)
	else:
		directory = "/".join(args.outprefix.split("/")[:-1])
		if not os.path.exists(directory):
			pathlib.Path(directory).mkdir(parents=True, exist_ok=True)

	# Add trailing slash to db path if missing
	if args.db_path[-1] != "/":
		args.db_path += '/'

	random.seed(args.seed)
	sts = load_sts(args)
	st_choices = choose_sts(args, sts)
	allele_seqs = load_allele_seqs(args)
	st_choices = fill_st_seq(st_choices, allele_seqs)
	reads1, reads2 = generate_reads(args.num_reads, st_choices, args.av_frag_len, args.av_read_len)
	write_reads_file(reads1, f"{args.outprefix}reads1.fastq")
	write_reads_file(reads2, f"{args.outprefix}reads2.fastq")

	
if __name__ == '__main__':
	args = cmdline_args()
	main(args)	
