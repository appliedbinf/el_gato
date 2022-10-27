#!/usr/bin/env python3

"""
Logic:

If one allele:
	blast allele

If 2 alleles:
	resolve two sequences
		if 1 SNP = easy
			create 2 seqs that differ at that site
		if multiple SNPS = hard
			assess phasing of SNPs
			are the SNPs on any of the same reads or pairs?
			create 2 seqs with SNPs incorporated
	blast alleles
	one or both have normal gene context
"""

import sys
import re
from collections import Counter

class SAM_data(object):
	"""stores columns of SAM entry as attributes"""
	def __init__(self, object):
		self.qname = object.split('\t')[0]
		self.flag = object.split('\t')[1]
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




insam = "reads_vs_momps_ref_filt_4.sam"
contig_name = "mompS_ref"
start_pos = 367
end_pos = 718

contig_dict ={}

with open(insam, 'r') as sam:
	for line in sam.readlines():
		if line[0] != "@":
			entry = SAM_data(line)
			if entry.rname in contig_dict.keys():
				contig_dict[entry.rname].append(entry)
			else:
				contig_dict[entry.rname] = [entry]

# reads_dict {bases: {pos:base}, qualities: {pos:qual}, readnames: {pos:rname}}
reads_dict ={
	'bases' : {k:[] for k in range(start_pos-1, end_pos)},
	'qualities' : {k:[] for k in range(start_pos-1, end_pos)},
	'readnames' : {k:[] for k in range(start_pos-1, end_pos)}
}


momps_reads = contig_dict[contig_name]
x = 0
for read in momps_reads:
	if ( 
		(read.pos <= start_pos and start_pos <= read.end) or
		(read.pos <= end_pos and end_pos <= read.end) or
		(start_pos <= read.pos and read.end <= end_pos) or
		(read.pos <= start_pos and end_pos <= read.end)
		):
		if read.pos < start_pos:
			pad = start_pos - read.pos
		else:
			pad = 0
		

		for read_idx in range(pad, min([read.ln - 1, 719-read.pos])):
			ref_seq_idx = read.pos + read_idx - 1
			reads_dict['bases'][ref_seq_idx].append(read.seq[read_idx])
			reads_dict['qualities'][ref_seq_idx].append(read.qual[read_idx])
			reads_dict['readnames'][ref_seq_idx].append(read.qname)

seq = []

for position in range(start_pos-1, end_pos):
	count = Counter(reads_dict['bases'][position])
	if len(count) == 1:
		seq.append([[b for b in count][0]])
	else:
		total = sum([i for i in count.values()])
		for base, num in count.items():
			if num > 0.9*total: ########## THRESHOLD #################
				seq.append([base])
			elif num > 0.1*total:
				seq.append([b for b in count.keys()])
				break
 
num_alleles_per_site = [len(i) for i in seq]

n_multiallelic = len([i for i in num_alleles_per_site if i > 1])

if n_multiallelic == 1:
	alleles = ['']*max(num_alleles_per_site)
	for base in seq:
		for i in range(max(num_alleles_per_site)):
			if len(base) == 1:
				alleles[i] += base[0]
			else:
				alleles[i] += base[i]
print("\n".join(alleles))


