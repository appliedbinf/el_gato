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
from collections import Counter, defaultdict

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





insam = sys.argv[1] #"reads_vs_momps_ref_filt_4.sam"
contig_name = "mompS_ref"
start_pos = 367
end_pos = 718

contig_dict = {}
read_info_dict = defaultdict(list) #{'readname': [SAM_data, SAM_data]}

with open(insam, 'r') as sam:
	for line in sam.readlines():
		if line[0] != "@":
			entry = SAM_data(line)
			read_info_dict[entry.qname].append(entry)
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
				seq.append([b for b,c in count.items() if c > 0.1*total])
				break
 
num_alleles_per_site = [len(i) for i in seq]

n_multiallelic = len([i for i in num_alleles_per_site if i > 1])

if n_multiallelic == 0:
	allele = "".join([b[0] for b in seq])
	print(allele)


elif n_multiallelic == 1:
	alleles = ['']*max(num_alleles_per_site)
	for base in seq:
		for i in range(max(num_alleles_per_site)):
			if len(base) == 1:
				alleles[i] += base[0]
			else:
				alleles[i] += base[i]
	print("\n".join(alleles))

elif n_multiallelic == 2:
	multi_allelic_idx = [n for n,i in enumerate(num_alleles_per_site) if i == 2]
	reads_at_a = reads_dict['readnames'][multi_allelic_idx[0]+start_pos-1]
	reads_at_b = reads_dict['readnames'][multi_allelic_idx[1]+start_pos-1]
	intersect = sorted([i for i in set(reads_at_a).intersection(set(reads_at_b))])
	conflicting_reads = []
	read_pair_base_calls = []
	for read_name in intersect:
		read_pair = read_info_dict[read_name]
		calls_list = []
		for mate in read_pair: # If read is supplementary alignment, ignore.
			if mate.flag > 2047:
				continue
			calls_list.append(mate.get_base_calls(
				[multi_allelic_idx[0]+start_pos-1, multi_allelic_idx[1]+start_pos-1]
				))

		agreeing_calls = []
		for k in range(2):
			calls = []
			
			for mate_calls in calls_list:
				if mate_calls[k] != '':
					calls.append(mate_calls[k])
			if len(set(calls)) > 1:
				conflicting_reads.append(calls_list)
				break
			else:
				agreeing_calls.append(list(set(calls))[0])

		if len(agreeing_calls) == 2:
			read_pair_base_calls.append("".join(agreeing_calls))


	if len(conflicting_reads) > 0.1 * len(reads_at_a):
		print("more than 10% of reads disagree with which variant bases"
			" are in the same gene")

	biallele_results_count = Counter(read_pair_base_calls)
	total_read_count = sum([v for v in biallele_results_count.values()])

	bialleles = [k for k,v in biallele_results_count.items() if v > 0.1*total_read_count]

	if len(bialleles) > 2:
		print(f"{len(bialleles)} well-supported mompS alleles identified"
		" and can't be resolved. Aborting.")
		sys.exit()

	alleles = ['']*max(num_alleles_per_site)
	biallic_count = 0
	for base in seq:
		for i in range(max(num_alleles_per_site)):
			if len(base) == 1:
				alleles[i] += base[0]
			else:
				alleles[i] += bialleles[i][biallic_count]
		if len(base) > 1:
			biallic_count+=1

	print("\n".join(alleles))



# M03199:169:000000000-B8KWC:1:1101:28160:15045
# [646, 706]
