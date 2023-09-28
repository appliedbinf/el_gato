#!/usr/bin/env python3

import sys
import json
import math
from dataclasses import dataclass
from datetime import datetime
from fpdf import FPDF
from ctypes import alignment

LOGO="""\
 /\\_/\\  
( o.o ) 
 > ^ <\
"""

summary_header = """\
El_gato utilizes either a genome assembly (.fasta) or Illumina paired-end reads (.fastq) to replicate \
Legionella pneumophila Sequence Based Typing (SBT). From the input, 7 loci (flaA, pilE, asd, mip, \
mompS, proA, and neuA/neuAh) are identified and compared to a database of subtypes. The sequence \
type provided for each input sample is based on the unique combination of the allelic identities of the 7 \
target loci. It is possible that multiple alleles for a specific locus were identified and unable to be resolved properly. \
This would lead to a ST that indicates that one or more loci contain multiple alleles.\
"""

reads_header = """\
The following sample was analyzed using the paired-end reads functionality. The tables below show the full \
MLST profile of the sample, the coverage data for each locus, and information regarding the primers used to \
identify the primary mompS allele. Highlighted rows illustrate the data that caused failure of allele identification.\
"""

assembly_header = """\
The following sample was analyzed using the assembly functionality. The tables below show the full \
MLST profile of the sample and the corresponding locus location information. Unless specified by the user, \
el_gato utilizes a default 30% (0.3) BLAST hit length threshold and a 95% (95.0) sequence identity threshold \
to identify multiple alleles on multiple contigs ('?'). Highlighted rows illustrate the data that caused \
failure of allele identification.\
"""

cdc_header = """\
National Center for Immunization and Respiratory \n
Diseases (NCIRD) \n
Division of Bacterial Diseases (DBD) / \n
Respiratory Diseases Branch (RDB) \n
Pneumonia Response and Surveillance Laboratory (PRSL)\
"""
abbrev_key = """\
Novel ST = The alleles for all 7 loci were identified, however their unique combination and corresponding ST has not been found in the database. \n
Novel ST* = One or more locus failed to amplify, which may indicate a novel allele. \n
MA? = One or more locus contain multiple alleles that could not be resolved, leading to an ambiguous ST. \n
MD- = One or more locus was unidentifiable, leading to an unidentifiable ST. \n
'-' = There is missing data for the locus and could not be identified. \n
'NAT' = Amplification of the locus was unsuccessful, possibly due to the presence of a novel allele type. \n
'?' = The locus contains multiple alleles that could not be resolved.\
"""

primer_footer = """\
The primary mompS allele is identified using the following criteria: \n
1. Only one sequence has associated reads with the correctly oriented primers. \n
2. One sequence has more than 3 times as many reads with the correctly oriented primer as the other. \n
3. One sequence has no associated reads with the primer in either orientation, but the other has reads with the primer only in the wrong direction. The sequence with no associated reads is considered the primary locus in this case.\
"""
@dataclass
class Report(FPDF):
	sample_id: str
	st: str
	flaA: str
	pilE: str
	asd: str
	mip: str
	mompS: str
	proA: str
	neuA_neuAH: str
	mode: str
	mode_specific: dict

	@classmethod
	def from_json(cls, json_data):
		sample_id = json_data["id"]
		st = json_data["mlst"]["st"]
		flaA = json_data["mlst"]["flaA"]
		pilE = json_data["mlst"]["pilE"]
		asd = json_data["mlst"]["asd"]
		mip = json_data["mlst"]["mip"]
		mompS = json_data["mlst"]["mompS"]
		proA = json_data["mlst"]["proA"]
		neuA_neuAH = json_data["mlst"]["neuA_neuAH"]
		mode = json_data["operation_mode"]
		mode_specific = json_data["mode_specific"]
		x = cls(
			sample_id,
			st,
			flaA,
			pilE,
			asd,
			mip,
			mompS,
			proA,
			neuA_neuAH,
			mode,
			mode_specific,
		)
		return x

	def list_mlst(self):
		return [
			self.sample_id,
			self.st,
			self.flaA,
			self.pilE,
			self.asd,
			self.mip,
			self.mompS,
			self.proA,
			self.neuA_neuAH
			]

	def sample_report(
		self,
		pdf,
		typeface='Courier',
		body_style='',
		body_size=11,
		head_style='B',
		head_size=16
		):

		pdf.add_page()
		pdf.set_font(typeface, head_style, head_size)

		if self.mode == "Assembly":
			pdf = self.assembly_report(pdf, typeface, body_style, body_size)
		elif self.mode == "Reads":
			pdf = self.reads_report(pdf, typeface, body_style, body_size)
		else:
			sys.exit(
				f"Unsupported operation mode identified for sample {self.sample_id}"
				)
		return pdf


	def reads_report(self, pdf, typeface, style, size):
		pdf.set_font(typeface, style, size)
		pdf.set_font('Courier', 'B', 10)
		pdf.cell(80)
		pdf.multi_cell(h=4,w=0,txt="Epidemiology of Legionella: Genome-based Typing (el_gato) Paired-End Reads Report", align="C")
		pdf.ln(10)
		pdf.set_font('Courier', '', 11)
		pdf.multi_cell(
			w=0,h=5,
			txt=reads_header,
			new_x="LMARGIN", new_y="NEXT"
			)
		pdf.ln(10)

		pdf = self.make_mlst_table(pdf, [self.list_mlst()])
		pdf.ln(10)

		pdf.set_font(style="BU")
		pdf.cell(
			w=0,h=10,
			txt=f"Locus Coverage Information",
			new_x="LMARGIN", new_y="NEXT", align="C"
			)

		pdf.set_font()
		pdf = self.read_coverage_table(pdf)

		pdf.add_page()
		pdf.ln(10)
		pdf.set_font(style="BU")
		pdf.cell(
			w=0,h=10,
			txt=f"mompS Primer Information",
			new_x="LMARGIN", new_y="NEXT", align="C"
			)
		pdf.set_font()
		pdf = self.mompS_primer_table(pdf)
		pdf.ln(0)

		pdf.multi_cell(
			w=0, h=3.5,
			txt=primer_footer,
			new_x="LMARGIN", new_y="NEXT"
		)

		return pdf


	def read_coverage_table(self, pdf):
		contents = [["Locus", "Proportion Covered", "Minimum Coverage"]]
		contents += [
			[
				k, v["Proportion_covered"], v["Minimum_coverage"]
			] for k, v in self.mode_specific["locus_coverage"].items()]
		col_widths = (50, 50, 50)
		alignment = ("CENTER", "CENTER", "CENTER")

		highlight_rows = set()
		for n, row in enumerate(contents[1:]):
			if (
				(float(row[1]) < 100 and row[0] != "neuA_neuAH")
				or (float(row[1]) < 99 and row[0] == "neuA_neuAH")
			):
				highlight_rows.add(n+1)
		pdf = self.make_table(
			pdf,
			contents,
			col_widths=col_widths,
			text_align=alignment,
			highlight_rows=highlight_rows
			)

		return pdf

	def mompS_primer_table(self, pdf):
		contents = [["Allele", "Reads Indicating Primary", "Reads Indicating Secondary"]]
		contents += self.mode_specific["mompS_primers"]
		col_widths = (50, 50, 50)
		alignment = ("CENTER", "CENTER", "CENTER")
		if self.mode_specific["mompS_primer_conclusion"] == "inconclusive":
			highlight_rows = {1,2}
		else:
			highlight_rows = set()

		pdf = self.make_table(
			pdf,
			contents,
			col_widths=col_widths,
			text_align=alignment,
			highlight_rows=highlight_rows
			)
		pdf.ln(4)

		return pdf

	def assembly_report(self, pdf, typeface, style, size):
		pdf.set_font(typeface, style, size)
		pdf.cell(80)
		pdf.set_font('Courier', 'B', 10)
		pdf.multi_cell(h=4,w=0, txt="Epidemiology of Legionella: Genome-based Typing (el_gato) Assembly Results", align="C")
		pdf.ln(10)
		pdf.set_font('Courier', '', 11)
		pdf.multi_cell(
			w=0,h=5,
			txt=assembly_header,
			new_x="LMARGIN", new_y="NEXT"
			)
		pdf.ln(10)
		pdf = self.make_mlst_table(pdf, [self.list_mlst()])
		pdf.ln(4)
		pdf.cell(
			w=0,h=2,
			txt="Length Identity and Sequence Identity: " + self.mode_specific['length_id'] +"; "+self.mode_specific['sequence_id'] + "%",
			new_x="LMARGIN", new_y="NEXT"
		)
		pdf.ln(10)

		pdf.set_font(style="BU")
		pdf.cell(
			w=0,h=10,
			txt=f"Locus Location Information",
			new_x="LMARGIN", new_y="NEXT", align="C"
			)
		pdf.set_font()
		pdf = self.locus_location_table(pdf)

		return pdf


	def locus_location_table(self, pdf):
		contents = [["locus", "allele", "contig", "start", "stop"]]
		highlight_rows = set()
		x = 1
		for k, v in self.mode_specific["BLAST_hit_locations"].items():
			contents.append([k] + v[0])
			if len(v) > 1:
				for _ in range(len(v)):
					highlight_rows.add(x)
					x+=1
				for row in v[1:]:
					contents.append([""] + row)
			else:
				x+=1
		col_widths = (25, 30, 70, 20, 20)
		alignment = ("CENTER", "CENTER", "CENTER", "CENTER", "CENTER")

		content = [i for i in contents]
		batches = self.fit_table(pdf, content, pdf.get_y(), 34)
		pdf = self.make_table(
			pdf,
			batches[0],
			col_widths=col_widths,
			text_align=alignment,
			highlight_rows=highlight_rows
			)
		if len(batches) > 1:
			for batch in batches[1:]:
				pdf.add_page()
				pdf.set_y(pdf.get_y() + 10)
				pdf = Report.make_table(
					pdf,
					batch,
					col_widths=col_widths,
					text_align=alignment,
					highlight_rows=highlight_rows
				)

		return pdf

	@staticmethod
	def make_table(pdf, data, col_widths=None, text_align=None, highlight_rows=set()):
		with pdf.table(
			col_widths=col_widths,
			text_align=text_align,
			#borders_layout="MINIMAL"
		) as table:
			for n, data_row in enumerate(data):
				row = table.row()
				if n in highlight_rows:
					pdf.set_fill_color(233, 79, 88)
				else:
					pdf.set_fill_color(0, 0, 0)
				for item in data_row:
					row.cell(item)
			return pdf
	
	@staticmethod
	def fit_table(pdf, data, initial_y, characters:int):
		font_size = pdf.font_size
		pdf_y = initial_y
		n = 0
		max_length = 0
		batches = []
		this_batch = []
		while n < len(data):
			row = data[n]
			for i in row:
				column_length = len(i)
				if column_length > max_length:
					max_length = column_length
			num_lines = math.ceil(max_length / characters)
			cell_height = 2* num_lines * font_size
			if pdf_y + cell_height > pdf.page_break_trigger:
				batches.append(this_batch)
				this_batch = [row]
				n+=1
				pdf_y = 60 + cell_height # Whatever we want the starting y position to be on a new page
				continue

			n+=1
			pdf_y += cell_height
			this_batch.append(row)
		batches.append(this_batch)
		return batches
		

	@staticmethod
	def make_mlst_table(pdf, data):
		contents = [["Sample ID","ST","flaA","pilE","asd","mip","mompS","proA","neuA"]]
		for sample in data:
			contents.append(sample)
		col_widths = (50, 20, 20, 20, 20, 20, 20, 20, 20)
		alignment = ("CENTER", "CENTER", "CENTER", "CENTER", "CENTER", "CENTER", "CENTER", "CENTER", "CENTER")
		pdf = Report.make_table(pdf, contents, col_widths=col_widths, text_align=alignment)
		return pdf

	@staticmethod
	def read_jsons(files):
		data = []
		for file in files:
			with open(file) as fin:
				json_data = json.load(fin)
				data.append(Report.from_json(json_data))
		return data
	
	@staticmethod
	def read_multi_json(files):
		data = []
		with open(files) as fin:
			json_data = json.load(fin)
			for i in json_data:
				data.append(Report.from_json(i))
				
		return data	

class PDF(FPDF):
	def header(self):
		self.image("https://www.drought.gov/sites/default/files/styles/i_square_240_240/public/hero/partners/CDC_logo.png.webp?itok=Yja6wdoM", 10, 8, 33, keep_aspect_ratio=True)
		self.set_font('Courier', '', 10)
		self.cell(80)
		self.multi_cell(h=2,w=0, txt=cdc_header, align="C")
		self.ln(2)

def main():
	with open(sys.argv[1]) as fin:
		if fin.read().startswith("["):
			data = Report.read_multi_json(sys.argv[1])
		else:
			data = Report.read_jsons(sys.argv[1:])
	pdf = PDF('P', 'mm', 'Letter')
	pdf.add_page()
	pdf.set_font('Courier', 'B', 10)
	pdf.cell(80)
	pdf.multi_cell(h=4,w=0,txt="Epidemiology of Legionella: Genome-based Typing (el_gato) Batch Results Report", align="C")
	pdf.ln(2)
	pdf.set_font('Courier', '', 16)
	pdf.multi_cell(w=0,h=6, txt=LOGO, new_x="LMARGIN", new_y="NEXT")
	pdf.ln(4)
	pdf.set_font('Courier', '', 11)
	pdf.multi_cell(w=0,h=5, txt=summary_header, new_x="LMARGIN", new_y="NEXT")
	pdf.ln(10)
	content = [i.list_mlst() for i in data]
	batches = Report.fit_table(pdf, content, pdf.get_y(), 19)
	for batch in batches:
		if batch != batches[-1]:
			pdf.set_font('Courier', '', 11)
			pdf = Report.make_mlst_table(pdf, batch)
			pdf.add_page()
			pdf.ln(10)
		else:
			pdf.set_font('Courier', '', 11)
			pdf = Report.make_mlst_table(pdf, batch)
			pdf.ln(5)
	pdf.set_font(style="U")
	pdf.cell(w=0,h=0, txt="Abbreviation Key", new_x="LMARGIN", new_y="NEXT")
	pdf.ln(5)
	pdf.set_font()
	pdf.multi_cell(w=0,h=3.5, txt=abbrev_key, new_x="LMARGIN", new_y="NEXT")

	
	for datum in data:
		pdf = datum.sample_report(pdf)

	pdf.output("report.pdf")

if __name__ == '__main__':
	main()
