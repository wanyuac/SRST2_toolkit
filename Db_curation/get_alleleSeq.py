"""
This script retrieves allele sequences from the ARG-Annot database according to a list of allele names.

Inputs:
	-db: the ARG-ANNOT database, which is a FASTA file
	-l: a text file containing a list of allele names, each line contains a single allele name.
		You can generate this list in R by using write.table(...).

Example:
	python get_alleleSeq.py -db <ARG-Annot database name> -l <a list of alleles> > seq.fna

Author: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Date of development: 9 July 2015

License: BSD
"""

from Bio import SeqIO
from argparse import ArgumentParser
import sys

def parse_arguments():
	# read arguments from the command line
	parser = ArgumentParser(description= "Read arguments: -db and -l")
	parser.add_argument("-db", type = str, required = True, help="The ARG-ANNOT database", default="")
	parser.add_argument("-l", type = str, required = True, help="The minimum % nucleotide identity to count as a hit")
	return parser.parse_args()

def get_field(descr, i):
	"""
	Extract a field from the header line (sequence description) of a FASTA record.
	Generally, the header should have the format: '>[cluster ID]__[gene name]__[allele name]__[sequence ID]'
	For example, >0__OqxBgb_Flq__OqxBgb__49
	"""
	fields = descr.split("__")
	return fields[i]  # return the allele name
	
def main():
	args = parse_arguments()
	db = SeqIO.parse(args.db, "fasta")  # read a multi-FASTA file and keep each sequence as an instance of SeqRecord class
	f = open(args.l, "rU")  # open the allele list
	allele_list = []
	for line in f:  # for each allele in the list
		allele_list.append(line.rstrip("\n"))  # remove the newline character and add the allele name to the list
	f.close()  # finish reading the allele list
	for seq in db:  # go through all records, searching for matches to the allele list
		allele = get_field(seq.description, 2)  # get the allele name from the sequence description
		if allele in allele_list:
			seq.id = allele  # override the sequence ID
			seq.description = ""  # remove the description
			SeqIO.write([seq], sys.stdout, "fasta")  # print this sequence to STDOUT
			"""
			For example: seq.id = "289__Aac3-Iva_AGly__Aac3-Iva__1489"
			seq.description = "289__Aac3-Iva_AGly__Aac3-Iva__1489 no;no;Aac3-Iva;AGly;X01385;244-1029;786"
			"""
	return

if __name__ == "__main__":
	main()