"""
This script finds discrepancies in SRST2's ARG-Annot database.
Usage: python evaluate_argannot_seqlen.py db.fasta > len.txt
Author: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Date: 11 July 2015
License: GNU GPL 2.1
"""

from Bio import SeqIO
import sys

print "seqid\treal_len\tdescr_len\tcalc_len"
for rec in SeqIO.parse(sys.argv[1], "fasta"):
	descr = rec.description
	list = descr.split(";")  # parse the description of the record
	descr_len = int(list[-1])  # take the last component
	coordinates = list[-2].split("-")  # parse the expression: "start-end"
	calc_len = abs(int(coordinates[1]) - int(coordinates[0])) + 1  # calculate the sequence length according to the described coordinates
	real_len = len(rec.seq)  # measure the actual sequence length
	if (real_len != descr_len) or (calc_len != descr_len):
		print "\t".join([rec.id, str(real_len), str(descr_len), str(calc_len)])  # print the result for manual investigation