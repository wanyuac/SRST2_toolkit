"""
This script removes redundant sequences and fixes genomic coordinates in SRST2's ARG-Annot database.

Inputs:
	-in: the ARG-ANNOT database, which is a FASTA file
	-remove: a text file containing a list of sequence identities whose sequences will be removed from the database
	-coord: whether to fix problems in genomic coordinates. Its value becomes True if this option is chosen.
	-out: the name of the output file

Outputs:
	<your_path>/
		<db name>.fasta: a new FASTA file. It is just a replica of the input database if -remove == "" and -coord is absent.
		ref.nsq...: files of an BLAST database named "ref"; delete them if you are going to add more sequences to the BLAST database.
		blast.txt: raw megaBLAST results in format 6
		/seq
			<accession numbers>.fna: nucleotide sequences downloaded from NCBI
			ref.fasta: the concatenated FASTA file for building a BLAST database
			accession_numbers.txt: a list of accession numbers of every sequence in the ARGannot database, may contain duplicates

Options of arguments:
	-remove
		"": skip this function
		[path/basename]: sequences on the list will be removed
	-coord
		all: sequence coordinates of all sequences in the input database will be examined and corrected
		[path/basename]: only sequences on the list will be investigated
		none: skip this function

Usage:
	python fix_ARGAnnotDB.py -in ARGannot_v130715.fasta -remove seqlist1.txt --coord seqlist2.txt -out ./ARGannot.fasta
	
Requirement: BLAST (both the local and internet based versions) must be accessible. For instance, to run "module load blast+-intel/2.2.29"
first on a Linux-based computational cluster.
Python version: 2.7.5

Copyright 2017 Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
Development history: 13-23 Jul 2015, 4-5 Sep 2017
"""

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from argparse import ArgumentParser
import sys, os, subprocess, tempfile

def parse_arguments():
	# read arguments of options
	parser = ArgumentParser(description="Fix problems in SRST2's ARG-Annot database")
	parser.add_argument("--in", "-i", dest="input", type=str, required=True, help="The ARG-ANNOT database to be fixed", default="")
	parser.add_argument("--remove", "-r", dest="remove", type=str, required=False, default="", help="A text file containing a list of sequence IDs whose sequences are to be removed")
	parser.add_argument("--coord", "-c", dest="coord", type=str, required=False, default="none", help="A text file containing IDs of sequences to be examined in terms of genomic coordinates")
	parser.add_argument("--out", "-o", dest="out", type=str, required=False, default="./ARGannot_fixed.fasta", help="File name of the output database, which includes the path")
	parser.add_argument("--email", "-e", dest="email", type=str, required=True, help="Your email address for retrieving records from NCBI")
	parser.add_argument("--downloader", "-d", dest = "downloader", type = str, required = False, default = "BINF_toolkit/download_NCBI_records.py",
						help = "Path to download_NCBI_records.py, which is used for downloading reference sequences")
	return parser.parse_args()

def read_entries(file):
	# read a text file and store each line as a component of a list
	with open(file, "rU") as f:
		lines = f.read().splitlines()  # read every line into a list and strip newline characters
		# The last line produced by the last allele name ended by a newline character, which is always empty, will be removed as well.
		
	return lines

def read_seqIDs(database):
	# returns identities of all sequences in the database
	ids = []
	for rec in database:
		ids.append(rec.id)
		
	return ids

def remove_redundancy(db, targets):
	print "Removing %d redundant sequences..." % len(targets)
	n = 0
	res = []  # the result list
	for rec in db:
		if rec.id in targets:
			print "Sequence %s is removed successfully." % rec.id
			n += 1
			continue
		res.append(rec)	
	print str(n) + " sequences have been removed."
	
	return res
	
def make_blastDB(db, selected, targets, outdir, email, downloader):
	# prepare the BLAST database
	print "Preparing the BLAST database containing %d sequences." % len(targets)
	path = os.path.join(outdir, "seq")  # in Linux: outdir/seq		
	if os.path.exists(outdir + "/ref.nsq"):
		print "BLAST database has been created already."
	else:
		acc_file = outdir + "/accession_numbers.txt"
		f = open(acc_file, "w")
		for rec in db:
			if selected:  # if only assess a given subset of sequences
				if not (rec.id in targets):
					continue  # skip this record
			# for all other cases
			f.write(rec.description.split(";")[-3] + "\n")  # get the accession number of its reference sequence and write it to a buffer
		f.close()  # write the buffer to the temporary file
		cmd = "python %s -r %s -f fasta -e %s -x fna -o %s -sk" % (downloader, "file:" + acc_file, email, path)
		print "Executing the command: " + cmd
		os.system(cmd)  # The downloader will make a new directory "seq" if it is absent.
		print "Concatenating FASTA files into a single one."
		os.system("cat " + path + "/*.fna > " + path + "/ref.fasta")  # concatenate sequences
		print "Making a local BLAST database."
		os.system("makeblastdb -in %s -dbtype nucl -out %s" % (path + "/ref.fasta", outdir + "/ref"))  # create a database named "ref"
		
	return

def run_blast(query, database):
	cmd = ["blastn", "-query", query, "-db", database, "-max_target_seqs", "1", "-max_hsps", "1", "-perc_identity", "100",
			"-task", "megablast", "-evalue", "0.001", "-outfmt", "6 qseqid sseqid sstart send qlen length gaps"]
	print " ".join(cmd)
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out = p.communicate()  # obtain the output of BLAST
	
	return out[0]  # stderr: out[1]

def parse_blastResult(val):
	print "Parsing the BLAST result line: " + val
	if val != "":
		fields = val.split("\t")[1:]  # [sseqid, sstart, send, qlen, length, gaps], drop the qseqid column
		#id = fields[0].split("|")[3][:-2]  # eg. only take AY649755 from the original sseqid "gi|532235082|gb|AY649755.2|" (NCBI does not use this format anymore)
		id = fields[0].split(".")[0]  # eg. EU370913.1 => EU370913
		dict = {"sseqid" : id, "sstart" : int(fields[1]), "send" : int(fields[2]), "qlen" : int(fields[3]),
				"length" : int(fields[4]), "gaps" : int(fields[5])}
	else:
		dict = None
		
	return dict
	
def parse_seqDescr(descr):
	# parses a SeqRecord's description and returns a dictionary of metrics
	list = descr.split(";")[-3:]  # parse the description of the record, only take the last three components
	coordinates = list[1].split("-")  # parse the expression: "start-end"
	dict = {"accession" : list[0], "sstart" : int(coordinates[0]), "send" : int(coordinates[1]), "length" : int(list[2])}
	
	return dict

def check_refError(nominal, actual):
	# check errors in coordinates and the accession number of the reference sequence
	return nominal["sstart"] != actual["sstart"] or nominal["send"] != actual["send"] or nominal["length"] != actual["length"] or nominal["accession"] != actual["sseqid"]

def override_seqDescr(original_descr, new_coord):
	fields = original_descr.split(";")
	fields[-2:] = [str(new_coord["sstart"]) + "-" + str(new_coord["send"]), str(new_coord["length"])]  # override the last two elements
	fields[-3] = new_coord["sseqid"]
	
	return ";".join(fields)
	
def search_NCBI(qseq, qseqid, dir):
	# search NCBI database using a nucleotide sequence
	path = dir + "/remote_BLAST"
	xml_filename = path + "/" + qseqid + ".xml"
	if not os.path.exists(path):
		os.system("mkdir " + path)
	if not os.path.exists(xml_filename):
		print "Searching " + qseqid + " in NCBI database through BLASTn..."
		handle = NCBIWWW.qblast(program="blastn", database="nt", sequence=qseq, alignments=1, expect=0.001, format_type="XML", hitlist_size=1, megablast=True, perc_ident=100)  # search in prokaryotae database
		xml = open(xml_filename, "w")
		xml.write(handle.read())  # The result can only be read once.
		xml.close()
		handle.close()
	# parse the result
	xml = open(xml_filename, "rU")
	blast_rec = NCBIXML.read(xml)  # for results generated by a single query sequence, returns a BLAST record class: type(blast_rec) = Bio.Blast.Record.Blast.
	if len(blast_rec.alignments) > 0:  # if there is at least one alignment in the result. You can also use: len(blast_rec.descriptions[0]) > 0
		hsp = blast_rec.alignments[0].hsps[0]  # I only take the first HSP in the first alignment, ignoring others.
		querylen = len(qseq)  # length of the query sequence
		matchlen = len(hsp.match)  # length of matched nucleotides
		"""
		Define an exact match:
			query coverage = 100% <=> querylen == matchlen
			identity of matched range = 100%, which includes situations in which gaps = 0. I already ensured perc_ident = 100% in megaBLAST search.
		"""
		if querylen == matchlen:  # 100% query coverage and 100% alignment identity
			info = blast_rec.alignments[0].title.split("|")[3:]  # an example of titles: gi|1245353|gb|U50278.1|ECU50278 Enterobacter cloacae class A carbepenem-hydrolyzing...
			dict = {"sseqid" : info[0][:-2], "sstart" : hsp.sbjct_start, "send" : hsp.sbjct_end, "sbjct_descr" : info[1], "length" : querylen}  # info[0][:-2]: drop the version number from the accession number
		else:  # query coverage < 100%
			dict = None  # no exact match as well. For example, the query sequence may has a mutation at its beginning or end regions, but identical in the matched region.
			"""
			In this case, I am not going to return the subject alignment as a replacement to the original sequence, because I think the user should check the problem manually if there is no exact hit in the whole NCBI nucleotide database,
			which could be a problem larger than just assigning a wrong accession number to a known sequence.
			"""
	else:  # no exact match returned from the NCBI server
		dict = None
	return dict

def fix_coordinates(db, selected, targets, outdir, email, downloader):
	print "Checking errors in genomic coordinates and fix them for %d sequences" % len(targets)
	corrected = 0  # number of corrected sequences
	failed = 0
	newDB = []
	make_blastDB(db, selected, targets, outdir, email, downloader)  # prepare the BLAST database
	f = open(outdir + "/blast.txt", "w")  # save results from megaBLAST
	f.write("qseqid\tsseqid\tsstart\tsend\tqlen\tlength\tgaps\n")  # write column names
	for rec in db:
		if selected:  # selected == True: only assess a given subset of sequences
			if not (rec.id in targets):
				newDB.append(rec)
				continue  # go to the next record
		# for others, recalculate genomic coordinates
		tmp = tempfile.NamedTemporaryFile(mode="w+t", delete=True)
		SeqIO.write(rec, tmp, "fasta")
		tmp.flush()
		blast_line = run_blast(tmp.name, outdir + "/ref")  # Get a single line from the BLAST output.
		f.write(blast_line)  # save raw BLAST result to a file, which is useful for debugging and validation
		blast = parse_blastResult(blast_line)  # blast becomes a dictionary variable if an alignment is found.
		if blast == None:  # no match is found in the reference
			remote_rec = search_NCBI(rec.seq, rec.id, outdir)  # try to find a perfect match in NCBI's database
			if remote_rec == None:
				print "Error: " + rec.id + " does not match to either its reference or any sequence in the NCBI database. Thus no correction was applied."
				newDB.append(rec)  # no correction is applied to this record
				failed += 1
			else:
				rec.description = override_seqDescr(original_descr=rec.description, new_coord=remote_rec)
				newDB.append(rec)
				corrected += 1
				print "Warning: " + rec.id + " does not match to its reference, but an exact match was found in NCBI database. A new sequence description was generated based on this record."
				print "\t>>New accession number: %s | %s" % (remote_rec["sseqid"], remote_rec["sbjct_descr"])
		else:
			slen_aligned = abs(blast["send"] - blast["sstart"]) + 1  # length of the aligned region on the subject sequence
			"""
			To define an exact match, as I said before, we only need to specify qlen == length. The identity percentage has already been guaranteed by perc_ident=100 in my command for BLAST.
			"""
			if blast["qlen"] == blast["length"]:  # This is a perfect match. So we check its coordinate and make a correction.
				nominal_coord = parse_seqDescr(rec.description)  # get nominal genomic coordinates
				if check_refError(nominal_coord, blast):  # if the position or accession number is wrong
					rec.description = override_seqDescr(original_descr=rec.description, new_coord=blast)
					newDB.append(rec)  # add the corrected record into the result
					corrected += 1
					print "Warning: " + rec.id + " has wrong/outdated reference information. Now it has been corrected."
				else:  # no error, skip
					#print rec.id + " has no error in its genomic coordinate or nominal length. Skip this sequence."
					newDB.append(rec)  # no correction is needed.
			else:  # query coverage < 100%
				remote_rec = search_NCBI(rec.seq, rec.id, outdir)  # try to search for a better match
				if remote_rec == None:
					print "Error: " + rec.id + " does not completely match to either its reference or any sequence in the NCBI database. No correction was applied."
					newDB.append(rec)  # do nothing to this record, just copy
					failed += 1
				else:  # if there is a perfect match found on NCBI's database
					rec.description = override_seqDescr(original_descr=rec.description, new_coord=remote_rec)
					newDB.append(rec)
					corrected += 1
					print "Warning: " + rec.id + "does not completely match to its reference, but an exact match was found in NCBI database. A new sequence description was generated based on this record."
					print "\t>>New accession number: %s | %s" % (remote_rec["sseqid"], remote_rec["sbjct_descr"])
		tmp.close()  # delete the temporary file
	f.close()
	print "Coordinates or lengths of %d sequences have been corrected." % corrected
	print "%d sequences failed correction." % failed
	return newDB

def write_multifasta(data, file):
	print "Writing sequences into " + file + "..."
	n = 0
	with open(file, "w") as f:
		for rec in data:
			SeqIO.write(rec, f, "fasta")
			n += 1
	print str(n) + " sequences has been written to the new database."
	return

def main():
	args = parse_arguments()
	db = list(SeqIO.parse(args.input, "fasta"))  # import the whole database as a list of SeqRecord instances
	(path, basename) = os.path.split(args.out)  # parse the file name
	if not os.path.exists(path):
		os.system("mkdir " + path)
	# Phase 1: remove sequence redundancy
	if args.remove != "":  # if the user want to remove redundant alleles
		db = remove_redundancy(db, read_entries(args.remove))
	# Phase 2: check and correct genomic coordinates
	if args.coord == "none":
		pass
	elif args.coord == "all":
		db = fix_coordinates(db=db, selected=False, targets=read_seqIDs(db), outdir=path, email=args.email, downloader = args.downloader)
	else:
		db = fix_coordinates(db=db, selected=True, targets=read_entries(args.coord), outdir=path, email=args.email, downloader = args.downloader)
	# Finally, write the processed database to a new FASTA file
	write_multifasta(db, args.out)
	print "Done."

if __name__ == "__main__":
	main()