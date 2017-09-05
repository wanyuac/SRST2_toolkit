"""
This script splits a FASTA file that contains consensus sequences of both MLST and other genes into two files, each only contains either
MLST alleles or other sorts of alleles. It is useful if you want to separate different kinds of sequences from SRST2 results when SRST2
were run for MLST and genotyping simulatenously.

Usage
    python split_mlst_from_others.py -i *_samples.all_consensus_alleles.fasta -om ~/mlst -og ~/genes

Options
    -i: input FASTA files
    -om, -og: output directories for MLST genes and other genes, respectively. Importantly, om must not equal og.

Author: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Development history: 5 Sep 2016
Python version: 3.4.3
License: GNU GPL 2.1
"""

from argparse import ArgumentParser
import sys, os

def parse_arguments():
    # read arguments of options
    parser = ArgumentParser(description = "Separate two kinds of consensus sequences")
    parser.add_argument("-i", nargs = '+', type = str, required = True, default = "", help = "A list of input FASTA files")
    parser.add_argument("-om", type = str, required = False, default = "mlst", help = "reference nucleotide database for BLAST")
    parser.add_argument("-og", type = str, required = False, default = "others", help = "(optional) Comma-delimited names of bacterial strains")
    return parser.parse_args()
 
def main():
    args = parse_arguments()
    if args.om == args.og:
        sys.exit("Output directories must be different for MLST and other genes.")
    
    # make new folders at any places where possible
    if not os.path.exists(args.om):
        os.mkdir(args.om)
    if not os.path.exists(args.og):
        os.mkdir(args.og)
    
    # loop through every FASTA file
    nonMLST = False  # a flag marking whether the current sequence is from a non-MLST gene
    for f in args.i:
        fasta = open(f, "rU")
        o_mlst = open(os.path.join(args.om, f), "w")
        o_gene = open(os.path.join(args.og, f), "w")
        for line in fasta:  # keep newline characters at the end of each line
            if line.startswith(">"):  # a header line
                nonMLST = line[1].isdigit()  # Name of a non-MLST sequence starts with a digit, such as ">86__AmpH_Bla__AmpH__634.consensus".
            if nonMLST:
                o_gene.write(line)
            else:
                o_mlst.write(line)
        fasta.close()
        o_mlst.close()
        o_gene.close()
    
if __name__ == "__main__":
    main()
