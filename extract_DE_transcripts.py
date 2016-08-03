#!/usr/bin/env python

from Bio import SeqIO
import sys

'''
This script takes an output of the analyze_diff_expr.pl script (upregulated transcripts in one sample)
and extracts those from the transdecoder.pep fasta file

Usage: python extract_DE_transcripts.py <.subset file> <.pep fasta> <output>
'''

outfile = open(sys.argv[3], "w")
DE_genes = open(sys.argv[1], "rU")
fasta_in = open(sys.argv[2], "rU")


genes = []

for line in DE_genes:
    info = line.split('\t')
    genes.append(info[0])

for record in SeqIO.parse(fasta_in, "fasta") :
    for x in genes:
        if x in record.id:
            outfile.write(record.format("fasta"))
fasta_in.close()


