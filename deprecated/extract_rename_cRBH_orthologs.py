#!/usr/bin/env python

'''
This script takes output from a pairwise (for now) cRBH run, extracts those orthologs from the respective 
Trinity.fasta output file and renames transcripts in the second fasta file to the corresponding ones in the first

usage: extract_rename_cRBH_orthologs.py <cRBH_output.txt> <first_fasta_input> <second_fasta_input> <first_fasta_output> <second_fasta_output>

Note: This assumes cRBH was run on the transdecoder.pep files, filtered to retain only the only ORF per transcript. The cRBH scripts of Salikos and Rochas 2011 
require no special characters in the fasta file names (besides "_"). The original file names should be returned in the cRBH file before running (i.e. TR28398|c0_g1_i1)
'''
from Bio import SeqIO
import sys


outfile_first = open(sys.argv[4], "w")
outfile_second = open(sys.argv[5], "w")
orthos = open(sys.argv[1], "rU")
first_fasta = open(sys.argv[2], "rU")
second_fasta = open(sys.argv[3], "rU")

genes_first = []
genes_second = []
ortho_second = {}

for line in orthos:
    info = line.split(',')
    gene_first = info[0]
    gene_second = info[1].strip('\n')
    genes_first.append(gene_first)
    genes_second.append(gene_second)
    ortho_second[gene_second] = gene_first


for record in SeqIO.parse(first_fasta, "fasta"):
    for x in genes_first:
        if x == record.id:
            SeqIO.write(record, outfile_first, "fasta")
        
for record in SeqIO.parse(second_fasta, "fasta"):
    for x in genes_second:
        if x == record.id:
            record.id = ortho_second[x]
            SeqIO.write(record, outfile_second, "fasta")

outfile_first.close()
outfile_second.close()
orthos.close()
first_fasta.close()
second_fasta.close()


