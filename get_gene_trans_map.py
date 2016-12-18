#!/usr/bin/env python

'''
Python version of Trinity's get gene-trans map
usage: python get_gene_trans_map.py <Trinity.fasta> <Trinity.fasta.gene_trans_map>
Requires BioPython
'''
import sys, os
import re
from Bio import SeqIO
 
def gene_trans(fasta, output):
    handle = open(fasta, "rU")
    outfile = open(output, "w")
    for record in SeqIO.parse(handle, "fasta"):
        transcript = record.id
        gene = '_'.join(transcript.split('_')[:-1])
        outfile.write(gene+'\t'+transcript+'\n')

if __name__ == "__main__":
    input = sys.argv[1]
    output = sys.argv[2]
    gene_trans(input, output)

'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''