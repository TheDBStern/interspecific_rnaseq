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
        if re.search(r"TR\d+\|c\d+_g\d+_\i\d+", record.id):
            m = re.search(r"(TR\d+\|c\d+_g\d+)(_\i\d+)", record.id)
            gene = m.group(1)
            trans = m.group(1)+m.group(2)
            outfile.write(gene+'\t'+trans+'\n')
        elif re.search(r"TRINITY_DN\d+\|c\d+_g\d+_i\d+", record.id):
            m = re.search(r"(TRINITY_DN\d+\|c\d+_g\d+)(_i\d+)", record.id)
            gene = m.group(1)
            trans = m.group(1)+m.group(2)
            outfile.write(gene+'\t'+trans+'\n')


if __name__ == "__main__":
    input = sys.argv[1]
    output = sys.argv[2]
    gene_trans(input, output)