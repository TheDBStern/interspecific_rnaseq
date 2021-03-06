#!/usr/bin/env python

'''
This script takes output from a pairwise (for now) RSD run, extracts those orthologs from the respective 
Trinity.fasta output file and renames transcripts in the second fasta file to the corresponding ones in the first
This will extract all isoforms of each gene.
The first line of the RSD output file (or any line that doesn't contain the two orthologs) should be removed
Each line should be in this format 'OR\tTR99751|c0_g1_i1\tTR157456|c0_g1_i2\t5.029'

usage: extract_rename_RSD_orthologs.py <RSD_output.txt> <first_fasta_input> <second_fasta_input> <first_fasta_output> <second_fasta_output>

'''
from Bio import SeqIO
import sys
import re



outfile_first = open(sys.argv[4], "w")
outfile_second = open(sys.argv[5], "w")
orthos = open(sys.argv[1], "rU")
first_fasta = open(sys.argv[2], "rU")
second_fasta = open(sys.argv[3], "rU")

genes_first = []
genes_second = []
ortho_second = {}

for line in orthos:
    info = line.split('\t')
    gene_first = info[1]
    gene_first = re.sub(r'_i\d+', '', gene_first)
    gene_second = info[2]
    gene_second = re.sub(r'_i\d+', '', gene_second)
    genes_first.append(gene_first)
    genes_second.append(gene_second)
    ortho_second[gene_second] = gene_first


for record in SeqIO.parse(first_fasta, "fasta"):
    gene = re.sub(r'_i\d+', '', record.id)
    if gene in genes_first: 
        print ('Writing: ' + record.id)
        SeqIO.write(record, outfile_first, "fasta")

gene_list = []
for record in SeqIO.parse(second_fasta, "fasta"):
    gene = re.sub(r'_i\d+', '', record.id)
    gene_list.append(gene)
    if gene in genes_second:
        record.id = ortho_second[gene]  + '_i' + str(gene_list.count(gene))
        print ('Writing: ' + record.id)
        SeqIO.write(record, outfile_second, "fasta")

outfile_first.close()
outfile_second.close()
orthos.close()
first_fasta.close()
second_fasta.close()


