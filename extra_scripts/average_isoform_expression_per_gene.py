#!/usr/bin/env python
''' 
This script takes the isoforms.results output of RSEM and creates a gene-level results file with each 'gene' level abundance
taken as the average across member 'transcripts'
usage: python average_isoform_expression_per_gene.py [input.isoforms.results]
'''
from __future__ import division
import os
import sys
import numpy


def get_gene_list(infile):
    genes = []
    with open(infile, 'r') as iso_input:
        next(iso_input)
        for line in iso_input:
            gene = line.split('\t')[1]
            if gene not in genes:
                genes.append(gene)
    return genes

def average_expression(gene, infile, output):
    lengths = []
    e_lengths = []
    counts = []
    TPMs = []
    FPKMs = []
    isos = []
    with open(infile, 'r') as iso_input:
        for line in iso_input:
            vals = line.split('\t')
            trans_id,gene_id,length,effective_length,expected_count,TPM,FPKM,IsoPct = vals[0:8]
            if gene == gene_id:
                isos.append(trans_id)
                length_w = float(length)*float(IsoPct)/100
                e_length_w = float(effective_length)*float(IsoPct)/100
                lengths.append(length_w)
                e_lengths.append(e_length_w)
                counts.append(float(expected_count))
                TPMs.append(float(TPM))
                FPKMs.append(float(FPKM))
    gene_out = gene +'\t'+','.join(isos)+'\t'+str(round(numpy.sum(lengths),2))+'\t'+str(round(numpy.sum(e_lengths),2))+'\t'+str(round(numpy.mean(counts),2))+'\t'+str(round(numpy.mean(TPMs),2))+'\t'+str(round(numpy.mean(FPKMs),2))+'\n'
    output.write(gene_out) 
        
if __name__ == "__main__":
    infile = sys.argv[1]
    output = open(infile.split('.')[0]+'.averages.results', 'w')
    output.write('gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM\n')
    gene_list = get_gene_list(infile)
    for gene in gene_list:
        average_expression(gene, infile, output)
