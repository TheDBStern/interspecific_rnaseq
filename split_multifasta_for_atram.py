#!/usr/bin/env python

'''
Takes a fasta file of target sequences for atram and creates a directory of
 separate fasta for each target as well as the tab-delimited TargetFile.txt
requires biopython
'''
import sys, os
from Bio import SeqIO

input, binsize = sys.argv[1:]
if not os.path.isdir('./Targets'):
    os.mkdir('./Targets')
cwd = os.getcwd()
outtext = open('TargetFile.txt', 'w')
handle = open(input, "rU")
for record in SeqIO.parse(handle, "fasta"):
    gene = record.id
    output = 'Targets/'+gene+'.fasta'
    SeqIO.write(record, output, "fasta")
    outtext.write(gene+'\t'+cwd+'/'+output+'\n')