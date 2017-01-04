#!/usr/bin/env python

'''
This script creates a fasta file for each species in the analysis. It loops through each orthogroup
fetching all isoforms of Trinity "gene" sequences from the concatenated Trinity.fasta file, renaming to the orthogroup name
'''
from Bio import SeqIO
import sys
import re
import os

infile, reference = sys.argv[1:]
print 'Indexing reference...'
ref_index = SeqIO.index(reference, "fasta")
print 'Extracting sequences'

count_rbbh = 0
blastout = open(infile, 'rU')
for line in blastout:
    seq1 = line.split('\t')[0]
    species1 = seq1.split('@')[0]
    seq2 = line.split('\t')[1]
    species2 = seq2.split('@')[0]
    output1 = open(species1+'.rbbh.fasta', 'a')
    output2 = open(species2+'.rbbh.fasta', 'a')
    rec1 = ref_index[seq1]
    rec2 = ref_index[seq2]
    rec1.id = 'RBBH'+str(count_rbbh)
    rec1.description = 'RBBH'+str(count_rbbh)
    SeqIO.write(rec1, output1, "fasta")
    rec2.id = 'RBBH'+str(count_rbbh)
    rec2.description = 'RBBH'+str(count_rbbh)
    SeqIO.write(rec2, output2, "fasta")
    count_rbbh += 1
            
        
'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''