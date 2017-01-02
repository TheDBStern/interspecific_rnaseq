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

blastout = open(infile, 'rU')
for line in blastout:
    count_rbbh == 0
    seq1 = line.split('\t')[0]
    species1 = sequence.split('@')[0]
    seq2 = line.split('\t')[1]
    species2 = sequence.split('@')[1]
    output1 = open(species1+'.rbbh.fasta', 'a')
    output2 = open(species2+'.rbbh.fasta', 'a')
    gene1 = '_'.join(sequence1.split('_')[:-1])
    gene2 = '_'.join(sequence2.split('_')[:-1])
    trans1 = []
    for record in ref_index.keys():
        if gene1 in record and record not in trans1: 
            trans1.append(record)
    for record in ref_index.keys():
        if gene2 in record and record not in trans2: 
            trans1.append(record)
    records1 = (ref_index[tran] for tran in trans1)
    records2 = (ref_index[tran] for tran in trans2)
    for record in records1:
        record.id = 'RBBH'+str(count_rbbh)+'_'+str(records1.index(record))
        record.description = 'RBBH'+str(count_rbbh)+'_'+str(records1.index(record))
        SeqIO.write(record, output1, "fasta")
    for record in records2:
        record.id = 'RBBH'+str(count_rbbh)+'_'+str(records2.index(record))
        record.description = 'RBBH'+str(count_rbbh)+'_'+str(records2.index(record))
        SeqIO.write(record, output1, "fasta")
    count_rbbh += 1
            
        
'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''