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
    rec_count = 0
    seq1 = line.split('\t')[0]
    species1 = seq1.split('@')[0]
    seq2 = line.split('\t')[1]
    species2 = seq2.split('@')[0]
    output1 = open(species1+'.rbbh.fasta', 'a')
    output2 = open(species2+'.rbbh.fasta', 'a')
    gene1 = '_'.join(seq1.split('_')[:-1])
    gene2 = '_'.join(seq2.split('_')[:-1])
    trans1 = []
    trans2 = []
    for record in ref_index.keys():
        if gene1 in record and record not in trans1: 
            trans1.append(record)
    for record in ref_index.keys():
        if gene2 in record and record not in trans2: 
            trans2.append(record)
    records1 = (ref_index[tran] for tran in trans1)
    records2 = (ref_index[tran] for tran in trans2)
    for record in records1:
        record.id = 'RBBH'+str(count_rbbh)+'_'+str(rec_count)
        record.description = 'RBBH'+str(count_rbbh)+'_'+str(rec_count)
        SeqIO.write(record, output1, "fasta")
        rec_count += 1
    for record in records2:
        record.id = 'RBBH'+str(count_rbbh)+'_'+str(rec_count)
        record.description = 'RBBH'+str(count_rbbh)+'_'+str(rec_count)
        SeqIO.write(record, output2, "fasta")
        rec_count += 1
    count_rbbh += 1
            
        
'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''