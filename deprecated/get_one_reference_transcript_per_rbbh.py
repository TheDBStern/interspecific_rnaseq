#!/usr/bin/env python

'''
'''
from Bio import SeqIO
import sys
import re
import os

infile, reference = sys.argv[1:]
print 'Indexing reference...'
ref_index = SeqIO.index(reference, "fasta")
print 'Extracting sequences'

def pairwise(blastout):
    count_rbbh = 0
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

def iterative(blastout):
    for line in blastout:
        seq1 = line.split('\t')[0]
        seq2 = line.split('\t')[1]
        if seq1.startswith('RBBH'):
            rbbh_num = seq1
            seq = seq2
            species = seq2.split('@')[0]
        else:
            rbbh_num = seq2
            seq = seq1
            species = seq1.split('@')[0]
        output = open(species+'.rbbh.fasta', 'a')
        rec = ref_index[seq]
        rec.id = rbbh_num
        rec.description = rbbh_num
        SeqIO.write(rec, output, "fasta")    

blastout = open(infile, 'rU')        
first_line = blastout.readline() 
sequence1 = first_line.split('\t')[0]     
sequence2 = first_line.split('\t')[1]  
if sequence1.startswith('RBBH') or sequence2.startswith('RBBH'):
    blastout.close()
    blastout = open(infile, 'rU')
    iterative(blastout)
else:
    blastout.close()
    blastout = open(infile, 'rU')
    pairwise(blastout)

'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''