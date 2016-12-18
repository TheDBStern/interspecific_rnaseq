#!/usr/bin/env python

'''
This script creates a fasta file for each species in the analysis. It loops through each orthogroup
fetching all isoforms of Trinity "gene" sequences from the concatenated Trinity.fasta file, renaming to the orthogroup name
'''
from Bio import SeqIO
from Bio import Phylo
import sys
import re
import os

infile, reference = sys.argv[1:]
if DIR[-1] != "/": DIR += "/"
print 'Indexing reference...'
ref_index = SeqIO.index(reference, "fasta")
print 'Extracting sequences'

blastout = open(infile, 'rU')
for line in blastout:
    counter = 0
    name = line.split('\t')[0]
    sequence = line.split('\t')[1]
    species = sequence('@')[0]
    output = open(species+'.blasthits.fasta', 'a')
    gene = '_'.join(sequence.split('_')[:-1])
    trans = []
    for record in ref_index.keys():
        if gene in record and record not in trans: 
            trans.append(record)
    records = (ref_index[tran] for tran in trans)
    for record in records:
        record.id = name+'_'+str(counter)
        record.description = name+'_'+str(counter)
        SeqIO.write(record, output, "fasta")
        counter += 1
            
        
'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''