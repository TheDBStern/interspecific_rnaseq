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
print 'Indexing reference...'
ref_index = SeqIO.index(reference, "fasta")
print 'Extracting sequences'

blastout = open(infile, 'rU')
for line in blastout:
    name = line.split('\t')[0]
    sequence = line.split('\t')[1]
    species = sequence.split('@')[0]
    output = open(species+'.blasthits.fasta', 'a')
    rec = ref_index[sequence]
    rec.id = name
    rec.description = name
    SeqIO.write(rec, output, "fasta")
            
        
'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''