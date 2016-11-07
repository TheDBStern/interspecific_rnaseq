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

DIR, reference, ending = sys.argv[1:]
if DIR[-1] != "/": DIR += "/"
print 'Indexing reference...'
ref_index = SeqIO.index(reference, "fasta")
print 'Extracting sequences'

for file in os.listdir(DIR):
    if file.endswith(ending):
        counter = 0
        clust = file.split('.')[0]
        print clust
        tree = Phylo.read(DIR+file, "newick")
        tips = []
        for leaf in tree.get_terminals():
            tips.append(leaf.name)
        for tip in tips:
            species = tip.split('@')[0]
            output = open(species+'.homologs.fasta', 'a')
            gene = '_'.join(tip.split('_')[:-1])
            trans = []
            for record in ref_index.keys():
                if gene in record and record not in trans: 
                    trans.append(record)
            records = (ref_index[tran] for tran in trans)
            for record in records:
                record.id = clust+'_'+str(counter)
                record.description = clust+'_'+str(counter)
                SeqIO.write(record, output, "fasta")
                counter += 1
            
        
'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''