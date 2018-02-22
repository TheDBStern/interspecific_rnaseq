#!/usr/bin/env python

'''
This script takes a text file with a list of fasta files with overlapping ortholog names 
and outputs the intersection of these for each input file
'''
from Bio import SeqIO
import sys
import re
import os

infile = sys.argv[1]
fasta_list = open(infile, 'rU')
ortho_list = []

for line in fasta_list:
    handle = open(line.strip('\n'), 'rU')
    orthos = []
    fasta = SeqIO.parse(handle, "fasta")    
    for record in fasta:
        if record.id not in ortho_list:
            orthos.append(record.id)
    ortho_list.append(orthos)
    
intersect = set.intersection(*map(set,ortho_list))
print 'Found '+str(len(intersect))+' common genes'
fasta_list.close()
fasta_list = open(infile, 'rU')

for line in fasta_list:
    handle = open(line.strip('\n'), 'rU')
    species = line.split('/')[-1].split('.')[0]
    fasta = SeqIO.parse(handle, "fasta")
    output = open(species+'.intersect.fasta', 'w')
    for record in fasta:
        if record.id in intersect:
            SeqIO.write(record, output, 'fasta')

            