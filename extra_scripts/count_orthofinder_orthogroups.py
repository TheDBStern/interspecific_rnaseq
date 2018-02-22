#!/usr/bin/env python

'''
This script creates a fasta file for each species in the analysis. It loops through each orthogroup
fetching the Trinity sequences from the concatenated Trinity.fasta file, renaming to the orthogroup name
'''
from Bio import SeqIO
import sys
import re
import os
import numpy

input, MinTax = sys.argv[1:]

lengths = []

with open(input) as f:
	for line in f:
		clust = line.split(' ')[0].strip(':')
		counter = 0
		#print clust
		species = []
		for leaf in line.split(' ')[1:]:
			sp = leaf.split('@')[0].strip('\n')
			if sp not in species:
				species.append(sp)
		#print species
		if len(species) >= int(MinTax):
			#print clust
			tips = []
			for leaf in line.split(' ')[1:]:
				tips.append(leaf.strip('\n'))
			lengths.append(len(tips))
			print len(tips)
		else:
			pass

print numpy.mean(lengths)
count =0
for i in lengths:
	if i == 14:
		count +=1
print count
'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''