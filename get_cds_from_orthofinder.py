#!/usr/bin/env python

'''
usage: python get_cds_from_orthofinder.py 
This script needs to be run from the output directory
'''
from Bio import SeqIO
import sys
import re
import os

reference, input, MinTax = sys.argv[1:]
ref_index = SeqIO.index(reference, "fasta")
print 'Extracting sequences'

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
		print species
		if len(species) >= int(MinTax):
			output = open(clust+'.fasta', 'w')
			print clust
			tips = []
			for leaf in line.split(' ')[1:]:
				tips.append(leaf.strip('\n'))
			for tip in tips:
				output.write(ref_index.get_raw(tip))
		else:
			pass
