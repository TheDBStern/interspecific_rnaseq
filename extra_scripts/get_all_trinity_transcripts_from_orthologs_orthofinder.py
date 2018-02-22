#!/usr/bin/env python

'''
This script creates a fasta file for each species in the analysis. It loops through each orthogroup
fetching all isoforms of Trinity "gene" sequences from the concatenated Trinity.fasta file, renaming to the orthogroup name
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
		if len(species) == int(MinTax):
			print clust
			tips = []
			for leaf in line.split(' ')[1:]:
				tips.append(leaf.strip('\n'))
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
		else:
			pass
		
'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''