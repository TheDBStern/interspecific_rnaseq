#!/usr/bin/env python

'''
'''
from Bio import SeqIO
import sys
import re
import os

DIR, reference, ending = sys.argv[1:]
ref_index = SeqIO.index(reference, "fasta")

output = open('orthogroups.fasta','w')
for cluster in os.listdir(DIR):
	if cluster.endswith(ending):
		clust = open(DIR+'/'+cluster, 'rU')
		og = cluster.split('/')[-1].strip('.fa')
		longest = ''
		length = 0
		for record in SeqIO.parse(clust, "fasta"):
			if len(record) > length:
				longest = record.id
				length = len(record)
		record = ref_index[longest]
		record.id = og
		record.description = og
		SeqIO.write(record, output, "fasta")