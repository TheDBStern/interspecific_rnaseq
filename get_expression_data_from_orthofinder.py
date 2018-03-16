#!/usr/bin/env python

'''
'''
from Bio import SeqIO
from Bio import Phylo
import sys
import re
import os

input, MinTax,datDIR, outDIR = sys.argv[1:]

if datDIR[-1] != "/": datDIR += "/"
if outDIR[-1] != "/": outDIR += "/"

with open(input,'r') as f:
	for line in f:
		counter = 0
		clust = line.split(' ')[0].strip(':')
		print clust
		species = []
		for leaf in line.split(' ')[1:]:
			sp = leaf.split('@')[0].strip('\n')
			if sp not in species:
				species.append(sp)
		if len(species) >= int(MinTax):
			print clust
			print species
			tips = []
			for leaf in line.split(' ')[1:]:
				tips.append(leaf.strip('\n'))
			for tip in tips:
				species = tip.split('@')[0]
				iso = tip.split('@')[1]
				for datfile in os.listdir(datDIR):
					if datfile.split('_')[0] == species:
						output = open(outDIR+datfile, 'a')
						with open(outDIR+datfile,'r') as f:
							head = f.readline()
							if head.split('\t')[0] != 'transcript_id':
								output.write('transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct\n')	
						for line in open(datDIR+datfile,'r'):
							if line.split('\t')[0] == iso:
								output.write(line.replace(iso,clust+'_'+str(counter)))	
				counter += 1
			
'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''