#!/usr/bin/env python

from Bio import SeqIO
from Bio import Phylo
import sys
import re
import os
import argparse

parser = argparse.ArgumentParser(description='This script pulls and renames isoform expression levels from RSEM output for orthogroups in the OrthoFinder Orthogroups.txt file')
parser.add_argument('-d', dest = 'input', type = str, required=True,  help = 'Orthogroups.txt file')
parser.add_argument('-m', dest = 'MinTax', type = int, required=True,  help = 'Minimum number of taxa')
parser.add_argument('-i', dest= 'datDIR', required=True, help ='Directory with RSEM output files (indiv.isoforms.results)')
parser.add_argument('-o', dest= 'outDIR', required=True, help ='Output directory')

args = parser.parse_args()

if args.datDIR[-1] != "/": args.datDIR += "/"
if args.outDIR[-1] != "/": args.outDIR += "/"

with open(args.input,'r') as f:
	for line in f:
		counter = 0
		clust = line.split(' ')[0].strip(':')
		print clust
		species = []
		for leaf in line.split(' ')[1:]:
			sp = leaf.split('@')[0].strip('\n')
			if sp not in species:
				species.append(sp)
		if len(species) >= int(args.MinTax):
			print clust
			print species
			tips = []
			for leaf in line.split(' ')[1:]:
				tips.append(leaf.strip('\n'))
			for tip in tips:
				species = tip.split('@')[0]
				iso = tip.split('@')[1]
				for datfile in os.listdir(args.datDIR):
					if datfile.split('_')[0] == species:
						output = open(args.outDIR+datfile, 'a')
						with open(args.outDIR+datfile,'r') as f:
							head = f.readline()
							if head.split('\t')[0] != 'transcript_id':
								output.write('transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct\n')	
						for line in open(args.datDIR+datfile,'r'):
							if line.split('\t')[0] == iso:
								output.write(line.replace(iso,clust+'_'+str(counter)))	
				counter += 1
			
'''
transcript = 'CACUM@TRINITY_DN_c0_g1_i1'
gene = '_'.join(transcript.split('_')[:-1])
'''