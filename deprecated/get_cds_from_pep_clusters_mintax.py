#!/usr/bin/env python

'''
usage: python get_cds_from_pep_clusters.py <directory of peptide cluster fastas> <all cds fasta file> <minimum taxa>
This script needs to be run from the output directory
'''
from Bio import SeqIO
import sys
import re
import os

DIR, reference, mintax = sys.argv[1:]
ref_index = SeqIO.index(reference, "fasta")

for cluster in os.listdir(DIR):
    clust = open(DIR+'/'+cluster, 'rU')
    taxa = []
    for record in SeqIO.parse(clust, "fasta"):
    	rec = record.id
    	species = rec.split('@')[0]
    	if species not in taxa:	
    		taxa.append(species)
    if len(taxa) >= int(mintax):
    	output = open(cluster, "w")
    	for record in SeqIO.parse(clust, "fasta"):
        	output.write(ref_index.get_raw(record.id))