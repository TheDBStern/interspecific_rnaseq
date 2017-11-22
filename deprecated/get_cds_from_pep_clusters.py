#!/usr/bin/env python

'''
usage: python get_cds_from_pep_clusters.py <directory of peptide cluster fastas> <all cds fasta file>
This script needs to be run from the output directory
'''
from Bio import SeqIO
import sys
import re
import os

DIR, reference = sys.argv[1:]
ref_index = SeqIO.index(reference, "fasta")

for cluster in os.listdir(DIR):
    output = open(cluster, "w")
    clust = open(DIR+'/'+cluster, 'rU')
    for record in SeqIO.parse(clust, "fasta"):
        output.write(ref_index.get_raw(record.id))