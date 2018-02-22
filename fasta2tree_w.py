#!/usr/bin/env python

import glob, os
import sys
from Bio import AlignIO
import argparse

parser = argparse.ArgumentParser(description='This script takes a directory of fasta files, aligns them using translatorX, and estimates a tree using raxml, partitioned by codon. Alternatively, can specify to use iqtree with modeltesting. Need to adjust translatorx script path and raxml command. Needs mafft in the PATH')
parser.add_argument('-D', dest = 'DIR', type = str, required=True,  help = 'Path to input directory of fasta files')
parser.add_argument('-E', dest = 'ending', type = str, required=True,  help = 'Fasta files ending')
parser.add_argument('--iqtree', dest= 'iqtree', action ='store_true', default= False, help ='use iqtree to estimate phylogeny with modeltesting (TESTNEW), default = False.')
args = parser.parse_args()

<<<<<<< HEAD
translatorx_path = '/Users/dbstern/Desktop/Phylogenetics_Programs/translatorx_vLocal.pl'
raxml_cmd = 'raxmlHPC-PTHREADS-AVX -T 4'
iqtree_cmd = 'iqtree'
=======
translatorx_path = '/home/dbstern/Programs/translatorx_vLocal.pl'
raxml_cmd = 'raxmlHPC-PTHREADS-AVX -T 16'
>>>>>>> 6b3501c7cdc319ab6692ee9fd42dd95c8cebd5a6

def transx(fasta_file,DIR):
	transout = fasta_file.split('.')[0]+".transx"
	command = "perl "+translatorx_path+" -i "+DIR+fasta_file+" -o "+DIR+transout+" -p F -c 1 -t T"
	print "executing: " + command
	os.system(command)
	
def partition(alignment,DIR):
	ali = AlignIO.read(DIR+alignment, "fasta")
	length = str(ali.get_alignment_length())
	output = open(DIR+"partition.PART", "w")
	output.write("DNA, p1 = 1-"+length+"\\3\nDNA, p2 = 2-"+length+"\\3\nDNA, p3 = 3-"+length+"\\3\n")

def seqcount(alignment):
	count = 0
	for fasta in alignment:
		count += 1
	return count
	
def tree(alignment,DIR):
	ali = AlignIO.read(DIR+alignment, "fasta")
	if seqcount(ali) < 500:
		cluster = alignment.split('.')[0]
		if args.iqtree:
			command = iqtree_cmd+' -s '+DIR+alignment+' -spp '+DIR+'partition.PART -m TESTNEWMERGE -nt AUTO -pre '+DIR+cluster
			os.system(command)
			tree = DIR+cluster+".iqtree.tre"
			raw_tree = DIR+cluster+'.treefile'
			try:
				os.rename(raw_tree,tree)
				os.system('rm '+DIR+cluster+'.best_scheme '+DIR+cluster+".best_scheme.nex "+DIR+cluster+".bionj "+DIR+cluster+".ckp.gz "+DIR+cluster+".iqtree "+DIR+cluster+".log "+DIR+cluster+".mldist "+DIR+cluster+".model.gz")
			except:pass 
		else:
			command = raxml_cmd+" -f d -m GTRCAT -p 1293049 -# 3 -q "+DIR+"partition.PART -s "+DIR+alignment+" -w "+os.path.abspath(DIR)+" -n "+cluster
			print "executing: " + command
			os.system(command)
			tree = DIR+cluster+".raxml.tre"
			raw_tree = DIR+"RAxML_bestTree."+cluster
			try:
				os.rename(raw_tree,tree)
				os.system('rm '+DIR+"RAxML_*")
				os.system('rm '+DIR+cluster+".reduced*")
			except:pass 
	else:
		cluster = alignment.split('.')[0]
		command = 'FastTree -gtr -nt '+DIR+alignment+' > '+DIR+cluster+'.fasttree.tre'
		os.system(command)

if __name__ == "__main__":
	DIR = args.DIR
	ending = args.ending
	if DIR[-1] != "/": DIR += "/"
	for fasta_file in os.listdir(DIR):
		if fasta_file.endswith(ending):
			cluster = fasta_file.split('.')[0]
			if os.path.isfile(DIR+cluster+'.transx.nt_ali.fasta'):
				print('Detected alignment for '+cluster)
			else:
				transx(fasta_file,DIR)
				os.system('rm '+DIR+'*.html '+DIR+'*nt1_ali.fasta '+DIR+'*nt2_ali.fasta '+DIR+'*nt3_ali.fasta '+DIR+'*nt12_ali.fasta '+DIR+'*.log '+DIR+'*.aaseqs '+DIR+'*.aaseqs.fasta '+DIR+'*aa_based_codon* '+DIR+'*.aa_ali.fasta')
	for fasta_file in os.listdir(DIR):
		if fasta_file.endswith(ending):
			cluster = fasta_file.split('.')[0]
<<<<<<< HEAD
			if os.path.isfile(DIR+cluster+'.raxml.tre') or os.path.isfile(DIR+cluster+'.fasttree.tre') or os.path.isfile(DIR+cluster+'.iqtree.tre'): 
				pass
=======
			if os.path.isfile(DIR+cluster+'.raxml.tre') or os.path.isfile(DIR+cluster+'.fasttree.tre'): 
				print('Detected tree file for '+cluster)
>>>>>>> 6b3501c7cdc319ab6692ee9fd42dd95c8cebd5a6
			else:
				partition(cluster+".transx.nt_ali.fasta",DIR)
				tree(cluster+".transx.nt_ali.fasta",DIR)