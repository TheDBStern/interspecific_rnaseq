#!/usr/bin/env python
''' 
This script takes a directory of fasta files, aligns them using translatorX, and estimates a tree using raxml, partitioned by codon
Need to adjust translatorx script path and raxml command
Needs mafft in the PATH
'''
import glob, os
import sys
from Bio import AlignIO

translatorx_path = '/home/dbstern/Programs/translatorx_vLocal.pl'
raxml_cmd = 'raxmlHPC-PTHREADS-AVX -T 16'

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
    
def raxml(alignment,DIR):
    cluster = alignment.split('.')[0]
    command = raxml_cmd+" -f d -m GTRGAMMA -p 1293049 -# 1 -q "+DIR+"partition.PART -s "+DIR+alignment+" -w "+os.path.abspath(DIR)+" -n "+cluster
    print "executing: " + command
    os.system(command)
    tree = DIR+cluster+".raxml.tre"
    raw_tree = DIR+"RAxML_bestTree."+cluster
    try:
        os.rename(raw_tree,tree)
        os.remove(DIR+"RAxML_info."+cluster)
        os.remove(DIR+"RAxML_log."+cluster)
        os.remove(DIR+"RAxML_parsimonyTree."+cluster)
        os.remove(DIR+"RAxML_result."+cluster)
        os.remove(DIR+cluster+".reduced")
    except:pass 


if __name__ == "__main__":
    DIR, ending = sys.argv[1:]
    if DIR[-1] != "/": DIR += "/"
    for fasta_file in os.listdir(DIR):
        if fasta_file.endswith(ending):
            cluster = fasta_file.split('.')[0]
            transx(fasta_file,DIR)
            os.system('rm '+DIR+'*.html '+DIR+'*nt1_ali.fasta '+DIR+'*nt2_ali.fasta '+DIR+'*nt3_ali.fasta '+DIR+'*nt12_ali.fasta '+DIR+'*.log '+DIR+'*.aaseqs '+DIR+'*.aaseqs.fasta '+DIR+'*aa_based_codon* '+DIR+'*.aa_ali.fasta')
    for fasta_file in os.listdir(DIR):
        if fasta_file.endswith(ending):
            cluster = fasta_file.split('.')[0]
            partition(cluster+".transx.nt_ali.fasta",DIR)
            raxml(cluster+".transx.nt_ali.fasta",DIR)
            