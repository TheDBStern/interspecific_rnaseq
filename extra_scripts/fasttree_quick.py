import glob, os
import sys
from Bio import AlignIO

def seqcount(alignment):
    count = 0
    for fasta in alignment:
        count += 1
    return count

def tree(alignment,DIR):
    ali = AlignIO.read(DIR+alignment, "fasta")
    if seqcount(ali) > 500:
        cluster = alignment.split('.')[0]
        command = 'FastTree -gtr -nt '+DIR+alignment+' > '+DIR+cluster+'.fasttree.tre'
        os.system(command)


if __name__ == "__main__":
    DIR, ending = sys.argv[1:]
    if DIR[-1] != "/": DIR += "/"
    for fasta_file in os.listdir(DIR):
        if fasta_file.endswith(ending):
            cluster = fasta_file.split('.')[0]
            tree(cluster+".transx.nt_ali.fasta",DIR)
