# interspecific_rnaseq

These are scripts used to assemble data and analyze gene expression evolution in 'The evolution of gene expression underlying vision loss in cave animals'

* fasta2tree_w.py Takes a directory of fasta files, aligns them using translatorX, and estimates a tree using raxml, partitioned by codon. Alternatively, can specify to use iqtree with modeltesting. Need to adjust translatorx script path and raxml command. Needs mafft in the PATH
* get_cds_from_orthofinder.py Creates a fasta file for each orthogroup. Requires a concatenated cds file of all species and the Orthogroups.txt output from Orthofinder
* get_expression_data_from_orthofiner.py pulls and renames isoform expression levels from RSEM output for orthogroups in the OrthoFinder Orthogroups.txt file
* get_gene_trans_map.py Creates a tab delimited gene-transcipt map file from a Trinity assembly
* get_longest_protein_seq_per_orthogroup.py Loops through orthogroup fasta files and pulls the longest protein sequence in each from a concatenated sequence file
* get_reference_transcripts_from_orthofinder_orthogroups.py Creates a 'transcriptome' assembly file for each species with sequences renamed to orthogroups
* make_data_matrices_ouwie.py Takes an expression matrix and a .csv file with species states and creates an OUwie input file for each orthogroup
* run_ouwie.R Scripts to run OUwie on each orthogroup

* The data directory contains a species phylogeny, normalized expression matrix and 'states.csv' file with 'blind/sighted' designation for the 14 species (34 individuals)

Species Key  
CCRYP - Cambarus cryptodytes  
CDUBI - Cambarus dubius  
CGRAY - Cambarus graysoni  
CHAMU - Cambarus hamulatus  
CNERT - Cambarus nerterius  
CRUST - Cambarus rusticiformis  
CSETO - Cambarus setosus  
CTENE - Cambarus tenebrosus  
OAUST - Orconectes australis  
OINCO - Orconectes incomptus  
PFALL - Procambarus fallax  
PHORS - Procambarus horsti  
PLUCI - Procambarus lucifugus  
PPALL - Procambarus pallidus  
