#!/usr/bin/env python

import sys, os
import argparse
import pandas as pd
import numpy

parser = argparse.ArgumentParser(description='Script to create data matrices for input to OUwie')
parser.add_argument('-s', dest = 'SpFile', type = str, required=True,  help = 'Path to comma-delimited file with the names of all the species in the analysis separated by the state of interest (0 or 1)')
parser.add_argument('-i', dest= 'In', type = str, required=True, help ='Expression matrix')
parser.add_argument('--log2', dest= 'lg2', action ='store_true', default= False, help ='log2 transform expression values, default = False.')
parser.add_argument('--sqrt', dest= 'sqr', action ='store_true', default= False, help ='sqrt transform expression values, default = False.')
args = parser.parse_args()

sp_list = []
lib_dict = {}
state_dict = {}

exp_mat = pd.read_csv(args.In, sep='\t', index_col=0)
for lib in list(exp_mat):
	if lib.split('_')[0] not in sp_list:
		sp_list.append(lib.split('_')[0])
for sp in sp_list:
	for lib in list(exp_mat):
		if sp in lib:
			lib_dict.setdefault(sp, []).append(lib)
with open(args.SpFile, 'rU') as f:
	for line in f:
		state_dict[line.split(',')[0]] = line.split(',')[1].strip('\n')

for index in exp_mat.iterrows():
	print index[0]
	if args.lg2:
		outfile = open('ouwie_dat_sims/'+str(index[0])+'.dat.log2.csv','w')
	if args.sqr:
		outfile = open('ouwie_dat_sims/'+str(index[0])+'.dat.sqrt.csv','w')
	else:
		outfile = open('ouwie_dat_sims/'+str(index[0])+'.dat.csv','w')
	outfile.write('Species,State,ExpSE\n')
	for sp in sp_list:
		if args.lg2:
			mymean = numpy.mean(numpy.log2(exp_mat.loc[index[0],lib_dict[sp]]+1))
			myse = numpy.std(numpy.log2(exp_mat.loc[index[0],lib_dict[sp]]+1)) / numpy.sqrt(len(exp_mat.loc[index[0],lib_dict[sp]]))
			outfile.write(','.join([sp,state_dict[sp],str(myse+0.000000001)+'\n']))
		if args.sqr:
			mymean = numpy.mean(numpy.sqrt(exp_mat.loc[index[0],lib_dict[sp]]))
			myse = numpy.std(numpy.sqrt(exp_mat.loc[index[0],lib_dict[sp]])) / numpy.sqrt(len(exp_mat.loc[index[0],lib_dict[sp]]))
			outfile.write(','.join([sp,state_dict[sp],str(myse+0.000000001)+'\n']))		
		else:
			mymean = numpy.mean(exp_mat.loc[index[0],lib_dict[sp]])
			myse = numpy.std(exp_mat.loc[index[0],lib_dict[sp]])/numpy.sqrt(len(exp_mat.loc[index[0],lib_dict[sp]]))
			outfile.write(','.join([sp,state_dict[sp],str(myse)+'\n']))