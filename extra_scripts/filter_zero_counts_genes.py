#!/usr/bin/env python
''' 
This script takes a matrix of samples x genes counts and outputs a new matrix filtering any gene with a zero count in 
any sample

'''
import os
import sys


with open(sys.argv[1], 'r') as f:
    output = open('.'.join(sys.argv[1].split('.')[0:2])+'.filtered.matrix', 'w')
    incount = 0
    outcount = 0
    for line in f:
        incount += 1
        if '\t0' in line or '\t0\.0' in line or '\t0\.00' in line:
            pass
        else:
            outcount += 1
            output.write(line)

    print "Kept %s of %s genes" % (outcount, incount)

