# /usr/bin/python
#

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

import argparse

# Read-in command-line optins
parser = argparse.ArgumentParser(description='Process some command-line options')

parser.add_argument("-f", "--file", help="File stem to use when opening input file.")
parser.add_argument("-o", "--output", help="File stem to use when opening output file.")
parser.add_argument("-i", "--index", help="Sweep index.")

args = parser.parse_args()

# Read in data
data = np.genfromtxt("data/data_4/{file}.csv".format(file=args.file), dtype=float, delimiter='\t')

# Define useful arrays
data_offset = 505
i = 0

with open("{output}.csv".format(output=args.output),"wb") as outfile:

    while ((data_offset*i + 1) < np.size(data)):

        config = data[:5 + data_offset*i]
        lower = 224.37 - 0.4
        upper = 224.37 + 0.4        

        delta = ((upper - lower)/501)
        
        x = [lower+j*delta for j in range(501)]
        y = data[5 + data_offset*i:505 + data_offset*i]
        
        y = np.insert(y, 0, args.index)
        
        np.savetxt(outfile, y.reshape(1, 501), delimiter="\t", fmt="%.5f")
        print("%f of %f" % (data_offset*i, np.size(data)) )
        i = i + 1

