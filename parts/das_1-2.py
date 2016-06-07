import numpy as np
import itertools as it
import pandas as pd
from scipy.spatial.distance import cdist
import argparse

# Setup input
parser = argparse.ArgumentParser(description='Performs geometry analysis \
for a given molecule in .xyz format')
parser.add_argument('-x', '--xyz_input', help='File containing geometry \
of molecule in .xyz format')

args = parser.parse_args()
f = args.xyz_input

# Build pandas dataframe from given .xyz file
xyz_df = pd.read_csv('methane.xyz', delim_whitespace=True, header=None, 
names=['x','y','z'], index_col=0, dtype={'x':np.float64, 
'y':np.float64, 'z':np.float64}, skiprows=2) 

# Update labels list
labels = xyz_df.index.values
uni_lab = labels
for i in xrange(len(uni_lab)):
    occs = 0
    for j in xrange(len(labels)):
        if uni_lab[i] == labels[j] and i != j:
            occs +=1
    uni_lab[i] += str(occs)
xyz_df.index = uni_lab

# Build interatomic distances dataframe
R_df = pd.DataFrame(cdist(xyz_df.values, xyz_df.values, 'euclidean'), 
index=xyz_df.index.values, columns=xyz_df.index.values)
print R_df
