# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 15:45:27 2016

@author: spolansky
"""

import itertools


#part 1 reading the file

coord_dict = {}
#dictionary with element symbols as keys and coordinates as values

with open('test.xyz', 'r') as f:
    n_atoms = int(f.readline())
    symbol_list = []    #this is the list of the elemental symbols
    coord_list = []     #this is the list of the coordinates
    for line in f.readlines()[1: n_atoms + 3]:
        atom_list = line.split()
        atom_float_list = [float(x) for x in atom_list[1: 4]]
        coord_list.append(atom_float_list)    #at this point the coordinates are ready
        Occ = 0
        for i in range(len(symbol_list)):  #this is where the the symbols are checked for repetitions
            if symbol_list[i][0] == atom_list[0]:
                Occ += 1
        atom_list[0] = atom_list[0] + str(Occ)
        symbol_list.append(atom_list[0])
        coord_dict[symbol_list[-1]] = coord_list[-1]  #this is where the dictionary is appended


#part 2 finding distances

def distance(m, n):
    zdiff = coord_dict[m][2] - coord_dict[n][2]
    ydiff = coord_dict[m][1] - coord_dict[n][1]
    xdiff = coord_dict[m][0] - coord_dict[n][0]
    return(((xdiff ** 2) + (ydiff ** 2) + (zdiff ** 2)) ** 0.5)

#itertools allows every unique combination of two keys from the dictionary to be used to find the distances

for m, n in itertools.combinations(coord_dict.keys(), 2):   
    print(m, n),
    print(distance(m, n))        


