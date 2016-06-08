# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 15:45:27 2016

@author: spolansky
"""

coord_dict = {}

with open('test.xyz', 'r') as f:
    n_atoms = int(f.readline())
    symbol_list = []
    coord_list = []
    for line in f.readlines()[1: n_atoms + 3]:
        atom_list = line.split()
        atom_float_list = [float(x) for x in atom_list[1: 4]]
        coord_list.append(atom_float_list)
        Occ = 0
        for i in range(len(symbol_list)):
            if symbol_list[i][0] == atom_list[0]:
                Occ += 1
        atom_list[0] = atom_list[0] + str(Occ)
        symbol_list.append(atom_list[0])
        coord_dict[symbol_list[-1]] = coord_list[-1]

print(coord_dict)



