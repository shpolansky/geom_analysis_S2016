# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 15:45:27 2016

@author: spolansky
"""

import itertools
import math
import numpy as np


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


#part 3 this is to calculate the angles between the atoms


def angle(m, n, s):
    list_m = np.array(coord_dict[m])
    list_n = np.array(coord_dict[n])
    list_s = np.array(coord_dict[s])
    r1 = list_m - list_n
    r2 = list_m - list_s
    dot_product = np.dot(r1, r2)
    magnitude_product = np.linalg.norm(r1) * np.linalg.norm(r2)
    return(math.acos(dot_product / magnitude_product))

for m, n, s in itertools.permutations(coord_dict.keys(), 3):
    print(m, n, s),
    print(angle(m, n, s))

#part A 1 finding the out of plane angles


def outofplaneangle(m, n, s, t):
    m_list = np.array(coord_dict[m])
    n_list = np.array(coord_dict[n])
    s_list = np.array(coord_dict[s])
    t_list = np.array(coord_dict[t])
    r_mn = m_list - n_list
    r_ms = m_list - s_list
    r_mt = m_list - t_list
    e_mn = r_mn / np.linalg.norm(r_mn)
    e_ms = r_ms / np.linalg.norm(r_ms)
    e_mt = r_mt / np.linalg.norm(r_mt)
    angle_v = angle(m, n, s)
    return(math.asin(np.dot(((np.cross(e_mn, e_ms)) / math.sin(angle_v)), e_mt)))

for m, n, s, t in itertools.permutations(coord_dict.keys(), 4):
    print(m, n, s, t),
    print(outofplaneangle(m, n, s, t))
