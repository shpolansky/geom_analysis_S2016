# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 11:24:52 2016

@author: spolansky
"""

coord_dict = {'H0': [1.00, 2.00, 1.00], 'O0': [2.00, 1.00, 1.00], 'H1': [1.00, 2.34, 1.00], 'H2': [3.45, 1.00, 1.00]}

#this is to calculate the angles between the atoms

import math
import itertools
import numpy as np


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
