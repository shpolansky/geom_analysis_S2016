# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 11:24:52 2016

@author: spolansky
"""

coord_dict = {'H0': [1.00, 2.00, 1.00], 'O0': [2.00, 1.00, 1.00], 'H1': [1.00, 2.34, 1.00], 'H2': [3.45, 1.00, 1.00]}

#this is to calculate the angles between the atoms

import math
import itertools

def angle(m, n, s):
    r1 = []
    r2 = []
    r1_xdiff = coord_dict[m][0] - coord_dict[n][0]
    r1_ydiff = coord_dict[m][1] - coord_dict[n][1]
    r1_zdiff = coord_dict[m][2] - coord_dict[n][2]
    r2_xdiff = coord_dict[m][0] - coord_dict[s][0]
    r2_ydiff = coord_dict[m][1] - coord_dict[s][1]
    r2_zdiff = coord_dict[m][2] - coord_dict[s][2]
    r1.append(r1_xdiff)
    r1.append(r1_ydiff)
    r1.append(r1_zdiff)
    r2.append(r2_xdiff)
    r2.append(r2_ydiff)
    r2.append(r2_zdiff)
    dot_product = (r1[0] * r2[0]) + (r1[1] * r2[1]) + (r1[2] * r2[2])
    magnitude_product = (((r1[0] ** 2) + (r1[1] ** 2) + (r1[2] ** 2)) ** 0.5) * (((r2[0] ** 2) + (r2[1] ** 2) + (r2[2] ** 2)) ** 0.5)
    angle = math.acos(dot_product / magnitude_product)
    return(angle)
    
for m, n, s in itertools.combinations(coord_dict.keys(), 3):
    print(m, n, s),
    print(angle(m, n, s))