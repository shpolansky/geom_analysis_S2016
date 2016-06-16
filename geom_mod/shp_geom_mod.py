"""Module containing functions for 
geometry analysis of a set of coordinates
for any molecule."""

import itertools
import math
import numpy as np


#part 1 reading the file

def readgeom(filename):
    """Returns a dictionary of indexed elements with xyz coordinates

    Parameters:
            an xyz file entered as a string

    Return:
            a dictionary of indexed elements with xyz coordinates"""
    coord_dict = {}
    with open(filename, 'r') as f:
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
    return(coord_dict{})


#part 2 finding distances

def distance(m, n):
    """Returns the distance between two atomic coordinates

    Parameters:
            two lists of xyz coordinates
   
    Return:
            the distance value between the coordinates"""
    zdiff = m[2] - n[2]
    ydiff = m[1] - n[1]
    xdiff = m[0] - n[0]
    return(((xdiff ** 2) + (ydiff ** 2) + (zdiff ** 2)) ** 0.5)

def alldistances(dictionary):
    """Returns all of the distances from a dictionary of coordinates

    Parameters:
            takes a dictionary of indexed atoms and their coordinates

    Return:
            a dictionary of the atoms and their distances"""
    distances = {}
    for m, n in itertools.combinations(dictionary.keys(), 2):   
        distances[m, n] = distance(dictionary[m], dictionary[n])
    return(distances{})

#part 3 this is to calculate the angles between the atoms


def angle(m, n, s):
    """Returns the angle between three coordinates

    Parameters:
            three lists of atomic xyz coordinates
   
    Return:
            the value of the angle between the coordinates"""
    list_m = np.array(m)
    list_n = np.array(n)
    list_s = np.array(s)
    r1 = list_m - list_n
    r2 = list_m - list_s
    dot_product = np.dot(r1, r2)
    magnitude_product = np.linalg.norm(r1) * np.linalg.norm(r2)
    return(math.acos(dot_product / magnitude_product))

def allangles(dictionary):
    """Returns a dictionary of all the angles between the atoms

    Parameters:
            a dictionary of indexed atoms and their xyz coordinates

    Return:
            a dictionary of the atoms and the angles between them"""
    angles = {}
    for m, n, s in itertools.permutations(dictionary.keys(), 3):
        angles[m, n, s] = angle(dictionary[m], dictionary[n], dictionary[s])
    return(angles{})  

#part A 1 finding the out of plane angles


def outofplaneangle(m, n, s, t):
    """Returns the angle coming out of the plane for four atoms

    Parameters:
            four lists of atomic xyz coordinates

    Return:
            the value of the out of plane angle"""
    m_list = np.array(m)
    n_list = np.array(n)
    s_list = np.array(s)
    t_list = np.array(t)
    r_mn = m_list - n_list
    r_ms = m_list - s_list
    r_mt = m_list - t_list
    e_mn = r_mn / np.linalg.norm(r_mn)
    e_ms = r_ms / np.linalg.norm(r_ms)
    e_mt = r_mt / np.linalg.norm(r_mt)
    angle_v = angle(m, n, s)
    return(math.asin(np.dot(((np.cross(e_mn, e_ms)) / math.sin(angle_v)), e_mt)))

def alloutofplaneangles(dictionary):
    """Returns all of the out of plane angles for a set of atoms

    Parameters:
            a dictionary of indexed atomic symbols and their xyz coordinates

    Return:
            a dictionary of the atoms and the out of plane angle for each set"""
    outofplaneangles = {}
    for m, n, s, t in itertools.permutations(dictionary.keys(), 4):
        outofplaneangles[m, n, s, t] = outofplaneangle(dictionary[m], dictionary[n], dictionary[s], dictionary[t])
    return(outofplaneangles{})





