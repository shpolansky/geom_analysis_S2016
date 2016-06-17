"""Module containing functions for 
geometry analysis of a set of coordinates
for any molecule."""

import itertools
import math
import numpy as np



def readgeom(filename):
    """Returns a dictionary of indexed elements with xyz coordinates

    Parameters:
            string filename -- name of .xyz file containing geometry to be read

    Return:
            dict coord_dict -- dictionary containing all atoms and xyz coordinates.
            Keys are in the format "SymbolNumber", where Symbol is the element symbol
            for the atom, and values are lists of the Cartesian coordinates, as floats
            in [x, y, z] format."""
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
    return(coord_dict)



def distance(m, n):
    """Returns the distance between two atomic coordinates

    Parameters:
            list m and n -- contains coordinates for atom M, N as floats in [x, y, z] format
   
    Return:
            the Cartesian distance value between the coordinates"""
    zdiff = m[2] - n[2]
    ydiff = m[1] - n[1]
    xdiff = m[0] - n[0]
    return(((xdiff ** 2) + (ydiff ** 2) + (zdiff ** 2)) ** 0.5)

def alldistances(dictionary):
    """Returns all of the distances from a dictionary of coordinates

    Parameters:
            dict dictionary -- a dictionary containing keys that are atomic
            symbols and values that are Cartesian coordinates as floats in [x, y, z] format

    Return:
            dict distances -- dictionary of the atoms and their distances.
            The keys are the atomic symbols of the two atoms "Symbol1", "Symbol2"
            and the values are the distance values between the atoms."""
    distances = {}
    for m, n in itertools.combinations(dictionary.keys(), 2):   
        distances[m, n] = distance(dictionary[m], dictionary[n])
    return(distances)



def angle(m, n, s):
    """Returns the angle between three coordinates

    Parameters:
            list m, n, s -- contains coordinates for atom M, N, and S
            as floats in [x, y, z] format
   
    Return:
            the Cartesian value of the angle between the coordinates"""
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
            dict dictionary -- dictionary of atoms and their Cartesian coordinates.
            The keys are element symbols, and the values are lists of Cartesian coordinates
            as floats in [x, y, z] format.

    Return:
            dict angles -- dictionary of the atoms and the angles between them. The keys
            are the atomic symbols in "Symbol1", "Symbol2", "Symbol3" format, and the values
            are the values of the angles between the atoms."""
    angles = {}
    for m, n, s in itertools.permutations(dictionary.keys(), 3):
        angles[m, n, s] = angle(dictionary[m], dictionary[n], dictionary[s])
    return(angles)  



def outofplaneangle(m, n, s, t):
    """Returns the angle coming out of the plane for four atoms

    Parameters:
            list m, n, s, t -- contains coordinates of atom M, N, S, and T in Cartesian
            coordinates as floats in [x, y, z] format

    Return:
            the Cartesian value of the out of plane angle"""
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
            dict dictionary -- dictionary of atoms and their Cartesian coordinates. The keys
            are the atomic symbols, and the values are lists of Cartesian coordinates as floats in 
            [x, y, z] format.

    Return:
            dict outofplaneangles -- dictionary of atoms and the out of plane angles between them.
            The keys are the atomic symbols in "Symbol1", "Symbol2", "Symbol3", "Symbol4", and the
            values are the values of the out of plane angles."""
    outofplaneangles = {}
    for m, n, s, t in itertools.permutations(dictionary.keys(), 4):
        outofplaneangles[m, n, s, t] = outofplaneangle(dictionary[m], dictionary[n], dictionary[s], dictionary[t])
    return(outofplaneangles)





