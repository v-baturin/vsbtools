#this program calculates coordinating nambers of cell with given lattice, coordinates and radii of atoms
#developer - Bushlanov Pavel
import numpy as np
from scipy.spatial.distance import cdist

def calculate_coordinating_numbers(radii,lattice,coordinates):
    sz = 27
    closest = np.array([[[0,0,0]],[[-1,0,0]],[[-1,0,-1]],[[-1,-1,-1]],[[-1,-1,0]],[[0,-1,0]],[[0,-1,-1]],[[0,0,-1]],[[-1,-1,1]],[[-1,0,1]],[[-1,1,1]],[[-1,1,0]],[[-1,1,-1]],[[0,-1,1]],[[0,0,1]],[[0,1,1]],[[0,1,0]],[[0,1,-1]],[[1,0,0]],[[1,0,-1]],[[1,-1,-1]],[[1,-1,0]],[[1,-1,1]],[[1,0,1]],[[1,1,1]],[[1,1,0]],[[1,1,-1]]])
    vertices = np.dot(np.concatenate((closest + coordinates),axis=0),lattice)
    radii_large = np.tile(radii,sz)
    dists = cdist(vertices,vertices)
    dists_no_self = np.delete(np.triu(dists,1),0,1) + np.delete(np.tril(dists,-1),dists.shape[1]-1,1)
    base_bond_legth = radii_large.reshape(1,radii_large.size)+radii_large.reshape(radii_large.size,1)
    base_bond_legth_no_self = np.delete(np.triu(base_bond_legth,1),0,1) + np.delete(np.tril(base_bond_legth,-1),base_bond_legth.shape[1]-1,1)
    order = np.exp(-(dists_no_self-base_bond_legth_no_self)/0.23)
    coord_number = order.sum(axis=1)/order.max(axis=1)
    return coord_number

