from __future__ import division

import itertools
import numpy as np
import os
import random
import scipy.io as sio
import sqlite3
import string

from math import sin, cos, radians, sqrt
from latConverter import latConverter

def check1DArray(param, type):
    assert param is not None and type
    return isinstance(param, np.ndarray) and len(param.shape) == 1 and param.shape[0] and isinstance(param[0], type)

def check2DArray(param, type):
    assert param is not None and type
    return isinstance(param, np.ndarray) and len(param.shape) == 2 and param.shape[0] and param.shape[1] and isinstance(param[0][0], type)

def partitions(set_):
    if not set_:
        yield []
        return
    for i in range(2**len(set_) // 2):
        parts = [set(), set()]
        for item in set_:
            parts[i&1].add(item)
            i >>= 1
        for b in partitions(parts[1]):
            yield [parts[0]] + b

def get_atom_indices_by_type(atom_type_permutation, nod_numbers_by_type, nod_partition):
    '''

    :param atom_type_permutation:
    :param nod_numbers_by_type:
    :param nod_partition:
    :return:
    '''
    atomIndices = []
    for i in np.argsort(atom_type_permutation):
        if i < len(nod_partition):
            indices = []
            for nod_type in nod_partition[i]:
                begin = np.sum(nod_numbers_by_type[:nod_type])
                indices.extend(range(begin, begin + nod_numbers_by_type[nod_type]))
            atomIndices.append(indices)
        else:
            atomIndices.append([])
    return atomIndices


class RandomCoordinatesFillFailed(BaseException):
    pass

def fillCoordinates(atom2Numbers, slots_by_type, bondSizes):
    assigned_bonds = []
    for number, slots in zip(atom2Numbers, slots_by_type):
        while number > 0:
            large_bond_args = np.nonzero(bondSizes >= slots)[0]
            if not len(large_bond_args):
                raise RandomCoordinatesFillFailed
            assigned_bond_arg = np.random.choice(large_bond_args)
            assigned_bonds.append(assigned_bond_arg)
            bondSizes[assigned_bond_arg] -= slots
            number -= 1
    return assigned_bonds

def match_slots(atoms, slots, bonds):
    assert check1DArray(atoms, int)
    assert check1DArray(slots, int)
    assert check1DArray(bonds, int)

    ind = np.argsort(slots)[::-1]
    slots = slots[ind]
    atoms = atoms[ind]
    bonds = np.sort(bonds)

    while len(atoms):
        a = atoms[0]
        s = slots[0]
        if np.sum(bonds // s) < a:
            return False
        elif np.sum(bonds//s) == a:
            bonds = bonds % s
            bonds = np.sort(bonds[np.nonzero(bonds)])
            atoms = atoms[1:]
            slots = slots[1:]
        else:
            b = bonds[0]
            while b > 0:
                for i in range(len(atoms)):
                    rem = np.array([b], dtype=int)
                    for s,a in zip(slots[i:],atoms[i:]):
                        rem = np.concatenate([rem,[max(rem[-1]%s, rem[-1] - s*a)]])
                    if rem[-1] == 0 or np.all((rem[:-1]-rem[1:]) // slots[i:] == atoms[i:]):
                        assigned = (rem[:-1] - rem[1:]) // slots[i:]
                        atoms[i:] = atoms[i:] - assigned
                        bonds[0] = b - np.sum(assigned * slots[i:])
                        b = 0
                        break
                b -= 1
            bonds = bonds[1:]
            ind = np.nonzero(atoms)
            atoms = atoms[ind]
            slots = slots[ind]
            ind = np.argsort(slots)[::-1]
            slots = slots[ind]
            atoms = atoms[ind]
    return True

class Fill2AtomsTimeout(BaseException):
    pass

class Multicell:
    '''

    '''

    def __init__(self, multiplicity, latticeVolume, cell, bonds, nodeCoordinates):
        '''

        :param multiplicity:
        :param latticeVolume:
        :param cell:
        :param bonds:
        :param nod_coordinates:
        :return:
        '''

        assert isinstance(multiplicity, int) and multiplicity in [1,2,4,8]
        assert isinstance(latticeVolume, float) and latticeVolume > 0
        assert check1DArray(cell, float)
        assert check2DArray(bonds, int)
        assert check2DArray(nodeCoordinates, float)

        self.multiplicity = multiplicity
        self.perm = np.random.permutation(np.arange(3))
        m = self.multiplicity - 1
        self.mults = np.dot(np.array([m&1,m&2,m&4]), np.identity(3)[self.perm]) + np.ones(3)

        cell[3] = radians(cell[3])
        cell[4] = radians(cell[4])
        cell[5] = radians(cell[5])

        vnf = self.mults.reshape([3,1]) * latConverter([cell])

        factor = np.power(latticeVolume / np.abs(np.linalg.det(vnf)), 1/3)
        self.latticeVectors = factor * vnf
        self.nodeCoordinates = nodeCoordinates
        bondVectorsBase = nodeCoordinates[bonds[1]] + bonds[2:].T - nodeCoordinates[bonds[0]]
        bond_origin_base = nodeCoordinates[bonds[0]]
        shift = np.identity(3)[self.perm]
        self.bond_vectors = np.tile(bondVectorsBase, (self.multiplicity,1))
        self.bond_origin = np.empty([0,3])
        for i in range(self.multiplicity):
            self.bond_origin = np.vstack((self.bond_origin, bond_origin_base + (i&1)*shift[0] + (i&2)*shift[1] + (i&4)*shift[2]))
        self.bond_lengths = np.sqrt(np.sum(np.dot(self.bond_vectors / self.mults, self.latticeVectors)**2, axis=1))

    def match_slots(self, atom_2_numbers_by_type, atom_2_slots_by_type, slot_size):
        return match_slots(atom_2_numbers_by_type, atom_2_slots_by_type, np.array(self.bond_lengths / slot_size, dtype=int))

    def coordinates_3(self, atom_indices_by_type):
        '''
        :param atom_indices_by_type:
        :return: coordinates of >=3-coordinated atoms
        '''

        shift = np.identity(3)[self.perm]
        coordinates = []
        offset_2_coordinated = 0
        for atom_type in range(len(atom_indices_by_type)):
            for i in range(self.multiplicity):
                coordinates.extend(self.nodeCoordinates[atom_indices_by_type[atom_type]] + (i&1)*shift[0] + (i&2)*shift[1] + (i&4)*shift[2])
        return np.modf(np.array(coordinates / self.mults) + np.array([3.0001,3.0001,3.0001]))[0]

    def coordinates_2(self, atom2Numbers, atom2Slots, slot_size):
        '''

        :param atom2Numbers:
        :param atom_2_slots_by_type:
        :param slot_size:
        :return:
        '''

        count = 0
        while count < 1000:
            try:
                assigned_bonds = fillCoordinates(atom2Numbers, atom2Slots, self.bond_lengths//slot_size)
            except RandomCoordinatesFillFailed:
                count +=1
            else:
                coordinates = []
                for bond_ind, atom_per_bond in zip(*np.unique(assigned_bonds,return_counts=True)):
                    coordinates.extend(self.bond_origin[bond_ind] + self.bond_vectors[bond_ind] / (atom_per_bond+1)*i for i in range(1,atom_per_bond+1))
                return np.array(coordinates/self.mults) if coordinates else np.empty([0,3])
        raise Fill2AtomsTimeout


#def match_coordination_numbers(atom_coordination_numbers,nod_numbers_by_type,nod_coordination_numbers,atom_indices_by_type):
#    nod_coordination_numbers_all = np.repeat(nod_coordination_numbers,nod_numbers_by_type)
#    for atom_type in range(len(atom_coordination_numbers)):
#        match_vector = np.array(list((nod_coordination_numbers_all[atom_index] in atom_coordination_numbers[atom_type]) for atom_index in atom_indices_by_type[atom_type]))
#        if (len(atom_coordination_numbers[atom_type]) != 0) and (not (np.all(match_vector) )):
#            return False
#    return True


class CanNotMatch3Atoms(BaseException):
    pass


class CanNotMatch2Atoms(BaseException):
    pass


class TopologicalStructure:
    '''
    Topological Structure is an object which relates an entry
    from the topological database to requested cristal structure.
    '''

    def __init__(self, atom_3_numbers_by_type, atom_2_numbers_by_type, latticeVolume, atom_2_slots_by_type, slot_size, \
                 name, size, nnbt, ncn, bonds_number, symmetry, cell_info, ncoord, bnds):
        '''

        :param atom_3_numbers_by_type:
        :param atom_2_numbers_by_type:
        :param latticeVolume:
        :param atom_2_slots_by_type:
        :param slot_size:
        :param name: structure name in database
        :param size: total number of nodes in structure
        :param nnbt: number of nodes by node type
        :param ncn: coordination numbers of nodes by node type
        :param bonds_number: number of bonds
        :param symmetry:
        :param cell_info: latice parameters in 3x3 form
        :param ncoord: list of node coordinates
        :param bnds: list of bonds
        :return:
        '''

        assert check1DArray(atom_3_numbers_by_type, int)
        assert check1DArray(atom_2_numbers_by_type, int)
        assert check1DArray(atom_2_slots_by_type, int)
        assert isinstance(slot_size, int) and slot_size > 0

        fixed_atom_number = atom_3_numbers_by_type.sum()
        multiplicity = fixed_atom_number // size
        if multiplicity not in [1,2,4,8] or fixed_atom_number % size:
            raise CanNotMatch3Atoms
        self.name = name
        self.atom_3_numbers_by_type = atom_3_numbers_by_type
        self.atom_2_numbers_by_type = atom_2_numbers_by_type
        self.atom_2_slots_by_type = atom_2_slots_by_type
        self.slot_size = slot_size
        self.nod_numbers_by_type = np.array(nnbt.translate({ord('['):None,ord(']'):None}).split(), dtype=int)
        self.appropriateNodePartitions = self.getAppropriateNodePartitions(atom_3_numbers_by_type, multiplicity * self.nod_numbers_by_type)
        if not self.appropriateNodePartitions:
            raise CanNotMatch3Atoms
        bonds = np.array([line.split() for line in bnds.translate({ord('['):u' ',ord(']'):u' '}).split('\n')], dtype=int).T
        nodeCoordinates = np.array([line.split() for line in ncoord.translate({ord('['):u' ',ord(']'):u' '}).split('\n')], dtype=float)

        self.multicell = Multicell(multiplicity, latticeVolume, np.array(cell_info.split(), dtype=float), bonds, nodeCoordinates)
        if not self.multicell.match_slots(atom_2_numbers_by_type, atom_2_slots_by_type, slot_size):
            raise CanNotMatch2Atoms

    @staticmethod
    def getAppropriateNodePartitions(atomNumbersByType, nodeNumbersByType):
        assert check1DArray(atomNumbersByType, int)
        assert check1DArray(nodeNumbersByType, int)

        appropriateNodePartitions = []
        if nodeNumbersByType.shape[0] < 13: #partitions number is combinatoricaly large, need to cut them off
            nod_types = range(nodeNumbersByType.shape[0])
            for nodePartition in partitions(nod_types):
                if len(nodePartition) <= atomNumbersByType.shape[0]:
                    nodeNumbersForNodePartition = np.fromiter((np.sum(nodeNumbersByType[list(nodeTypeGroup)]) for nodeTypeGroup in nodePartition), dtype=int)
                    if np.all(np.sort(nodeNumbersForNodePartition)[::-1] <= np.sort(atomNumbersByType)[::-1][:nodeNumbersForNodePartition.shape[0]]):
                        appropriateNodePartitions.append(nodePartition)
        return appropriateNodePartitions

    def populate(self):
        for nod_partition_index in np.random.permutation(np.arange(len(self.appropriateNodePartitions))):
            nod_partition = self.appropriateNodePartitions[nod_partition_index]
            nodeNumbersForNodePartition = np.fromiter((np.sum(self.nod_numbers_by_type[list(nod_type_group)]) for nod_type_group in nod_partition), dtype=int)
            for atom_type_permutation in itertools.permutations(np.arange(self.atom_3_numbers_by_type.shape[0])):
                if np.all(self.multicell.multiplicity * nodeNumbersForNodePartition <= \
                            self.atom_3_numbers_by_type[list(atom_type_permutation)][:nodeNumbersForNodePartition.shape[0]]):
                    atom_indices_by_type = get_atom_indices_by_type(atom_type_permutation, self.nod_numbers_by_type, nod_partition)
                    coordinates = np.concatenate([self.multicell.coordinates_3(atom_indices_by_type),\
                                                  self.multicell.coordinates_2(self.atom_2_numbers_by_type, self.atom_2_slots_by_type, self.slot_size)])
                    #struct_name = name[5:name.find('-')]
                    return self.name, self.multicell.latticeVectors, coordinates


def getStructs(dbPath):
    assert isinstance(dbPath, str) and os.path.exists(dbPath)
    conn = sqlite3.connect(dbPath)
    c = conn.cursor()
    structs = list(c.execute("SELECT * FROM idealnets")) # WHERE size<=?",(total_atom_number,)))
    conn.close()
    return structs

def splitAtomNumbers(atomNumbersByType):
    atom_2_numbers_by_type = np.random.randint(3, size=atomNumbersByType.shape[0])
    atom_3_numbers_by_type = atomNumbersByType - atom_2_numbers_by_type
    return atom_3_numbers_by_type, atom_2_numbers_by_type


def generate_structure(latticeVolume, atomNumbersByType):
    assert isinstance(latticeVolume, float) and latticeVolume >= 0
    assert isinstance(atomNumbersByType, np.ndarray) and atomNumbersByType.shape[0] > 0
    assert isinstance(atomNumbersByType[0], int)

    atom_3_numbers_by_type = np.copy(atomNumbersByType)
    atom_2_numbers_by_type = atom_3_numbers_by_type - atomNumbersByType

    atom_2_slots_by_type = np.ones(atomNumbersByType.shape[0], dtype=int)
    slot_size = 1

    print(os.getcwd())
    structs = getStructs('idealnets.db')

    count = 0
    while count < 1000:
        if count > 500:
            atom_3_numbers_by_type, atom_2_numbers_by_type = splitAtomNumbers(atomNumbersByType)

        try:
            topologicalStructure = TopologicalStructure(atom_3_numbers_by_type, atom_2_numbers_by_type, \
                                                        latticeVolume, atom_2_slots_by_type, slot_size,\
                                                        *random.choice(structs))
        except CanNotMatch3Atoms:
            count += 1
        except CanNotMatch2Atoms:
            count += 1
        else:
            return topologicalStructure.populate()

    return u'data_0-', np.zeros(9), np.zeros(3 * np.sum(atomNumbersByType))
