#!/usr/bin/env python3

import numpy as np
from ase.io import read

supercell_min_size = 10.

atoms = read('POSCAR')
cell = atoms.cell
supercell_shape = []
vol = np.abs(np.linalg.det(cell))
for i in range(3):
    j, k = np.roll(np.arange(2), i)[0:2]
    area = np.linalg.norm(np.cross(cell[j], cell[k]))
    supercell_shape.append(int(np.ceil(supercell_min_size / (vol/area))))

print(" ".join(map(str, supercell_shape)))




