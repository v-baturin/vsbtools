from __future__ import print_function
import sys
from randomTopology import generate_structure

try:
    import numpy as np
    import scipy.io as sio
except ImportError:
    print('<CALLRESULT>')
    print((u''.join('0 ' for i in range(13))))
    exit()


individual_number = int(sys.argv[1])#not used
lattice_volume = float(sys.argv[2])
symmetry_simbol = sys.argv[3]#not used
sz = int(sys.argv[4])
atom_numbers_by_type = np.array(sys.argv[5:5+sz], dtype = int)
#coordinationNumbers = np.array(sys.argv[5+sz:],np.int32).reshape(sz,2)

struct_name,lattice_vectors,coordinates = generate_structure(lattice_volume,atom_numbers_by_type)

np.set_printoptions(threshold=sys.maxsize)
print('<CALLRESULT>')

coordinates_str = np.array_str(lattice_vectors) + u' ' + np.array_str(coordinates)

name = struct_name[5:struct_name.find('-')]

print(name + ' ' + coordinates_str.translate({ord('['):None,ord(']'):None,ord('\n'):None}))
