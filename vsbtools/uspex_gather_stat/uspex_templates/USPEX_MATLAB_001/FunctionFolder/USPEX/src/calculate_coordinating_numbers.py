#this is simple adapter to proper module coord_number.py
from __future__ import print_function
import sys
import string
try:
    import numpy as np
except ImportError:
    print('<CALLRESULT>')
    print(0)
    exit()
import coord_number

type_number = int(sys.argv[1])
numIons = np.array(sys.argv[2:2+type_number],np.int32)
R_val = np.array(sys.argv[2+type_number:2+2*type_number],np.float64)
lattice = np.array(sys.argv[2+2*type_number:11+2*type_number],np.float64).reshape((3,3))
coordinates = np.array(sys.argv[11+2*type_number:],np.float64).reshape(numIons.sum(),3)
radii = np.repeat(R_val,numIons)
cn = coord_number.calculate_coordinating_numbers(radii,lattice,coordinates)

np.set_printoptions(threshold=np.nan)
print('<CALLRESULT>')
print((np.array_str(cn)+u' ').translate({ord('['):None,ord(']'):None,ord('\n'):None}))
