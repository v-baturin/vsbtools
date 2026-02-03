import sys, io
from ase.io import read
from tensorpotential.calculator import grace_fm
calc = grace_fm('GRACE-2L-OMAT')

"""
Read an ASE JSON structure from stdin, print the total energy to stdout.
"""
atoms = read(io.StringIO(sys.stdin.read()), format="json")
atoms.calc = calc
e = atoms.get_potential_energy()
print(e)