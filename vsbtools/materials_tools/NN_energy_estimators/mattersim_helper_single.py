import torch
import sys, io
from ase.io import read
from mattersim.forcefield import MatterSimCalculator

device = "cuda" if torch.cuda.is_available() else "cpu"
"""
Read an ASE JSON structure from stdin, print the total energy to stdout.
"""

calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth", device=device)

atoms = read(io.StringIO(sys.stdin.read()), format="json")
atoms.calc = calc
e = atoms.get_potential_energy()
print(e)