import torch
import sys, io
import os
from ase.io import read
from mattersim.forcefield import MatterSimCalculator

force_gpu_index = os.getenv("VSB_FORCE_GPU_INDEX")
if force_gpu_index is not None:
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = force_gpu_index
device = "cuda" if (force_gpu_index is not None or torch.cuda.is_available()) else "cpu"
"""
Read an ASE JSON structure from stdin, print the total energy to stdout.
"""

calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth", device=device)

atoms = read(io.StringIO(sys.stdin.read()), format="json")
atoms.calc = calc
e = atoms.get_potential_energy()
print(e)
