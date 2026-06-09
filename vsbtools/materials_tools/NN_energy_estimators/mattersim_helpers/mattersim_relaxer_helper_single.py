# mattersim_relaxer_helper_single.py
import sys
import io
import os

import numpy as np
import torch
from ase.io import read, write

# --------------------------------------------------------------------------- #
#  Redirect ALL normal stdout to stderr (imports + runtime)                   #
# --------------------------------------------------------------------------- #

_orig_stdout = sys.stdout   # real stdout – only used for final JSON


class _StdoutToStderr(io.TextIOBase):
    def write(self, s):
        sys.stderr.write(s)
    def flush(self):
        sys.stderr.flush()


# from here on, any plain print() goes to stderr
sys.stdout = _StdoutToStderr()

from mattersim.forcefield.potential import MatterSimCalculator
from mattersim.applications.relax import Relaxer

force_gpu_index = os.getenv("VSB_FORCE_GPU_INDEX")
if force_gpu_index is not None:
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = force_gpu_index
device = "cuda" if (force_gpu_index is not None or torch.cuda.is_available()) else "cpu"

# mattersim looks up the checkpoint in its own directory
calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth", device=device)
relaxer = Relaxer(
    optimizer="BFGS",
    filter="FrechetCellFilter",
    constrain_symmetry=False,
)

# --------------------------------------------------------------------------- #


def relax_one(ats):
    atoms = ats.copy()
    atoms.positions += 0.1 * np.random.randn(len(atoms), 3)
    atoms.calc = calc

    flag, result = relaxer.relax(atoms, steps=500)
    # print("relaxed" if flag else "relaxation failed")  # goes to stderr
    return result


def main():
    data = sys.stdin.read()
    atoms = read(io.StringIO(data), format="json")
    relaxed = relax_one(atoms)

    buf = io.StringIO()
    write(buf, relaxed, format="json")
    out = buf.getvalue()

    # Write JSON ONLY via the original stdout
    _orig_stdout.write(out)
    _orig_stdout.flush()


if __name__ == "__main__":
    main()
