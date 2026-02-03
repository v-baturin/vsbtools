# mattersim_relaxer_helper_batch.py
import sys
import io
import traceback

import numpy as np
import torch
from ase.io import read, write

# --------------------------------------------------------------------------- #
#  Redirect ALL normal stdout to stderr (imports + runtime)                   #
# --------------------------------------------------------------------------- #

_orig_stdout = sys.stdout   # real stdout – used only for protocol lines


class _StdoutToStderr(io.TextIOBase):
    def write(self, s):
        sys.stderr.write(s)
    def flush(self):
        sys.stderr.flush()


sys.stdout = _StdoutToStderr()

from mattersim.forcefield.potential import MatterSimCalculator
from mattersim.applications.relax import Relaxer

device = "cuda" if torch.cuda.is_available() else "cpu"

calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth", device=device)
relaxer = Relaxer(
    optimizer="BFGS",
    filter="ExpCellFilter",
    constrain_symmetry=True,
)

# --------------------------------------------------------------------------- #


def relax_one(ats):
    ats = ats.copy()
    ats.positions += 0.1 * np.random.randn(len(ats), 3)
    ats.calc = calc

    flag, result = relaxer.relax(ats, steps=500)
    if not flag:
        result = ats
    print("relaxed with mattersim" if flag else "relaxation failed, original structure returned")  # to stderr
    return result


for i, line in enumerate(sys.stdin):
    line = line.strip()
    if not line:
        continue

    try:
        atoms = read(io.StringIO(line), format="json")
        relaxed = relax_one(atoms)
        print(f"{i}")
        buf = io.StringIO()
        write(buf, relaxed, format="json")
        json_str = buf.getvalue().replace("\n", " ")  # one JSON per line

        _orig_stdout.write(json_str + "\n")
        _orig_stdout.flush()

    except Exception as err:
        # Signal error to caller via stdout, details to stderr
        _orig_stdout.write(f"ERR {err}\n")
        _orig_stdout.flush()
        traceback.print_exc(file=sys.stderr)
