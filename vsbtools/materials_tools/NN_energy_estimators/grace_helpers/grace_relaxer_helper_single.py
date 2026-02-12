# mattersim_relaxer_helper_single.py
import sys
import io
import numpy as np
from ase.io import read, write
from ase.optimize import BFGS
from ase.filters import EXpCellFilter

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
from tensorpotential.calculator import grace_fm
calc = grace_fm('GRACE-2L-OMAT-large-ft-AM')



# --------------------------------------------------------------------------- #


def relax_one(ats):
    atoms = ats.copy()
    atoms.positions += 0.1 * np.random.randn(len(ats), 3)
    atoms.calc = calc

    atoms_ucf = EXpCellFilter(atoms)
    dyn = BFGS(atoms_ucf, logfile="relax.log")
    dyn.run(fmax=0.02, steps=500)

    flag = dyn.converged()
    if not flag:
        result = ats
    else:
        result = atoms
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
