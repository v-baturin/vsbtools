import torch
import sys, io
from ase.io import read
# --- temporary hijack of stdout so import-time banners go to stderr ----------
class _Stdout2Stderr(io.TextIOBase):
    def write(self, s):  sys.stderr.write(s)
    def flush(self):     sys.stderr.flush()

_orig_stdout = sys.stdout
sys.stdout   = _Stdout2Stderr()      # anything printed during imports → stderr
# ---------------------------------------------------------
from mattersim.forcefield import MatterSimCalculator
# restore clean stdout for our own protocol
sys.stdout = _OrigStdout = _orig_stdout

# energy_worker.py  –  read JSON atoms on stdin, emit energy per line
device = "cuda" if torch.cuda.is_available() else "cpu"
calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth", device=device)

for line in sys.stdin:                     # newline-delimited protocol
    line = line.strip()
    if not line:                           # ignore blank lines
        continue
    try:
        atoms = read(io.StringIO(line), format="json")
        atoms.calc = calc                  # attach the calculator
        e = atoms.get_potential_energy()  # get the energy
        print(e, flush=True)               # one float, one line
    except Exception as err:
        print(f"ERR {err}", flush=True)