from __future__ import annotations

import socket
from pathlib import Path
import subprocess
import io
import sys

from ase.io import write, read

# --------------------------------------------------------------------------- #
#  Resolve Python interpreter for mattersim by hostname                       #
# --------------------------------------------------------------------------- #

MATTERSIM_PYTHON_PATHS = {
    "nina":   "/home/vsbat/work/mattergen/mattergenbis/.venv/bin/python",
    "taurus": "/home/vsbat/work/venvs/mattersim_venv/bin/python",
    "serpens": "/home/vsbat/work/venvs/mattersim/bin/python",
}


def _resolve_mattersim_python() -> Path:
    host = socket.gethostname()
    if host in MATTERSIM_PYTHON_PATHS:
        return Path(MATTERSIM_PYTHON_PATHS[host])
    venv_root = input(
        "Enter full path to mattersim-containing virtual environment: "
    ).strip()
    return Path(venv_root) / "bin" / "python"


OTHER_PYTHON = _resolve_mattersim_python()

HERE = Path(__file__).resolve().parent
HELPER_SINGLE = HERE / "mattersim_relaxer_helper_single.py"
HELPER_BATCH = HERE / "mattersim_relaxer_helper_batch.py"


# --------------------------------------------------------------------------- #
#  Single-shot relaxation                                                     #
# --------------------------------------------------------------------------- #

def relax(atoms):
    """
    Relax a single ASE Atoms object using the external mattersim venv.
    """
    buf = io.StringIO()
    write(buf, atoms, format="json")
    json_atoms = buf.getvalue()

    proc = subprocess.run(
        [str(OTHER_PYTHON), "-u", str(HELPER_SINGLE)],
        input=json_atoms,
        text=True,
        capture_output=True,  # inspect stdout / stderr ourselves
    )

    if proc.returncode != 0:
        sys.stderr.write(
            f"[mattersim_relaxer] helper_single failed (code {proc.returncode})\n"
        )
        if proc.stdout:
            sys.stderr.write("--- helper stdout ---\n" + proc.stdout + "\n")
        if proc.stderr:
            sys.stderr.write("--- helper stderr ---\n" + proc.stderr + "\n")
        raise RuntimeError("mattersim_relaxer_helper_single.py failed")

    if not proc.stdout.strip():
        sys.stderr.write(
            "[mattersim_relaxer] helper_single returned success but stdout is empty.\n"
        )
        if proc.stderr:
            sys.stderr.write("--- helper stderr ---\n" + proc.stderr + "\n")
        raise RuntimeError("helper_single produced no output")

    return read(io.StringIO(proc.stdout), format="json")


# --------------------------------------------------------------------------- #
#  Streaming relaxation for many structures                                  #
# --------------------------------------------------------------------------- #

class RelaxStream:
    """
    Stream many Atoms through a single external mattersim process.

    Protocol:
      - stdin:  one ASE-JSON structure per line (newlines flattened to spaces)
      - stdout: one ASE-JSON relaxed structure per line
    """

    def __init__(self):
        self.proc = subprocess.Popen(
            [str(OTHER_PYTHON), "-u", str(HELPER_BATCH)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
        )

    def relax(self, atoms):
        # serialise to JSON and send a single line
        buf = io.StringIO()
        write(buf, atoms, format="json")
        json_str = buf.getvalue().replace("\n", " ")
        self.proc.stdin.write(json_str + "\n")
        self.proc.stdin.flush()

        # receive one line back
        while True:
            line = self.proc.stdout.readline()
            if not line:
                raise RuntimeError("mattersim_relaxer_helper_batch terminated")
            line = line.strip()
            if not line:
                continue
            if line.startswith("ERR "):
                raise RuntimeError(line)
            return read(io.StringIO(line), format="json")

    def close(self):
        if self.proc.stdin and not self.proc.stdin.closed:
            self.proc.stdin.close()
        self.proc.wait()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()


if __name__ == "__main__":
    from ase.build import bulk  # smoke test

    ats = bulk("Cu", "fcc", a=3.6)
    print("Relaxed structure (single):", relax(ats))

    atoms_batch = [
        bulk("Cu", "fcc", a=3.6),
        bulk("Cu", "fcc", a=3.6, orthorhombic=True),
        bulk("Cu", "fcc", a=3.6, cubic=True),
    ]

    with RelaxStream() as es:
        structures = [es.relax(at) for at in atoms_batch]

    print("Got", len(structures), "relaxed structures from batch helper.")
