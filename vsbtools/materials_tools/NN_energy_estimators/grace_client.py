import sys
import socket
import io, subprocess
from pathlib import Path
from ase.io import read, write

host = socket.gethostname()
GRACE_PYTHON_PATHS = {'nina': "/home/vsbat/work/python_venvs/grace_potential/bin/python",
                          # "taurus": "/home/vsbat/work/venvs/GRACE_venv/bin/python",
                      "serpens": "/home/vsbat/work/venvs/grace/bin/python"
                          }

if host not in GRACE_PYTHON_PATHS:
    other_python = Path(input("Enter full path to grace-containing virtual environment")) / "bin/python"
else:
    other_python = Path(GRACE_PYTHON_PATHS[host])

HERE = Path(__file__).resolve().parent
helpers_path = HERE / "grace_helpers"
energy_cli_single   = (helpers_path / "grace_helper_single.py").as_posix()
energy_worker = (helpers_path / "grace_helper_batch.py").as_posix()
relax_cli_single = helpers_path / "grace_relaxer_helper_single.py"
relax_worker = helpers_path / "grace_relaxer_helper_batch.py"



def get_energy(atms, grace_python=None):
    # ---- serialise the structure to an in-memory JSON string -------------
    if grace_python is None:
        grace_python = other_python
    assert Path(grace_python).exists(), "Path to grace python virtual environment not found"
    buf = io.StringIO()
    write(buf, atms, format="json")
    json_atoms = buf.getvalue()

    # ---- call the calculator python and feed it via stdin ----------------
    proc = subprocess.run(
        [other_python, energy_cli_single],
        input=json_atoms,
        text=True,
        capture_output=True,
        check=True
    )

    return float(proc.stdout.strip().split('\n')[-1])      # one number per call

class EnergyStream:
    """Context-manager that streams Atoms objects to the other v-env."""
    def __init__(self, grace_python=None):
        if grace_python is None: grace_python = other_python
        self.other_python = grace_python
        assert Path(self.other_python).exists(), "Path to grace python virtual environment not found"
        self.proc = subprocess.Popen(
            [self.other_python, energy_worker],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
        )

    def calc(self, atoms):
        # serialise to a JSON line
        buf = io.StringIO()
        write(buf, atoms, format="json")
        json_line = buf.getvalue().replace("\n", " ")  # keep it single-line
        self.proc.stdin.write(json_line + "\n")
        self.proc.stdin.flush()

        # receive the answer
        out = self.proc.stdout.readline().strip()
        if out.startswith("ERR"):
            raise RuntimeError(out)
        # print(out)
        return out

    def close(self):
        self.proc.stdin.close()
        self.proc.wait()

    # optional: make it a context-manager
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc, tb):
        self.close()

def relax(ats):
    """
    Relax a single ASE Atoms object using the external grace venv.
    """
    buf = io.StringIO()
    write(buf, ats, format="json")
    json_atoms = buf.getvalue()

    proc = subprocess.run(
        [str(other_python), "-u", str(relax_cli_single)],
        input=json_atoms,
        text=True,
        capture_output=True,  # inspect stdout / stderr ourselves
    )

    if proc.returncode != 0:
        sys.stderr.write(
            f"[grace_relaxer] helper_single failed (code {proc.returncode})\n"
        )
        if proc.stdout:
            sys.stderr.write("--- helper stdout ---\n" + proc.stdout + "\n")
        if proc.stderr:
            sys.stderr.write("--- helper stderr ---\n" + proc.stderr + "\n")
        raise RuntimeError("grace_relaxer_helper_single.py failed")

    if not proc.stdout.strip():
        sys.stderr.write(
            "[grace_relaxer] helper_single returned success but stdout is empty.\n"
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
    Stream many Atoms through a single external grace process.

    Protocol:
      - stdin:  one ASE-JSON structure per line (newlines flattened to spaces)
      - stdout: one ASE-JSON relaxed structure per line
    """

    def __init__(self):
        self.proc = subprocess.Popen(
            [str(other_python), "-u", str(relax_worker)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
        )

    def relax(self, ats):
        # serialise to JSON and send a single line
        buf = io.StringIO()
        write(buf, ats, format="json")
        json_str = buf.getvalue().replace("\n", " ")
        self.proc.stdin.write(json_str + "\n")
        self.proc.stdin.flush()

        # receive one line back
        while True:
            line = self.proc.stdout.readline()
            if not line:
                raise RuntimeError("grace_relaxer_helper_batch terminated")
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



# ------------------ example usage ------------------------------------------
if __name__ == "__main__":
    from ase.build import molecule, bulk
    atoms = molecule("H2O")
    print("Energy =", get_energy(atoms), "eV")

    atoms_batch = [bulk('Cu', 'fcc', a=3.6),
                   bulk('Cu', 'fcc', a=3.6, orthorhombic=True),
                   bulk('Cu', 'fcc', a=3.6, cubic=True)]  # pretend = 10 000

    with EnergyStream() as es:
        energies = [es.calc(at) for at in atoms_batch]

    print("Energies:", energies)

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
