import io, subprocess
from pathlib import Path
from ase.io import write


other_python = "/home/vsbat/work/venvs/mattersim_venv/bin/python"
energy_cli_single   = (Path(__file__).parent / "mattersim_helper_single.py").as_posix()
energy_worker = (Path(__file__).parent / "mattersim_helper_batch.py").as_posix()

def get_energy(atms, mattersim_python=None):
    # ---- serialise the structure to an in-memory JSON string -------------
    if mattersim_python is None: mattersim_python = other_python
    assert Path(mattersim_python).exists(), "Path to mattersim python virtual environment not found"
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
    def __init__(self, mattersim_python=None):
        if mattersim_python is None: mattersim_python = other_python
        self.other_python = mattersim_python
        assert Path(self.other_python).exists(), "Path to mattersim python virtual environment not found"
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
