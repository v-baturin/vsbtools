import sys
import io, subprocess
import os
from pathlib import Path
from ase.io import read, write
from ..external_paths import (
    ExternalPathResolutionError,
    python_executable_from_venv,
    python_import_validator,
    resolve_external_path,
)

def _resolve_mattersim_python(*, prompt: bool, explicit_path=None) -> Path | None:
    return resolve_external_path(
        name="MatterSim Python environment",
        config_key="mattersim_python",
        env_var="MATTERSIM_PYTHON",
        explicit_path=explicit_path,
        normalizer=python_executable_from_venv,
        validator=python_import_validator("mattersim.forcefield"),
        prompt=prompt,
        required=prompt,
        prompt_text="Enter full path to MatterSim virtual environment or Python executable: ",
    )


other_python = _resolve_mattersim_python(prompt=False)


def _require_mattersim_python(mattersim_python=None) -> Path:
    global other_python
    if mattersim_python is not None:
        resolved = _resolve_mattersim_python(prompt=False, explicit_path=mattersim_python)
    else:
        resolved = other_python
        if resolved is None:
            resolved = _resolve_mattersim_python(prompt=True)
            other_python = resolved
    if resolved is None:
        raise ExternalPathResolutionError(
            "MatterSim Python environment is not configured. Set MATTERSIM_PYTHON "
            "or configure external_paths.json."
        )
    return Path(resolved)

HERE = Path(__file__).resolve().parent
helpers_path = HERE / "mattersim_helpers"
energy_cli_single   = (helpers_path / "mattersim_helper_single.py").as_posix()
energy_worker = (helpers_path / "mattersim_helper_batch.py").as_posix()
relax_cli_single = helpers_path / "mattersim_relaxer_helper_single.py"
relax_worker = helpers_path / "mattersim_relaxer_helper_batch.py"



def _helper_env(force_gpu: int | None = None) -> dict[str, str]:
    env = os.environ.copy()
    gpu_index = _normalize_force_gpu(force_gpu)
    if gpu_index is None:
        env.pop("VSB_FORCE_GPU_INDEX", None)
    else:
        env["VSB_FORCE_GPU_INDEX"] = str(gpu_index)
    return env


def _normalize_force_gpu(force_gpu: int | None) -> int | None:
    if force_gpu is None:
        return None
    if isinstance(force_gpu, bool):
        raise TypeError("force_gpu must be int | None, bool is not supported")
    idx = int(force_gpu)
    if idx < 0:
        raise ValueError(f"force_gpu must be >= 0 or None, got {force_gpu}")
    return idx


def get_energy(atms, mattersim_python=None, force_gpu: int | None = None):
    # ---- serialise the structure to an in-memory JSON string -------------
    mattersim_python = _require_mattersim_python(mattersim_python)
    buf = io.StringIO()
    write(buf, atms, format="json")
    json_atoms = buf.getvalue()

    # ---- call the calculator python and feed it via stdin ----------------
    proc = subprocess.run(
        [mattersim_python.as_posix(), energy_cli_single],
        input=json_atoms,
        text=True,
        capture_output=True,
        check=True,
        env=_helper_env(force_gpu=force_gpu),
    )

    return float(proc.stdout.strip().split('\n')[-1])      # one number per call

class EnergyStream:
    """Context-manager that streams Atoms objects to the other v-env."""
    def __init__(self, mattersim_python=None, force_gpu: int | None = None):
        self.other_python = _require_mattersim_python(mattersim_python)
        self.proc = subprocess.Popen(
            [self.other_python.as_posix(), energy_worker],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            env=_helper_env(force_gpu=force_gpu),
        )

    def calc(self, atoms):
        if self.proc.poll() is not None:
            raise RuntimeError(f"MatterSim energy worker terminated with code {self.proc.returncode}")

        # serialise to a JSON line
        buf = io.StringIO()
        write(buf, atoms, format="json")
        json_line = buf.getvalue().replace("\n", " ")  # keep it single-line
        try:
            self.proc.stdin.write(json_line + "\n")
            self.proc.stdin.flush()
        except BrokenPipeError as err:
            raise RuntimeError(f"MatterSim energy worker terminated with code {self.proc.poll()}") from err

        # receive the answer
        out = self.proc.stdout.readline()
        if out == "":
            raise RuntimeError(f"MatterSim energy worker closed stdout with code {self.proc.poll()}")
        out = out.strip()
        if out.startswith("ERR"):
            raise RuntimeError(out)
        # print(out)
        return out

    def close(self):
        if self.proc.stdin and not self.proc.stdin.closed:
            self.proc.stdin.close()
        if self.proc.stdout and not self.proc.stdout.closed:
            self.proc.stdout.close()
        self.proc.wait()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

def relax(ats, force_gpu: int | None = None):
    """
    Relax a single ASE Atoms object using the external mattersim venv.
    """
    mattersim_python = _require_mattersim_python()
    buf = io.StringIO()
    write(buf, ats, format="json")
    json_atoms = buf.getvalue()

    proc = subprocess.run(
        [mattersim_python.as_posix(), "-u", str(relax_cli_single)],
        input=json_atoms,
        text=True,
        capture_output=True,  # inspect stdout / stderr ourselves
        env=_helper_env(force_gpu=force_gpu),
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

    def __init__(self, force_gpu: int | None = None):
        mattersim_python = _require_mattersim_python()
        self.proc = subprocess.Popen(
            [mattersim_python.as_posix(), "-u", str(relax_worker)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            env=_helper_env(force_gpu=force_gpu),
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
