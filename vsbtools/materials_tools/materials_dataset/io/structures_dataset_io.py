from __future__ import annotations

"""StructureDatasetIO – read/write utilities for CrystalDataset

Features
--------
* **Load** structures from a directory tree that may contain nested zip archives.
  Archive exploration is optional via *expand_archives* (default **True**).
* **Load** structures from a multi‑image POSCARS file.
* **Dump** individual POSCAR files for each entry.
* **Dump** a combined multi‑image POSCARS file (plus optional plain‑text ID list).

Zip handling (copy‑then‑explode with directory‑hierarchy preservation, akin to
`find … -name '*.zip' -exec cp {} --parents …`) is isolated in
`zip_handling.py` so that archive specifics stay out of this module.
"""

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Sequence, Collection

from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp import Poscar

from ..crystal_dataset import CrystalDataset, CrystalEntry
from .zip_handling import exploded_zip_tree         # archive helper
from ...ext_software_io.ext_software_io import read_poscars               # existing helper

LOG = logging.getLogger(__name__)

###############################################################################
# Internal helpers                                                           #
###############################################################################

def _safe_structure_from_file(path: Path, merge_tol: float = 0.01) -> Structure | None:
    """Return Structure or *None* on parser failure; logs the error."""
    try:
        return Structure.from_file(path, merge_tol=merge_tol)
    except Exception as exc:  # ValueError, RuntimeError, etc.
        LOG.warning("Skipping %s – %s", path, exc)
        return None


def _entry_id_sequence(prefix: str) -> Iterator[str]:
    """Yield id strings prefix0, prefix1, … endlessly."""
    counter = 0
    while True:
        yield f"{prefix}{counter}"
        counter += 1


def _iter_structure_files(root: Path, patterns: Sequence[str]) -> Iterator[Path]:
    """Yield files in *root* matching any of *patterns* (via rglob)."""
    for pattern in patterns:
        yield from root.rglob(pattern)

def _parse_extxyz_structures(extxyz_file: Path) -> list[Structure]:
    lines = extxyz_file.read_text().splitlines()
    idx = 0
    structures: list[Structure] = []
    lattice_re = re.compile(r'Lattice="([^"]+)"')

    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            continue
        try:
            natoms = int(line)
        except ValueError as exc:
            raise ValueError(f"Invalid atom count line in {extxyz_file}: {line}") from exc
        if idx + 1 >= len(lines):
            raise ValueError(f"Missing Lattice line after atom count in {extxyz_file}")

        lattice_line = lines[idx + 1]
        match = lattice_re.search(lattice_line)
        if not match:
            raise ValueError(f"Missing Lattice definition in {extxyz_file}: {lattice_line}")
        lattice_vals = [float(x) for x in match.group(1).split()]
        if len(lattice_vals) != 9:
            raise ValueError(f"Expected 9 lattice values in {extxyz_file}, got {len(lattice_vals)}")
        lattice = Lattice([
            lattice_vals[0:3],
            lattice_vals[3:6],
            lattice_vals[6:9],
        ])

        start = idx + 2
        end = start + natoms
        if end > len(lines):
            raise ValueError(f"Not enough atom lines in {extxyz_file} for natoms={natoms}")
        species: list[str] = []
        coords: list[list[float]] = []
        for atom_line in lines[start:end]:
            parts = atom_line.split()
            if len(parts) < 4:
                raise ValueError(f"Invalid atom line in {extxyz_file}: {atom_line}")
            species.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

        structures.append(Structure(lattice, species, coords, coords_are_cartesian=True))
        idx = end

    return structures

def get_batch_metadata(root: Path, prov_file):
    prov_metadata = None
    for found_provdata in root.rglob(prov_file):
        with open(found_provdata, 'rt') as pmdf:
            prov_metadata = pmdf.read()
            return prov_metadata

###############################################################################
# Public façade                                                              #
###############################################################################

@dataclass(slots=True)
class StructureDatasetIO:
    """Read and write `CrystalDataset` objects.

    Parameters
    ----------
    root : Path
        Source directory that may contain structure files and/or zip archives.
    patterns : Sequence[str], default ("*POSCAR*", "*.cif")
        Glob patterns that identify structure files **both** outside and inside
        archives.  Provide as many patterns as needed for your workflow.
    id_prefix : str, default ""
        Prefix for auto‑generated entry IDs.
    expand_archives : bool, default True
        If *True* the loader copies every ``*.zip`` (preserving directory tree)
        to a temporary area, recursively explodes them, and reads structures
        from that expanded tree **in addition to** structures already present
        on disk.  If *False* only files on disk are inspected.
    """

    root: Path
    patterns: Sequence[str] = ("*POSCAR*", "*.cif")
    id_prefix: str = ""
    expand_archives: bool = True
    source_name: str = "NA"
    batch_metadata_file: str | Path | None = None
    _id_iter: Iterator[str] = field(init=False, repr=False)

    # --------------------------------------------------------------------- #
    # Initialisation                                                        #
    # --------------------------------------------------------------------- #
    def __post_init__(self) -> None:
        if not self.root.is_dir():
            raise FileNotFoundError(self.root)
        self._id_iter = _entry_id_sequence(self.id_prefix)

    # --------------------------------------------------------------------- #
    # Reader methods                                                        #
    # --------------------------------------------------------------------- #
    def load_from_directory(
        self,
        *,
        elements: Collection[str] | None = None,
        message: str | None = None,
        # expand_archives: bool = True
    ) -> CrystalDataset:
        """Return a dataset built from *root*.

        Parameters
        ----------
        message : str | None
            Optional dataset message.
        expand_archives : bool | None
            Overrides the instance‑level ``expand_archives`` flag.
        """
        # use_archives = self.expand_archives if expand_archives is None else expand_archives

        entries: list[CrystalEntry] = []

        # 1) structures already on disk (outside archives)
        total = 0
        bad = 0
        for file in _iter_structure_files(self.root, self.patterns):
            batch_metadata = get_batch_metadata(self.root, self.batch_metadata_file) if \
                self.batch_metadata_file else None
            struct = _safe_structure_from_file(file)
            total += 1
            if struct:
                entries.append(CrystalEntry(id=next(self._id_iter), structure=struct,
                                            metadata={"source": self.source_name, "file": file.name}))
            else:
                bad += 1

        if elements:
            entries = [e for e in entries if not (set(e.composition.as_data_dict()["elements"]) - set(elements))]
        msg = message or f"{total - bad} structures out of {total} files collected from {self.root}"
        return CrystalDataset(entries, message=msg, supplementary_metadata={'batch_metadata': batch_metadata})

    def load_from_multiimage_poscar(
        self,
        poscars_file: Path,
        entries_ids_file: Path | None = None,
        *,
        message: str | None = None,
    ) -> CrystalDataset:
        """Return a dataset from a multi‑image POSCARS file."""
        from pymatgen.io.ase import AseAtomsAdaptor

        atoms_list = read_poscars(poscars_file)
        if entries_ids_file:
            ids = Path(entries_ids_file).read_text().splitlines()
            if len(ids) != len(atoms_list):
                raise ValueError("Number of IDs does not match number of structures.")
        else:
            ids = [next(self._id_iter) for _ in range(len(atoms_list))]

        entries = [
            CrystalEntry(id=id_, structure=AseAtomsAdaptor.get_structure(ats))
            for id_, ats in zip(ids, atoms_list)
        ]
        msg = message or f"Structures loaded from {poscars_file.as_posix()}"
        return CrystalDataset(entries, message=msg)

    def load_from_extxyz(
        self,
        extxyz_file: Path,
        entries_ids_file: Path | None = None,
        *,
        message: str | None = None,
    ) -> CrystalDataset:
        """Return a dataset from an extended XYZ file."""
        structures = _parse_extxyz_structures(extxyz_file)
        if entries_ids_file:
            ids = Path(entries_ids_file).read_text().splitlines()
            if len(ids) != len(structures):
                raise ValueError("Number of IDs does not match number of structures.")
        else:
            ids = [next(self._id_iter) for _ in range(len(structures))]

        entries = [
            CrystalEntry(id=id_, structure=struct)
            for id_, struct in zip(ids, structures)
        ]
        msg = message or f"Structures loaded from {extxyz_file.as_posix()}"
        return CrystalDataset(entries, message=msg)

    # --------------------------------------------------------------------- #
    # Writer methods                                                        #
    # --------------------------------------------------------------------- #
    @staticmethod
    def dump_poscar_files(
        dataset: Iterable[CrystalEntry],
        out_dir: Path,
    ) -> None:
        """Write each entry to its own POSCAR file in *out_dir*."""
        out_dir.mkdir(parents=True, exist_ok=True)
        for entry in dataset:
            filename = out_dir / entry.poscarname
            entry.structure.to(fmt="poscar", filename=filename, comment=entry.id)

    @staticmethod
    def dump_multiimage_poscar(
        dataset: Iterable[CrystalEntry],
        filename: Path,
        *,
        sort_key=None,
    ) -> None:
        """Write a POSCARS multi‑image file; optional plain‑text ID list."""
        sequence = sorted(dataset, key=sort_key) if sort_key else list(dataset)
        with open(filename, "wt") as fid:
            for entry in sequence:
                poscar_str = Poscar(entry.structure, comment=entry.id).get_str()
                fid.write(poscar_str)
