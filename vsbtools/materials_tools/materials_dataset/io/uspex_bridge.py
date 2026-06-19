from pathlib import Path
from functools import lru_cache
import re
import numpy as np
from pymatgen.core.periodic_table import Element
from ..crystal_dataset import CrystalDataset, CrystalEntry
from ...external_paths import (
    add_sys_path,
    import_from_path_validator,
    resolve_external_path,
)
from .structures_dataset_io import StructureDatasetIO


def _normalize_uspex_path(path: Path) -> Path:
    if path.name == "USPEX" and path.is_dir():
        return path.parent
    return path

Atomistic = None
Engine = None
Flavour = None
Entry = None
RadialDistributionUtility = None
TOLERANCE_DEFAULT = None
atomistic = None
_USPEX_IMPORT_ERROR = None


def _import_uspex_modules() -> bool:
    global Atomistic, Engine, Flavour, Entry, RadialDistributionUtility, TOLERANCE_DEFAULT, atomistic
    global _USPEX_IMPORT_ERROR
    if Atomistic is not None:
        return True
    try:
        from USPEX.components import Atomistic as _Atomistic
        from USPEX.DataModel.Engine import Engine as _Engine
        from USPEX.DataModel.Flavour import Flavour as _Flavour
        from USPEX.DataModel.Entry import Entry as _Entry
        from USPEX.Atomistic.RadialDistributionUtility import (
            RadialDistributionUtility as _RadialDistributionUtility,
            TOLERANCE_DEFAULT as _TOLERANCE_DEFAULT,
        )
    except ImportError as exc:
        _USPEX_IMPORT_ERROR = exc
        return False

    Atomistic = _Atomistic
    Engine = _Engine
    Flavour = _Flavour
    Entry = _Entry
    RadialDistributionUtility = _RadialDistributionUtility
    TOLERANCE_DEFAULT = _TOLERANCE_DEFAULT
    Engine.createEngine(":memory:")
    atomistic = Atomistic()
    _USPEX_IMPORT_ERROR = None
    return True


def _ensure_uspex_available(*, prompt: bool = True) -> None:
    if _import_uspex_modules():
        return
    USPEX_PYTHON_PATH = resolve_external_path(
        name="USPEX Python package",
        config_key="uspex_python_path",
        env_var="USPEX_PYTHON_PATH",
        normalizer=_normalize_uspex_path,
        validator=import_from_path_validator("USPEX.components"),
        prompt=prompt,
        required=prompt,
        prompt_text="Enter path to USPEX Python root or USPEX package directory: ",
    )
    if USPEX_PYTHON_PATH is not None:
        add_sys_path(USPEX_PYTHON_PATH)
    if not _import_uspex_modules():
        raise ModuleNotFoundError(
            "USPEX is not configured. Set USPEX_PYTHON_PATH or configure external_paths.json."
        ) from _USPEX_IMPORT_ERROR


def _sanitize_species_symbol(raw_symbol: str) -> str:
    """Convert species labels like Fe2+, O- to plain element symbols."""
    cleaned = re.sub(r"[\d\+\-]", "", raw_symbol).strip()
    cleaned = re.sub(r"[^A-Za-z]", "", cleaned)
    if not cleaned:
        raise ValueError(f"Cannot sanitize empty species symbol from '{raw_symbol}'")
    normalized = cleaned[0].upper() + cleaned[1:].lower()
    if not Element.is_valid_symbol(normalized):
        raise ValueError(f"Sanitized symbol '{normalized}' from '{raw_symbol}' is not a valid element")
    return normalized


class USPEXBridge:
    def __init__(self, elements, legacy=True, tol_FP=None):
        _ensure_uspex_available()
        self.tol_FP = tol_FP or TOLERANCE_DEFAULT
        self.rdu = RadialDistributionUtility(symbols=elements,
                                             suffix='origin', legacy=legacy, tolerance=self.tol_FP, storeDistances=False)
        self.uspex_entry_extensions = dict(atomistic=(atomistic, atomistic.propertyExtension.propertyTable),
                          radialDistributionUtility=(self.rdu, self.rdu.propertyExtension.propertyTable))
        self._sig = (tuple(elements), bool(legacy), self.tol_FP)
        self.id=-1

    def __hash__(self):
        return hash(self._sig)

    def __eq__(self, other):
        return isinstance(other, USPEXBridge) and self._sig == other._sig

    @lru_cache(maxsize=15000)
    def uspex_entry_from_de(self, de_entry: CrystalEntry) -> "Entry":
        types, coords, cell = ([atomistic.atomType(_sanitize_species_symbol(s.species_string)) for s in de_entry.structure],
                               de_entry.structure.cart_coords,
                               atomistic.cellType(de_entry.structure.lattice.matrix, pbc = (1, 1, 1)))
        uspex_structure = atomistic.AtomicStructureRepresentation.structureType(atomTypes=types, coordinates=coords,
                                                                                cell=cell)
        self.id += 1
        return Entry.newEntry(Flavour(extensions=self.uspex_entry_extensions,
                                      **{'.howCome': 'Seeds', '.parent': None, '.label': self.id},
                                      **atomistic.atomicDisassemblerType(
                                          np.arange(len(uspex_structure)).reshape((-1, 1))).disassemble(
                                          uspex_structure)))

    def fp_dist(self, de_entry_1: CrystalEntry, de_entry_2: CrystalEntry) -> float:
        return self.rdu.dist(self.uspex_entry_from_de(de_entry_1), self.uspex_entry_from_de(de_entry_2))

    @staticmethod
    def prepare_seeds(ds: CrystalDataset, seeds_file_path: str | Path | None = None, id_list_file: str | Path | None = None):
        StructureDatasetIO.dump_multiimage_poscar(ds, seeds_file_path)
        if id_list_file:
            with open(id_list_file, 'rt') as id_list_h:
                id_list_h.writelines([e.id for e in ds])

    @staticmethod
    def read_idlist(idlist_path: str | Path | None):
        return idlist_path.read_text().split('\n')
