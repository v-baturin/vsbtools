from __future__ import annotations

import copy
import os
import warnings
from datetime import datetime

"""
unified thermodynamic entry and single‑source dataset
Structures live in memory; POSCARs are written only when
`export_poscars()` is called.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Set, Callable
import logging
from prettytable import PrettyTable
import numpy as np
import pandas as pd
from pymatgen.core import Composition, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Poscar
from USPEX.components import Atomistic
from USPEX.DataModel.Engine import Engine
from USPEX.DataModel.Flavour import Flavour
from USPEX.DataModel.Entry import Entry
from my_packages.materials_tools.uspex_toolkit.remove_duplicates import remove_duplicates, prepare_dist_function


datetime.today().strftime('%Y%m%d')
logger = logging.getLogger(__name__)

# USPEX infrastructure objects
Engine.createEngine(':memory:')
atomistic = Atomistic()
uspex_entry_extensions = dict(atomistic=(atomistic, atomistic.propertyExtension.propertyTable))

# Representation preset
ENTRY_ATTRIBUTE_LABELS = {"origin_id": "ID", "composition": "Composition",
                         "energy_total": "Total energy (eV)",
                          "e_above_hull":  "Energy above hull (eV/atom)",
                          "sym_group_no": "Symmetry group No",
                          "sym_group_symbol": "Symmetry group Symbol",
                          "nonequivalent_sites": "Non-equivalent sites"}

# DATASET_ATTRIBUTE_LABELS = {"e_above_hull_list": "Energy above hull (eV/atom)"}
DATASET_ATTRIBUTE_LABELS = {}
DEFAULT_SYMPREC = 1e-4

class CrystalEntry:

    def __init__(self,
                 origin_id: str,
                 structure: Structure | None = None,
                 composition: Composition | None = None,
                 energy_total: float | None = None,
                 origin: str = 'NA',
                 calc_settings: Mapping[str, Any] | None = None,
                 parent_set: CrystalDataset | None = None,
                 symprec: float = None, **kwargs):
        self.origin_id = origin_id
        self.structure = structure               # in‑memory structure
        self.composition = composition
        self.energy_total = energy_total         # eV per formula unit
        self.origin = origin                       # provenance tag ("oqmd", "mp", …)
        self.calc_settings = calc_settings
        self.parent_set = parent_set
        self._sga = None
        self.symprec = symprec or DEFAULT_SYMPREC


    # ---------------------------------------------------------------------
    # Derived helpers
    # ---------------------------------------------------------------------
    @property
    def natoms(self) -> int:
        return len(self.structure)

    @property
    def formula(self) -> str:
        """Integer formula string, e.g. "Al4Fe2Ni2"."""
        return self.composition.formula.replace(" ","")

    @property
    def e_above_hull(self) -> float:
        return self.parent_set.e_above_hull_pa()[self]

    @property
    def formation_energy(self) -> float:
        return self.parent_set.phase_diagram.get_form_energy(self.as_pd_entry())

    @property
    def formation_energy_pa(self) -> float:
        return self.parent_set.phase_diagram.get_form_energy(self.as_pd_entry()) / self.natoms

    @property
    def sga(self):
        if not hasattr(self, '_sga'):
            self._sga = None
        self._sga = self._sga or SpacegroupAnalyzer(self.structure, symprec=self.symprec)
        return self._sga

    @property
    def sym_group_no(self) -> int:
        return self.sga.get_space_group_number()

    @property
    def sym_group_symbol(self) -> str:
        return self.sga.get_space_group_symbol()

    @property
    def nonequivalent_sites(self) -> dict:
        symm_struct = self.sga.get_symmetrized_structure()
        nonequivalent_positions = {}
        for group in symm_struct.equivalent_sites:
            element = group[0].specie.symbol
            nonequivalent_positions[element] = nonequivalent_positions.get(element, 0) + 1
        return nonequivalent_positions

    def get_symmetrized_structure(self, symprec=1e-4, update_sga=False, primitive=True):
        # debug-------------
        sg_before = self.sym_group_no
        # end debug---------
        sga = SpacegroupAnalyzer(self.structure, symprec=symprec)
        if update_sga:
            self._sga = sga
        # debug-------------
        sg_after = self.sym_group_no
        if sg_before != sg_after:
            print(f"{self.origin_id}: sg_before {sg_before} != sg_after {sg_after}")
        # end debug---------
        if primitive:
            return sga.get_primitive_standard_structure()
        else:
            return sga.get_refined_structure()



    # ------------------------------------------------------------------
    # Conversions
    # ------------------------------------------------------------------
    def as_pd_entry(self) -> PDEntry | None:
        """Return a PhaseDiagram-compatible entry (total energy)."""
        if self.energy_total is not None:
            return PDEntry(self.composition, self.energy_total, attribute=self.origin_id)
        else:
            warnings.warn("No energy provided, thermodynamic properties unavailable")
            return None

    def as_uspex_entry(self):
        atoms = AseAtomsAdaptor.get_atoms(self.structure)
        uspex_structure = atomistic.AtomicStructureRepresentation.fromAtoms(atoms)
        return Entry.newEntry(Flavour(extensions=uspex_entry_extensions,
                              **{'.howCome': 'Seeds', '.parent': None, '.label': self.origin_id},
                              **atomistic.atomicDisassemblerType(np.arange(len(uspex_structure)).reshape((-1, 1))).disassemble(uspex_structure)))
    # ------------------------------------------------------------------
    # Hash / equality
    # ------------------------------------------------------------------
    def __hash__(self):  # noqa: D401
        return hash((self.origin, self.origin_id))
    def __eq__(self, other: object):  # noqa: D401
        if not isinstance(other, CrystalEntry):
            return NotImplemented
        return (
            self.origin == other.origin
            and self.origin_id == other.origin_id
        )

    # ------------------------------------------------------------------
    # Factory from DB row
    # ------------------------------------------------------------------
    @classmethod
    def from_row(cls, row: Mapping[str, Any] | pd.Series | Any, source: str) -> "CrystalEntry":
        """Create a `ThermoEntry` from a DataFrame row / dict / namedtuple."""
        # helper that works for dict, Series, or namedtuple
        def _get(obj: Any, key: str, default: Any = None):
            if isinstance(obj, Mapping):
                return obj.get(key, default)
            return getattr(obj, key, default)

        sid = str(_get(row, "id"))

        # ---------------structure and composition ---------------------
        structure: Structure | None = _get(row, "structure")
        composition: Composition | None = Composition(_get(row, "formula"))
        if structure is None:
            spath = _get(row, "structure_path")
            if spath is None:
                warnings.warn("Row lacks structure object or path")
            else:
                structure = Structure.from_file(Path(spath))

        if structure and composition:
            assert composition == structure.composition, f"{composition} != {structure.composition}"
        elif structure:
            composition = structure.composition
        else:
            raise ValueError("Row lacks both composition and structure")

        # ----------------------- energy -------------------------------
        if (et := _get(row, "e_total")) is not None:
            energy_total = float(et)
        elif (epa := _get(row, "e_pa")) is not None:
            energy_total = float(epa) * len(structure)
        else:
            raise ValueError("Row lacks energy information")

        calc_settings = _get(row, "settings") or _get(row, "calc_settings")

        return cls(
            origin_id=sid,
            energy_total=energy_total,
            composition=composition,
            structure=structure,
            origin=source,
            calc_settings=calc_settings,
        )

    @classmethod
    def from_uspex_entry(cls, uspex_entry: Entry, **kwargs):
        atoms = atomistic.AtomicStructureRepresentation.toAtoms(
            uspex_entry["atomistic.structure.origin"])
        return cls(
            structure=AseAtomsAdaptor.get_structure(atoms),
            **kwargs
        )

    @classmethod
    def from_struc_file(cls, file_path, **kwargs):
        struct = Structure.from_file(file_path)
        return cls(structure=struct, **kwargs)




# ---------------------------------------------------------------------------
# ThermoDataset
# ---------------------------------------------------------------------------


class CrystalDataset(list[CrystalEntry]):
    """Collection of entries from one provenance (total energies)."""

    # ------------------------------------------------------------------
    # Validation / cache helpers
    # ------------------------------------------------------------------
    def _reset_caches(self):
        self._pd_cache: PhaseDiagram | None = None
        self._ehull_cache: Dict[CrystalEntry, float] | None = None
        self._rdf_utility = None
        self._uspex_entries = None

    def _reset_parent(self):
        for e in self:
            e.parent_set = self

    def _validate_single_source(self):
        if not self:
            return
        origins = {e.origin for e in self}
        if len(origins) > 1:
            raise ValueError(f"Dataset must be single-source, got {origins}")

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------
    @classmethod
    def from_client(cls, client: Any, elements: Set[str], *, label: str | None = None) -> "CrystalDataset":
        if hasattr(client, "query"):
            df = client.query(elements)
        else:
            raise AttributeError("Client lacks 'query' method")
        tag = label or client.__class__.__name__.replace("Client", "").lower()
        entries = [CrystalEntry.from_row(r, tag) for r in df.itertuples(index=False)]
        dataset = cls(entries)
        dataset._reset_caches()
        dataset._reset_parent()
        return dataset

    @classmethod
    def from_struct_folder(cls,
                           structures_path: Path| str,
                           search_pattern: str = '*',
                           parameters_dict: Dict[str, list]| None = None):
        structures_path = Path(structures_path)
        found_files = sorted(structures_path.rglob(search_pattern))
        if parameters_dict is not None:
            for val in parameters_dict.values():
                assert len(val) == len(list(found_files))
            params = [dict(zip(parameters_dict.keys(), values)) for values in zip(*parameters_dict.values())]
        else:
            params = [{'origin_id': str(i) } for i in range(len(found_files))]
        entry_list = [CrystalEntry.from_struc_file(f, origin=f.name, **params[i]) for i, f in enumerate(found_files)]
        dataset = cls(entry_list)
        dataset._reset_parent()
        return dataset


    # ------------------------------------------------------------------
    # Properties & heavy operations
    # ------------------------------------------------------------------
    @property
    def source(self) -> Optional[str]:
        return self[0].origin if self else None

    @property
    def uspex_entry_list(self) -> list[Entry]:
        if getattr(self, "_uspex_entries", None) is None:
            self._uspex_entries = [e.as_uspex_entry() for e in self]
        return self._uspex_entries

    @property
    def phase_diagram(self) -> PhaseDiagram | None:
        if not any([e.energy_total for e in self]):
            warnings.warn("No energies provided, thermodynamics unavailable")
            return None
        if getattr(self, "_pd_cache", None) is None:
            self._pd_cache = PhaseDiagram([e.as_pd_entry() for e in self])
        return self._pd_cache  # type: ignore[return-value]

    def e_above_hull_pa(self, refresh: bool = False) -> Dict[CrystalEntry, float | None]:
        if getattr(self, "_ehull_cache", None) is None or refresh:
            pd = self.phase_diagram
            if pd is None:
                return {e: None for e in self}
            self._ehull_cache = {e: pd.get_e_above_hull(e.as_pd_entry()) for e in self}
        return self._ehull_cache  # type: ignore[return-value]

    @property
    def e_above_hull_list(self):
        return [self.e_above_hull_pa()[e] for e in self]

    #------------------------------------------------------------------
    # Filtering
    #------------------------------------------------------------------
    def filter(self, predicate_fn: Callable[[CrystalEntry], bool],
               reset_caches: bool = True,
               reset_parent: bool = True) -> CrystalDataset:
        filtered_set = self.__class__([e for e in self if predicate_fn(e)])
        if reset_caches:
            filtered_set._reset_caches()
        if reset_parent:
            filtered_set._reset_parent()
        return filtered_set

    # ------------------------------------------------------------------
    # Conversion
    # ------------------------------------------------------------------
    def __repr__(self):  # noqa: D401
        return f"<ThermoDataset n={len(self)} source='{self.source}'>"

    def as_pandas_df(self,
                     entry_attribute_labels: dict | None = None,
                     dataset_attribute_labels: dict | None = None):
        # labels = {parameter_label: parameter_name}
        if entry_attribute_labels is None:
            entry_attribute_labels = ENTRY_ATTRIBUTE_LABELS
        if dataset_attribute_labels is None:
            dataset_attribute_labels = DATASET_ATTRIBUTE_LABELS
        data={v:[getattr(e,k) for e in self] for k,v in entry_attribute_labels.items()}
        data.update({v:getattr(self, k) for k,v in dataset_attribute_labels.items()})
        df = pd.DataFrame(data)
        # Identify and drop all-None columns
        all_none_cols = [col for col in df.columns if df[col].isna().all()]
        df.drop(columns=all_none_cols, inplace=True)
        # Report dropped columns
        if all_none_cols:
            print(f"Dropped columns with only None values: {', '.join(all_none_cols)}")
        return df

    # ------------------------------------------------------------------
    # Representation
    # ------------------------------------------------------------------
    def present_as_table(self, dump_path: Path | str | None = None, sort_by=None, columns=None):
        df = self.as_pandas_df().sort_values(by=sort_by)
        if not columns:
            columns = list(df.columns)
        else:
            # Validate that each requested column exists
            missing = set(columns) - set(df.columns)
            if missing:
                raise ValueError(f"Columns {missing} are not in the DataFrame ({df.columns})")
        pt = PrettyTable()
        pt.field_names = columns

        # Slice df to just the chosen columns, then iterate
        for row in df[columns].itertuples(index=False, name=None):
            pt.add_row(row)

        # Write the table string to a file
        if dump_path:
            with open(dump_path, "w") as file:
                file.write(pt.get_string())
        else:
            print(pt)

    # ------------------------------------------------------------------
    # File I/O
    # ------------------------------------------------------------------
    def write_poscar_files(self, directory: Path | str | None = None):
        if directory is None:
            directory = os.getcwd()
        directory = Path(directory).expanduser().resolve()
        directory.mkdir(parents=True, exist_ok=True)
        for e in self:
            fname = f"{e.formula}-{e.origin_id.replace('POSCAR', '')}_POSCAR".split("/")[-1]
            e.structure.to(fmt="poscar", filename=directory / fname, comment=e.origin_id)

    def write_multistructure_poscar(self, filename: Path | str | None = None, sort_key=None):
        if filename is None:
            filename = Path(os.getcwd()) / "POSCARS"
        with open(filename, 'wt') as fid:
            sequence = self if not sort_key else sorted(self, key=sort_key)
            for e in sequence:
                poscar_str = Poscar(e.structure, comment=e.origin_id).get_str()
                fid.write(poscar_str)

    def to_energy_csv(self, path: Path | str | None = None, *, unique: bool = False) -> str:
        df = pd.DataFrame({
            "id": [e.origin_id for e in self],
            "composition": [e.formula for e in self],
            "e_total": [e.energy_total for e in self],
        })
        if unique:
            df = df.loc[df.groupby("composition").e_total.idxmin()].reset_index(drop=True)
        csv_txt = df.to_csv(index=False)
        if path is not None:
            Path(path).write_text(csv_txt)
        return csv_txt

    def __getstate__(self):
        memo = {}
        data_copy = []
        for e in self:
            entry_dict = e.__dict__.copy()
            if "structure" in entry_dict:
                entry_dict["structure_dict"] = entry_dict["structure"].as_dict()
            if "composition" in entry_dict:
                entry_dict["composition_dict"] = entry_dict["composition"].as_dict()
            for to_drop in ["parent_set"]:
                entry_dict.pop(to_drop)
            data_copy.append(entry_dict)
        attrs_copy = self.__dict__.copy()
        return data_copy, attrs_copy

    def __setstate__(self, state):
        data_copy, attrs_copy = state
        self.clear()
        for e in data_copy:
            if "structure_dict" in e:
                e["structure"] = Structure.from_dict(e["structure_dict"])
                e.pop("structure_dict")
            if "composition" in e:
                e["composition"] = Composition.from_dict(e["composition_dict"])
                e.pop("composition_dict")
            e["parent_set"] = self
        self.extend([CrystalEntry(**dct) for dct in data_copy])
        self.__dict__.update(attrs_copy)

    # ------------------------------------------------------------------
    # USPEX-based structural filtration
    # ------------------------------------------------------------------
    def deduplicated(self,
                     check_clusters_file=False,
                     check_dist_matrix_file=True,
                     tolFp=0.08, **kwargs):
        fitness_list = [e.energy_total for e in self]
        if not any(fitness_list):
            fitness_list = None
        best_representatives, clusters, best_idx = remove_duplicates(self.uspex_entry_list, fitness_list,
                                                                     check_clusters_file=check_clusters_file,
                                                                     check_dist_matrix_file=check_dist_matrix_file,
                                                                     tolFp=tolFp, **kwargs)
        filtered_list = [self[i] for i in best_idx]
        return CrystalDataset(filtered_list), clusters, best_idx

    def contains_structure(self, crystal_entry: CrystalEntry, tol_FP = 0.08) -> list:
        uspex_entry_ref = crystal_entry.as_uspex_entry()
        rho = prepare_dist_function(self.uspex_entry_list + [uspex_entry_ref])
        true_idcs = [i for i, uspex_entry in enumerate(self.uspex_entry_list) if rho(uspex_entry_ref, uspex_entry) <= tol_FP ]
        return [self[i].origin for i in true_idcs]





if __name__ == "__main__":

    # from uspexOutput_client import USPEXOutputClient
    # ds = CrystalDataset.from_client(USPEXOutputClient(Path('goodStructures'), mode="goodStructures"),
    #                                 elements={"Mo", "Si", "B"})
    # df = ds.as_pandas_df().sort_values(by="Energy above hull (eV)")
    # pt = PrettyTable()
    # pt.field_names = df.columns.tolist()
    #
    # for row in df.itertuples(index=False):
    #     pt.add_row(row)
    #
    # # Print the table
    # print(pt)


    ############ TERNARY PHASE DIA
    # from db_clients.mp_client import MPClient
    # ds = CrystalDataset.from_client(MPClient(), elements={"Al", "Fe", "Ni"})
    # phase_dia = ds.phase_diagram
    # df = ds.as_pandas_df().sort_values(by="Energy above hull (eV/atom)")
    # pt = PrettyTable()
    # pt.field_names = df.columns.tolist()
    #
    # for row in df.itertuples(index=False):
    #     pt.add_row(row)
    #
    # # Print the table
    # print(pt)
    # plot_ternary_pd(phase_dia, ["Al", "Fe", "Ni"], show_unstable=0.06,
    #             unstable_alpha=0.2)
    ################################

    # ds = CrystalDataset.from_struct_folder("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/STRUCTURE_GENERATION/Mattergen_Auguste/Re__Liste_de_systemes/Vladimir/Mo-Si-B")
    # ref_struct = CrystalEntry.from_struc_file("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/Reference_structures/Mo5SiB2/POSCAR", origin_id="0")
    # print(ds.contains_structure(ref_struct))

    # print(ddc[1])
    # ds.to_energy_csv("MoSiB.csv")
    # print(ds.phase_diagram)
    # for t,v in ds.e_above_hull().items():
    #     print(t.origin_id, v)
    # ds = ThermoDataset.from_client()
    # ds = ThermoDataset([ThermoEntry.from_row(r, "uspex") for r in df.itertuples(index=False)])
    # ds.to_energy_csv("MoSiB.csv")
    # ds.export_poscars("mosib_poscars")
    # print(ds.phase_diagram)
    # for t, v in ds.e_above_hull().items():
    #     print(t.origin_id, v)

    # TESTING write_multistructure_poscar
    import pickle as pkl
    with open("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/Scripts/alexMoSiBP_filtered.pkl", 'rb') as pkldump_fid:
        alex_filtered = pkl.load(pkldump_fid)
    alex_filtered.write_multistructure_poscar(sort_key=lambda e: e.e_above_hull)