from __future__ import annotations
import os
from pathlib import Path
from typing import Iterable, Iterator, List, Sequence, Callable, NewType, Tuple, Set
from my_packages.genutils.misc import describe_predicate

from .crystal_entry import CrystalEntry
from datetime import datetime

DatasetID = NewType("DatasetID", str)

class CrystalDataset(Sequence[CrystalEntry]):

    def __init__(
            self,
            entries: Iterable[CrystalEntry],
            dataset_id: DatasetID | None = None,
            parent_ids: Tuple[DatasetID, ...] = (),
            message='',
            base_path=None,
            provenance=None,  # optional handle to provenance store/registry
    ) -> None:
        self._entries = list(entries)
        self.metadata = {"message": message, "created_on": datetime.today().strftime('%Y-%m-%d %H:%M')}
        self.dataset_id = dataset_id or self._generate_id()
        self.parent_ids = parent_ids
        self.base_path = base_path or provenance.base_path if provenance else Path(os.getcwd())
        self._prov = provenance  # can be None
        self._elements = None

    """Thin, read‑oriented container around CrystalEntry objects."""

    # --- Sequence protocol (read‑only operations) -----------------
    def __len__(self) -> int:
        return len(self._entries)

    def __iter__(self) -> Iterator[CrystalEntry]:
        return iter(self._entries)

    def __getitem__(self, idx) -> CrystalEntry:
        return self._entries[idx]

    def _generate_id(self):
        return hex(hash((id(self), datetime.today(), self.metadata["message"])))[2:]

    @property
    def provenance(self):
        return self._prov

    @property
    def elements(self):
        """Set of elements in the dataset."""
        if not self._elements:
            self._elements = set(
                el for e in self if e.structure for el in e.composition.as_data_dict()['elements'])
        return self._elements

    @classmethod
    def from_parents(cls, entries, parents, message, **kwargs):
        parent_ids = tuple(ds.dataset_id for ds in parents)
        if all([ds.base_path == parents[0].base_path for ds in parents]):
            base_path = parents[0].base_path
        else:
            base_path=kwargs.get("base_path", None)
        if all([ds.provenance == parents[0].provenance for ds in parents]):
            provenance = parents[0].provenance
        else:
            provenance=kwargs.get("provenance", None)
        return cls(entries, parent_ids=parent_ids, message=message, provenance=provenance, base_path=base_path)

    # --- Minimal *pure* helpers that return new datasets ----------
    def merge(self, other: "CrystalDataset", message=None) -> "CrystalDataset":
        message = message or 'Simple merge'
        return CrystalDataset.from_parents(self._entries + list(other), parents=(self, other), message=message)

    def filter(self, predicate: Callable[[CrystalEntry], bool], message=None) -> "CrystalDataset":
        message = message or f"Filtered by {describe_predicate(predicate)}"
        return CrystalDataset.from_parents([e for e in self if predicate(e)], parents=(self,), message=message)

