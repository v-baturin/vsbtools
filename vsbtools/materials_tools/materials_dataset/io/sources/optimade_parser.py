from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
import json
from typing import Any, Callable, Dict, Iterable, Iterator, Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urljoin, urlparse, urlunparse
from urllib.request import Request, urlopen

import pandas as pd
from pymatgen.core import Structure

from ...analysis.similarity_tools import SimilarityTools
from ...analysis.structural_distance.dscribe_bridge import DScribeBridge
from ...crystal_dataset import CrystalDataset
from ...crystal_entry import CrystalEntry


@dataclass(slots=True)
class OptimadeProvider:
    name: str
    base_url: str
    source: str
    energy_fields: Sequence[str]


@dataclass(slots=True)
class OptimadeClient:
    """Query remote OPTIMADE endpoints and return rows ready for CrystalDataset."""

    providers: Sequence[str] | str = ("alexandria", "materials_project", "oqmd")
    alexandria_functional: str = "pbe"
    energy_mode: str = "formation"
    mp_thermo_type: str = "gga_gga+u"
    page_limit: int = 500
    timeout: float = 60.0
    skip_missing_energy: bool = True
    strict_structures: bool = True
    do_deduplication: bool = True
    show_progress: bool = True
    progress_provider_no: int | None = None
    progress_total_providers: int | None = None
    tol_FP: float | None = None
    similarity_tk: SimilarityTools | None = None
    similarity_bridge_kwargs: Dict[str, Any] = field(default_factory=dict)

    _provider_defs: Dict[str, OptimadeProvider] = field(init=False, repr=False)

    def __post_init__(self):
        if self.energy_mode not in {"formation", "total"}:
            raise ValueError("energy_mode must be 'formation' or 'total'")

        alex_base = {
            "pbe": "https://alexandria.icams.rub.de/pbe",
            "pbesol": "https://alexandria.icams.rub.de/pbesol",
        }.get(self.alexandria_functional.lower())
        if alex_base is None:
            raise ValueError("alexandria_functional must be 'pbe' or 'pbesol'")

        self._provider_defs = {
            "materials_project": OptimadeProvider(
                name="materials_project",
                base_url="https://optimade.materialsproject.org",
                source="MaterialsProject/OPTIMADE",
                energy_fields=("_mp_stability",),
            ),
            "mp": OptimadeProvider(
                name="materials_project",
                base_url="https://optimade.materialsproject.org",
                source="MaterialsProject/OPTIMADE",
                energy_fields=("_mp_stability",),
            ),
            "oqmd": OptimadeProvider(
                name="oqmd",
                base_url="https://oqmd.org/optimade/",
                source="OQMD/OPTIMADE",
                energy_fields=("_oqmd_delta_e",),
            ),
            "alexandria": OptimadeProvider(
                name="alexandria",
                base_url=alex_base,
                source=f"Alexandria/{self.alexandria_functional.upper()}/OPTIMADE",
                energy_fields=("_alexandria_formation_energy_per_atom", "_alexandria_energy"),
            ),
        }

    def query(self, elements: Iterable[str]) -> pd.DataFrame:
        elements = set(elements)
        rows = []
        similarity_tk = self._similarity_tools(elements) if self.do_deduplication else None
        providers = self._providers()
        subspaces = list(self._subspaces(elements))
        total_queries = len(providers) * len(subspaces)
        status_len = 0

        def progress(status: str):
            nonlocal status_len
            if not self.show_progress:
                return
            print(f"\r{status}{' ' * max(0, status_len - len(status))}", end="", flush=True)
            status_len = len(status)

        for provider_no, provider in enumerate(providers, start=1):
            provider_rows = []
            display_provider_no = self.progress_provider_no or provider_no
            display_total_providers = self.progress_total_providers or len(providers)
            provider_label = f"provider {display_provider_no}/{display_total_providers} {provider.name}"
            progress(f"OPTIMADE {provider_label}: starting")
            for subset_no, subset in enumerate(subspaces, start=1):
                query_no = (provider_no - 1) * len(subspaces) + subset_no
                chemsys = "-".join(subset)
                progress(
                    f"OPTIMADE {provider_label}: querying {chemsys} "
                    f"subspace {subset_no}/{len(subspaces)}, accepted {len(provider_rows)}"
                )
                subset_rows = self._query_provider(
                    provider,
                    subset,
                    progress=progress,
                    query_no=query_no,
                    total_queries=total_queries,
                    provider_label=provider_label,
                    provider_accepted_offset=len(provider_rows),
                )
                provider_rows.extend(subset_rows)
                progress(
                    f"OPTIMADE {provider_label}: finished {chemsys}, "
                    f"accepted {len(provider_rows)}"
                )

            if similarity_tk is not None and rows:
                progress(
                    f"OPTIMADE {provider_label}: comparing against previous providers, "
                    f"accepted {len(provider_rows)}"
                )
                provider_rows = self._unseen_rows(
                    provider_rows,
                    rows,
                    similarity_tk,
                    progress=progress,
                    progress_prefix=f"OPTIMADE {provider_label}: comparing against previous providers",
                )
            rows.extend(provider_rows)
            progress(
                f"OPTIMADE {provider_label}: finished, kept {len(provider_rows)}, "
                f"total accepted {len(rows)}"
            )

        if self.show_progress and status_len:
            print()

        if not rows:
            return pd.DataFrame(columns=["id", "formula", "energy", "structure", "metadata"])

        return pd.DataFrame(rows).reset_index(drop=True)

    # legacy alias
    def lowest_energy(self, elements: set[str]) -> pd.DataFrame:
        return self.query(elements)

    def _similarity_tools(self, elements: Iterable[str]) -> SimilarityTools:
        if self.similarity_tk is not None:
            return self.similarity_tk

        bridge_kwargs = dict(self.similarity_bridge_kwargs)
        if self.tol_FP is not None:
            bridge_kwargs.setdefault("tol_FP", self.tol_FP)
        bridge = DScribeBridge(elements, **bridge_kwargs)
        return SimilarityTools(bridge.fp_dist, bridge.tol_FP)

    def _unseen_rows(
            self,
            rows: list[dict],
            reference_rows: list[dict],
            similarity_tk: SimilarityTools,
            progress: Callable[[str], None] | None = None,
            progress_prefix: str | None = None,
    ) -> list[dict]:
        def mark_duplicate(row: dict, duplicate_of: dict):
            metadata = dict(row.get("metadata") or {})
            metadata["duplicate_of"] = duplicate_of.get("id")
            row["metadata"] = metadata

        return similarity_tk.get_unseen_successively(
            rows,
            reference_rows,
            entry_factory=OptimadeClient._entry_from_row,
            progress=progress,
            progress_prefix=progress_prefix,
            on_duplicate=mark_duplicate,
            skip_errors=True,
        )

    @staticmethod
    def _entry_from_row(row: dict) -> CrystalEntry:
        return CrystalEntry(
            id=row["id"],
            formula=row.get("formula"),
            energy=row.get("energy"),
            structure=row.get("structure"),
            metadata=row.get("metadata") or {},
        )

    @classmethod
    def _dataset_from_rows(cls, rows: list[dict], message: str = "") -> CrystalDataset:
        return CrystalDataset([cls._entry_from_row(row) for row in rows], message=message)

    def _providers(self) -> list[OptimadeProvider]:
        names = [self.providers] if isinstance(self.providers, str) else list(self.providers)
        if any(str(name).lower() == "all" for name in names):
            names = ["materials_project", "oqmd", "alexandria"]

        providers = []
        for name in names:
            key = str(name).lower()
            try:
                providers.append(self._provider_defs[key])
            except KeyError as err:
                known = ", ".join(sorted(self._provider_defs))
                raise ValueError(f"Unknown OPTIMADE provider '{name}'. Known providers: {known}") from err
        return providers

    @staticmethod
    def _subspaces(elements: set[str]) -> Iterator[tuple[str, ...]]:
        for r in range(1, len(elements) + 1):
            yield from combinations(sorted(elements), r)

    def _query_provider(
            self,
            provider: OptimadeProvider,
            elements: tuple[str, ...],
            progress: Callable[[str], None] | None = None,
            query_no: int | None = None,
            total_queries: int | None = None,
            provider_label: str | None = None,
            provider_accepted_offset: int = 0,
    ) -> list[dict]:
        fields = [
            "id",
            "chemical_formula_reduced",
            "lattice_vectors",
            "cartesian_site_positions",
            "species_at_sites",
            "species",
            "structure_features",
            "nsites",
            *provider.energy_fields,
        ]
        if provider.name == "oqmd":
            fields.append("_oqmd_entry_id")
        params = {
            "filter": self._elements_filter(elements),
            "page_limit": self.page_limit,
            "response_fields": ",".join(dict.fromkeys(fields)),
        }
        url = self._endpoint(provider.base_url, "v1/structures") + "?" + urlencode(params)

        rows = []
        page_no = 0
        query_label = f"{query_no}/{total_queries}" if query_no is not None and total_queries is not None else "?"
        chemsys = "-".join(elements)
        provider_label = provider_label or provider.name
        while url:
            page_no += 1
            if progress is not None:
                progress(
                    f"OPTIMADE {provider_label}: {chemsys} query {query_label} "
                    f"page {page_no}: requesting, accepted {provider_accepted_offset + len(rows)}"
                )
            payload = self._get_json(url)
            page_items = len(payload.get("data", []))
            accepted_before = len(rows)
            for item in payload.get("data", []):
                row = self._row_from_item(provider, item)
                if row is not None:
                    rows.append(row)
            if progress is not None:
                progress(
                    f"OPTIMADE {provider_label}: {chemsys} query {query_label} "
                    f"page {page_no}: fetched {page_items}, "
                    f"accepted +{len(rows) - accepted_before}, "
                    f"provider accepted {provider_accepted_offset + len(rows)}"
                )
            url = self._next_url(url, payload)
        return rows

    @staticmethod
    def _elements_filter(elements: tuple[str, ...]) -> str:
        quoted = ", ".join(f'"{el}"' for el in elements)
        return f"elements HAS ALL {quoted} AND nelements={len(elements)}"

    @staticmethod
    def _endpoint(base_url: str, path: str) -> str:
        return urljoin(base_url.rstrip("/") + "/", path)

    def _get_json(self, url: str) -> dict:
        req = Request(url, headers={"Accept": "application/json", "User-Agent": "vsbtools-optimade"})
        try:
            with urlopen(req, timeout=self.timeout) as fh:
                return json.loads(fh.read().decode("utf-8"))
        except HTTPError as err:
            raise RuntimeError(f"OPTIMADE request failed with HTTP {err.code}: {url}") from err
        except URLError as err:
            raise RuntimeError(f"OPTIMADE request failed: {url}") from err

    @staticmethod
    def _next_url(current_url: str, payload: dict) -> str | None:
        next_url = (payload.get("links") or {}).get("next")
        if not next_url:
            return None

        current = urlparse(current_url)
        next_url = urljoin(current_url, next_url)
        parsed = urlparse(next_url)
        if current.scheme == "https" and parsed.scheme == "http":
            parsed = parsed._replace(scheme="https")
        return urlunparse(parsed)

    def _row_from_item(self, provider: OptimadeProvider, item: dict) -> dict | None:
        attrs = item.get("attributes") or {}
        structure = self._structure_from_attributes(attrs)
        if structure is None:
            return None

        energy, energy_meta = self._energy(provider, attrs, len(structure))
        if energy is None:
            if self.skip_missing_energy:
                return None
            raise ValueError(f"No usable {self.energy_mode} energy for {provider.name}:{item.get('id')}")

        optimade_id = item.get("id")
        source_entry_id = self._source_entry_id(provider, attrs, optimade_id)
        metadata = {
            "source": provider.source,
            "optimade_provider": provider.name,
            "optimade_id": optimade_id,
            "energy_mode": self.energy_mode,
            **energy_meta,
        }
        if provider.name == "oqmd" and attrs.get("_oqmd_entry_id") is not None:
            metadata["oqmd_entry_id"] = attrs.get("_oqmd_entry_id")
        return {
            "id": source_entry_id,
            "formula": attrs.get("chemical_formula_reduced") or structure.composition.formula,
            "energy": energy,
            "structure": structure,
            "metadata": metadata,
        }

    @staticmethod
    def _source_entry_id(provider: OptimadeProvider, attrs: dict, optimade_id: Any) -> str:
        if provider.name == "oqmd" and attrs.get("_oqmd_entry_id") is not None:
            return f"oqmd_{attrs['_oqmd_entry_id']}"
        return f"optimade_{provider.name}_{optimade_id}"

    def _structure_from_attributes(self, attrs: dict) -> Structure | None:
        features = set(attrs.get("structure_features") or [])
        if self.strict_structures and features:
            return None

        lattice = attrs.get("lattice_vectors")
        positions = attrs.get("cartesian_site_positions")
        species_at_sites = attrs.get("species_at_sites")
        if not lattice or not positions or not species_at_sites:
            return None

        species_map = self._species_map(attrs.get("species") or [])
        species = [species_map.get(name, name) for name in species_at_sites]
        if any(symbol in {"X", "vacancy"} for symbol in species):
            return None

        try:
            return Structure(lattice, species, positions, coords_are_cartesian=True)
        except Exception:
            if self.strict_structures:
                return None
            raise

    @staticmethod
    def _species_map(species: Sequence[dict]) -> dict[str, str]:
        mapped = {}
        for spec in species:
            symbols = spec.get("chemical_symbols") or []
            concentrations = spec.get("concentration") or []
            if len(symbols) == 1 and (not concentrations or concentrations[0] == 1):
                mapped[spec.get("name")] = symbols[0]
        return mapped

    def _energy(self, provider: OptimadeProvider, attrs: dict, natoms: int) -> tuple[float | None, dict]:
        if provider.name == "materials_project":
            return self._mp_energy(attrs, natoms)
        if provider.name == "oqmd":
            return self._per_atom_energy(attrs, "_oqmd_delta_e", natoms)
        if provider.name == "alexandria":
            if self.energy_mode == "total":
                return self._raw_energy(attrs, "_alexandria_energy")
            return self._per_atom_energy(attrs, "_alexandria_formation_energy_per_atom", natoms)
        return None, {}

    def _mp_energy(self, attrs: dict, natoms: int) -> tuple[float | None, dict]:
        if self.energy_mode == "total":
            return None, {}

        stability = attrs.get("_mp_stability") or {}
        thermo = stability.get(self.mp_thermo_type)
        if thermo is None and stability:
            thermo_key, thermo = next(iter(stability.items()))
        else:
            thermo_key = self.mp_thermo_type
        if not thermo:
            return None, {}

        e_pa = thermo.get("formation_energy_per_atom")
        if e_pa is None:
            return None, {}
        return float(e_pa) * natoms, {
            "energy_field": "_mp_stability.formation_energy_per_atom",
            "mp_thermo_type": thermo_key,
        }

    @staticmethod
    def _per_atom_energy(attrs: dict, field: str, natoms: int) -> tuple[float | None, dict]:
        e_pa = attrs.get(field)
        if e_pa is None:
            return None, {}
        return float(e_pa) * natoms, {"energy_field": field}

    @staticmethod
    def _raw_energy(attrs: dict, field: str) -> tuple[float | None, dict]:
        energy = attrs.get(field)
        if energy is None:
            return None, {}
        return float(energy), {"energy_field": field}
