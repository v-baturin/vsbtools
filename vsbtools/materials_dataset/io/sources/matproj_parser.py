from dataclasses import dataclass, field
from itertools import combinations
from typing import Any, List, Set, Iterable, Generator
import getpass
import json
import os
from pathlib import Path

import pandas as pd
from pymatgen.core import Structure

_BaseRester = None

_CONFIG_ENV_VAR = "VSBTOOLS_MP_API_KEY_FILE"


def _default_key_path() -> Path:
    """Return the per-user location for the Materials Project API key."""
    configured = os.environ.get(_CONFIG_ENV_VAR)
    if configured:
        return Path(configured).expanduser()
    config_home = os.environ.get("XDG_CONFIG_HOME") or os.environ.get("APPDATA")
    root = Path(config_home).expanduser() if config_home else Path("~/.config").expanduser()
    return root / "vsbtools" / "materials_project.json"


def _load_saved_key(path: Path) -> str | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8")).get("api_key")
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None
    return value.strip() if isinstance(value, str) and value.strip() else None


def _save_key(path: Path, key: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"api_key": key}, indent=2) + "\n", encoding="utf-8")
    if os.name != "nt":
        path.chmod(0o600)


def _is_authentication_error(error: Exception) -> bool:
    text = str(error).casefold()
    return any(marker in text for marker in ("401", "403", "unauthorized", "forbidden", "api key", "apikey"))


@dataclass(slots=True)
class MPClient:
    """Materials Project v2 client (total-energy rows for ThermoDataset)."""

    api_key: str | None = None
    key_path: Path | None = None
    prompt_for_key: bool = True
    _mpr: Any = field(init=False, default=None, repr=False)

    # ------------------------------------------------------------------ #
    # Internals                                                          #
    # ------------------------------------------------------------------ #

    def _key_file(self) -> Path:
        return self.key_path or _default_key_path()

    def _prompt_for_key(self, *, reason: str | None = None) -> str:
        if not self.prompt_for_key:
            detail = f" ({reason})" if reason else ""
            raise RuntimeError(
                "Materials Project API key is required" + detail
                + ". Set MAPI_KEY or PMG_MAPI_KEY, pass api_key, or save one in "
                + str(self._key_file())
            )
        if reason:
            print(f"Materials Project authentication failed: {reason}", file=sys.stderr)
        try:
            key = getpass.getpass("Materials Project API key (input hidden): ").strip()
        except EOFError as exc:
            raise RuntimeError(
                "Materials Project API key is required but this process cannot accept input. "
                "Set MAPI_KEY or PMG_MAPI_KEY before running non-interactively."
            ) from exc
        if not key:
            raise RuntimeError("Materials Project API key was not provided.")
        _save_key(self._key_file(), key)
        return key

    def _resolve_api_key(self) -> str:
        key = self.api_key or os.getenv("MAPI_KEY") or os.getenv("PMG_MAPI_KEY") or _load_saved_key(self._key_file())
        return key if key else self._prompt_for_key()

    def _connect(self):
        global _BaseRester
        if _BaseRester is None:
            try:
                from mp_api.client import MPRester as _BaseRester
            except ImportError as err:
                raise ImportError("Install `mp-api` to use MPClient.") from err
        if self._mpr is None:
            self._mpr = _BaseRester(api_key=self._resolve_api_key())
        return self._mpr

    def _retry_after_authentication_failure(self, error: Exception) -> bool:
        if self.api_key is not None or not _is_authentication_error(error):
            return False
        self._mpr = None
        self.api_key = self._prompt_for_key(reason=str(error))
        return True

    @staticmethod
    def _subspaces(elements: Set[str]) -> Generator[str, Any, None]:
        for r in range(1, len(elements) + 1):
            for combo in combinations(sorted(elements), r):
                yield "-".join(combo)

    # ------------------------------------------------------------------ #
    # Public API                                                         #
    # ------------------------------------------------------------------ #

    def query(self, elements: Iterable[str]) -> pd.DataFrame:
        """
        Return every calculation that contains *elements*.

        Columns: id, formula (full), e_total (eV), natoms, structure, source
        """
        try:
            return self._query_once(elements)
        except Exception as error:
            if self._retry_after_authentication_failure(error):
                return self._query_once(elements)
            raise

    def _query_once(self, elements: Iterable[str]) -> pd.DataFrame:
        rester = self._connect()

        thermo_docs = []
        for space in self._subspaces(set(elements)):
            thermo_docs.extend(
                rester.thermo.search(
                    chemsys=space,
                    thermo_types=['GGA_GGA+U'],  # TODO: check other xc-functionals
                    fields=["material_id", "formula_pretty", "energy_per_atom"],
                )
            )

        if not thermo_docs:
            return pd.DataFrame(
                columns=["id", "formula", "energy", "structure", "metadata"]
            )

        # One bulk structures call
        mids = {d.material_id for d in thermo_docs}
        struct_map = {
            s.material_id: s.structure
            for s in rester.summary.search(
                material_ids=list(mids), fields=["material_id", "structure"]
            )
        }

        rows = []
        for doc in thermo_docs:
            struct: Structure | None = struct_map.get(doc.material_id)
            if struct is None:
                continue
            natoms = len(struct)
            e_total = doc.energy_per_atom * natoms
            rows.append(
                {
                    "id": str(doc.material_id),
                    "formula": struct.composition.formula,
                    "energy": e_total,
                    "structure": struct,
                    "metadata": {"source": "MaterialsProject"},
                }
            )

        return pd.DataFrame(rows)

    # Legacy alias for old code paths
    def lowest_energy(self, elements: Set[str]) -> pd.DataFrame:
        return self.query(elements)
