from pathlib import Path
import json
import hashlib
from typing import Callable, Any, Dict
from ..io.yaml_csv_poscars import read, write
from ..converters import df2ds
from .sources.matproj_parser import MPClient
from .sources.alexandria_parser import AlexandriaClient
from .sources.uspex_output_parser import USPEXOutputClient
from .sources.oqmd_parser import OQMDClient
from .sources.structure_and_energy_files_reader import CSV_and_POSCARS_client

mp_client = MPClient()
oqmd_client = OQMDClient()
HOME = Path.home()
CACHE_DIR = HOME / ".cache" / "vsbtools" / "DB_caches"
CACHE_DIR.mkdir(parents=True, exist_ok=True)  # ensure tree exists


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def _normalize_elements(elements) -> str:
    """
    Canonical, hyphen-joined (sorted) element key.
    E.g., 'Si-O', ['O','Si'] -> 'O-Si'
    """
    if isinstance(elements, str):
        s = elements.strip().replace(",", "-").replace("_", "-")
        parts = [t for t in s.split("-") if t]
    else:
        try:
            parts = [str(x).strip() for x in elements]
        except TypeError:
            parts = [str(elements)]
    parts = sorted({p for p in parts if p})
    return "-".join(parts) if parts else "UNSPEC"

def _kwargs_hash(kwargs: Dict[str, Any]) -> str:
    """
    Stable short hash of kwargs, excluding non-semantic controls.
    """
    filtered = {k: v for k, v in kwargs.items() if k not in ("force_refresh", "message")}
    if not filtered:
        return "none"
    payload = json.dumps(filtered, sort_keys=True, default=repr).encode("utf-8")
    return hashlib.sha1(payload).hexdigest()[:10]


# --------------------------------------------------------------------------------------
# Decorator: uses your read(dir/manifest.yaml) and write(dataset, enforce_base_path=dir, ...)
# --------------------------------------------------------------------------------------
def cache_loader(read_fn, write_fn, *, cache_root: Path = CACHE_DIR, manifest_name: str = "manifest.yaml"):
    """
    Decorator for loaders returning a CrystalDataset.
    - read_fn(manifest_path: str) -> CrystalDataset
    - write_fn(dataset: CrystalDataset, enforce_base_path: str|Path|None = None, comment=None, **kwargs) -> None

    Special kwarg on wrapped loader:
      force_refresh=True -> bypass cache
    """
    def decorator(fn):
        source = fn.__name__.split("load_from_", 1)[-1] if "load_from_" in fn.__name__ else fn.__name__

        def wrapper(elements, message=None, **kwargs):
            force_refresh = bool(kwargs.get("force_refresh", False))
            elements_key = _normalize_elements(elements)
            kw_hash = _kwargs_hash(kwargs)
            folder_name = f"{source}_{elements_key}_{kw_hash}"
            folder = Path(cache_root) / folder_name
            manifest_path = folder / manifest_name

            # 1) Try cache
            if not force_refresh and manifest_path.exists():
                return read_fn(str(manifest_path))

            # 2) Miss: compute dataset via underlying loader
            #    (Do not pass control-only kwarg to the client)
            client_kwargs = {k: v for k, v in kwargs.items() if k not in ("force_refresh",)}
            dataset = fn(elements, message=message, **client_kwargs)

            # 3) Persist using your writer into the cache folder
            folder.mkdir(parents=True, exist_ok=True)
            write_fn(dataset, enforce_base_path=folder, comment=message)

            return dataset

        wrapper.__name__  = fn.__name__
        wrapper.__doc__   = fn.__doc__
        wrapper.__module__= fn.__module__
        return wrapper
    return decorator


@cache_loader(read, write)
def load_from_materials_project(elements, message=None, **kwargs):
    df = mp_client.query(elements)
    message = message or f"Full {elements} system from Materials Project"
    return df2ds(df, message=message)

@cache_loader(read, write)
def load_from_alexandria(elements, message=None, **kwargs):
    alex_client = AlexandriaClient(**kwargs)
    df = alex_client.query(elements)
    message = message or f"Full {elements} system from Alexandria database"
    return df2ds(df, message=message)

@cache_loader(read, write)
def load_from_oqmd(elements, message=None, **kwargs):
    df = oqmd_client.query(elements)
    message = message or f"Full {elements} system from OQMD"
    return df2ds(df, message=message)

def load_from_uspex_calc_folders(calcfolds_path: Path, stage=None, message=None, elements=None):
    uspex_client = USPEXOutputClient(calcfolds_path, stage=stage, mode='calcFolds')
    df = uspex_client.query(elements=elements)
    message = message or f"CalcFolds" +  (f' with elements={elements}' if elements else '.')
    return df2ds(df, message=message)

def load_from_uspex_goodstructures(goodstructures_path: Path, message=None, elements=None):
    uspex_client = USPEXOutputClient(goodstructures_path, mode='goodStructures')
    df = uspex_client.query(elements=elements)
    message = message or f"goodStructures" +  (f' with elements={elements}' if elements else '.')
    return df2ds(df, message=message)

def load_mattersim_estimated_set(csv_path: Path, poscars_base_path: Path, message=None):
    energy_struct_client = CSV_and_POSCARS_client(results_csv=csv_path, poscars_parent_path=poscars_base_path)
    df = energy_struct_client.query()
    message = message or f"Mattersim-estimated structures from {poscars_base_path.resolve()}"
    return df2ds(df, message=message)