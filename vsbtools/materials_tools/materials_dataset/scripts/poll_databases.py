import os
from typing import Dict
from collections import defaultdict
from pathlib import Path
from ....genutils.misc import serialize_structure, is_substructure
from ..crystal_dataset import CrystalDataset
from ..io.yaml_csv_poscars import read, write
from ..io.preset_loaders import load_from_alexandria, load_from_oqmd, load_from_materials_project
from ..analysis.phase_diagram_tools import PhaseDiagramTools

HOME = Path.home()
CACHE_DIR = HOME / ".cache" / "vsbtools" / "DB_caches"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
MAX_EHULL = 0.02  # Maximum energy above hull in eV/atom for filtering
TOL_FP = 0.008
LOADERS = {"al": load_from_alexandria, "oq": load_from_oqmd, "mp": load_from_materials_project,
           "ma": load_from_materials_project}  # List of database clients to fetch data from
# ------------------------------------------------------------------ #

def poll_databases(elements,
                   database_names=None,
                   pref_db: str = 'al',
                   do_ehull_filtering=True,
                   max_ehull=None,
                   do_deduplication=True,
                   similarity_tk=None,
                   tol_FP: float | None = None,
                   loader_kwargs: Dict | None =None,
                   cache_root_path: Path | None = None,
                   **kwargs):
    """
    Fetch data from Alexandria, OQMD, and Materials Project databases.
    Data is first taken from the reference database (Alexandria) then only the structures unseen in the reference DB
    are added.
    """

    parameters_dict: Dict[str, str|float|int] = {'e_hull_filtering': do_ehull_filtering, 'deduplication': do_deduplication}
    if do_deduplication:
        parameters_dict['tol_FP'] = similarity_tk.tol_FP if tol_FP is None else tol_FP

    if do_ehull_filtering:
        max_ehull = MAX_EHULL if max_ehull is None else max_ehull
        parameters_dict['max_ehull'] = max_ehull

    cache_root_path = cache_root_path if cache_root_path is not None else CACHE_DIR
    cache_name = f"{'-'.join(sorted(elements))}__{serialize_structure(parameters_dict)}"
    cache_base_path = cache_root_path / cache_name

    if cache_base_path is not None and (cache_base_path / "manifest.yaml").exists():
        polled_db = read(cache_base_path / "manifest.yaml")
        try:
            if set(polled_db.elements) == set(elements) and \
                is_substructure(polled_db.metadata, parameters_dict):
                print(f"Data for elements {' '.join(elements)} read from cache")
                return polled_db
        except (AssertionError, KeyError):
            print("Incorrect db file")
        print("Failed to read cached DB")


    if loader_kwargs is None: loader_kwargs = dict()

    pref_db = pref_db[:2].casefold()
    database_names = database_names or ['alexandria', 'oqmd', 'MatProj']  # Default databases to fetch data from
    short_names_dict = {name[:2].casefold(): name for name in database_names}
    loader_kw = defaultdict(dict, {k.casefold()[:2]: v for k, v in loader_kwargs.items()})
    for db in short_names_dict:
        assert isinstance(db, str), "Database names must be strings."
        assert db in LOADERS, f"Database '{short_names_dict[db]}' is not recognized. Available databases: {list(LOADERS.keys())}."
    # other_loaders = [LOADERS[short_name] for short_name in short_names_dict if not short_name.startswith(ref_db)]

    reference_loader = LOADERS[pref_db]

    print(f"Fetching data from preferred DB: {short_names_dict.pop(pref_db)}...")
    reference_data = reference_loader(elements, **loader_kw[pref_db])

    if do_ehull_filtering:
        max_ehull = max_ehull or MAX_EHULL
        pd_tools = PhaseDiagramTools(reference_data)
        reference_data = reference_data.filter(lambda e: pd_tools.height_above_hull_pa(e) < max_ehull)

    for loader in short_names_dict: # Skip the first loader as it's already processed
        print(f"Fetching data from {short_names_dict[loader]}...")
        loaded_ds = LOADERS[loader](elements, **loader_kw[loader])
        pd_tools = PhaseDiagramTools(loaded_ds)
        if do_ehull_filtering:
            loaded_ds = loaded_ds.filter(lambda e: pd_tools.height_above_hull_pa(e) < max_ehull)
        if do_deduplication:
            unseen = similarity_tk.get_unseen_in_ref(loaded_ds, reference_data, tol_FP=tol_FP)
            loaded_ds = unseen
        reference_data = reference_data.merge(loaded_ds)

    if do_deduplication:
        os.makedirs(cache_base_path, exist_ok=True)
        reference_data.override_base_path(cache_base_path)
        reference_data, _, _ = similarity_tk.deduplicate(reference_data, tol_FP=tol_FP)

    msg = (f"Gathered from {', '.join(database_names)} databases "
           f"with elements: {', '.join(elements)}")
    ds = CrystalDataset([e.copy_with(**{'energy': None}) for e in reference_data], message=msg)
    ds.metadata = {**ds.metadata, **parameters_dict}
    write(ds, enforce_base_path=cache_base_path)
    return ds
