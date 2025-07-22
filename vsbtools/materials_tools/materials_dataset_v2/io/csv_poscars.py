import os
import warnings
from dataclasses import dataclass
from pathlib import Path
import yaml
import pandas as pd
from numba.core.types import double
from pymatgen.core import Structure
from ..crystal_dataset import CrystalDataset, CrystalEntry


def _get_str_value_for_csv(e: CrystalEntry, attr_name):
    if attr_name == 'structure':
        val = f"{e.id}POSCAR"
    else:
        val = getattr(e, attr_name)
    if isinstance(val, float):
        val = f"{val:.4f}"
    else:
        val = str(val)
    return val

def _remove_private_keys(dct, recursive=True):
    for key in [k for k in dct if k.startswith('_')]:
        if isinstance(dct[key], dict) and recursive:
            _remove_private_keys(dct[key])
        else:
            dct.pop(key)

def _prepare_for_yaml_write(dct: dict):
    for k, v in dct.items():
        if isinstance(v, dict):
            _prepare_for_yaml_write(v)
        elif not bool(v):
            dct[k] = None
        elif isinstance(v, Path):
            dct[k] = v.as_posix()
        elif isinstance(v, float):
            dct[k] = f"{v:.4f}"



def read(dataset_description: str | Path) -> CrystalDataset:
    def _rebase_path(pth: Path, base: Path):
        return pth if pth.is_absolute() else base / pth
    with open(dataset_description, 'rt') as ds_h:
        dataset_info = yaml.safe_load(ds_h.read())
    dataset_info.pop("comment")
    _remove_private_keys(dataset_info)
    entry_description_file = _rebase_path(Path(dataset_info.pop('entries_csv')), dataset_description.parent)
    poscars_dir = _rebase_path(Path(dataset_info.pop('poscars')), dataset_description.parent)
    dataset = read_csv_poscars(entry_description_file, poscars_dir=poscars_dir, use_fname_as_id=None)
    for k in dataset_info:
        if k == "parent_ids":
            dataset_info["parent_ids"] = tuple(dataset_info["parent_ids"]) if dataset_info["parent_ids"] is not None \
                else tuple()
        elif k == 'base_path':
            dataset_info["base_path"] = Path(dataset_info["base_path"])
        setattr(dataset, k, dataset_info[k])
    return dataset

def write(dataset: CrystalDataset, dataset_description: str | Path | None, comment=None):
    dataset_file = dataset_description or Path(os.getcwd()) / dataset.dataset_id + '.yaml'
    csv_file = dataset_file.with_name(dataset_file.stem + '.csv')
    poscars_path = dataset_file.with_name(dataset_file.stem + 'POSCARS')
    write_csv_poscars(dataset, csv_file, poscars_path)
    ds_dict = dataset.__dict__.copy()
    _remove_private_keys(ds_dict)
    _prepare_for_yaml_write(ds_dict)
    ds_dict["entries_csv"] = csv_file.relative_to(dataset_file.parent).as_posix()
    ds_dict["poscars"] = poscars_path.relative_to(dataset_file.parent).as_posix()
    with open(dataset_description, 'wt') as ds_desc:
        yaml.dump(ds_dict, ds_desc)








def read_csv_poscars(csv_file, poscars_dir=None, use_fname_as_id=True, provenance=None) -> CrystalDataset:
    csv_file = Path(csv_file)
    poscars_dir = poscars_dir or csv_file.with_name(csv_file.stem + 'POSCARS')
    dataset_ID = csv_file.stem if fname_as_id else None
    with open(csv_file, 'rt') as csv_h:
        labels = [h.strip() for h in csv_h.read().split(',')]
    entries = []
    for i, line in enumerate(csv_h):
        entry_data = {k: v for (k,v) in zip(labels, [h.strip() for h in line.split(',')])}
        if not "id" in labels:
            entry_data["id"] = str(i)
        if 'energy' in entry_data:
            entry_data['energy'] = float(entry_data['energy'])
        if 'structure' in entry_data:
            entry_data['structure'] = Structure.from_file(poscars_dir / entry_data['structure'])
            del entry_data['structure']
        entries.append(CrystalEntry(**entry_data))
    return CrystalDataset(entries,
                          dataset_id=dataset_ID,
                          message=f'Loaded from {csv_file}',
                          base_path=csv_file.parent,
                          provenance=provenance
                          )

def write_csv_poscars(ds: CrystalDataset,
          csv_file: str | Path | None = None,
          poscars_path: str | Path | None = None,
          labels: list[str] | None = None):

    labels = labels or ["id", "composition", "energy", "structure"]

    base_path = getattr(ds, "basepath", Path(os.getcwd()))
    base_path.mkdir(exist_ok=True)
    csv_file = csv_file or base_path / ds.dataset_id + '.csv'
    poscars_path = poscars_path or csv_file.with_name(csv_file.stem + 'POSCARS')
    poscars_path.mkdir(exist_ok=True)
    with open(csv_file, 'wt') as csv_h:
        csv_h.write(','.join(labels) + '\n')
        for e in ds:
            csv_h.write(",".join([_get_str_value_for_csv(e, label) for label in labels]) + '\n')
            e.structure.to(fmt="poscar", filename=poscars_path / _get_str_value_for_csv(e, "structure"), comment=e.id)






#
# def append_info(self, df, csv):
#     df_temp = pd.read_csv(csv, dtype={'Energy': float})
#     current_path = csv.parent
#     df_temp["id"] = df_temp["structure"].apply(lambda x: f"{current_path.relative_to(csv.parent)}/{x}")
#     df_temp["structure"] = df_temp["structure"].apply(lambda x: Structure.from_file(current_path /
#                                                                                poscars_parent_path / x))
#     df_temp["formula"] = df_temp["structure"].apply(lambda x: x.composition.formula)
#     df_temp["natoms"] = df_temp["structure"].apply(lambda x: len(x))
#     df_temp["energy"] = df_temp["Energy"] * (df_temp["natoms"] if source in per_atom_energy else 1)
#     df_temp["metadata"] = {"source": source}
#     df = pd.concat([df, df_temp], ignore_index=True)
#     return df


# def query(self) -> pd.DataFrame:
#     assert results_csv or results_dir_path
#     df = pd.DataFrame(columns=['id', 'energy', 'structure', 'metadata'])
#     if results_csv:
#         df = append_info(df, results_csv)
#     elif results_dir_path:
#         for current_csv in results_dir_path.rglob('*.csv'):
#             df = append_info(df, current_csv)
#     df.drop(columns=["Energy", "structure", "natoms"], inplace=True)
#     return df
#
# def write(self, ds: CrystalDataset, write_struct_info=True, columns=None):
#     columns = columns or ["id", "composition", "energy", "structure"]
#     if write_struct_info:
#         assert all([e.structure for e in ds]), "Error writing POSCARs: Structures missing"
