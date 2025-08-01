import pandas as pd
from .crystal_dataset import CrystalDataset
from .crystal_entry import CrystalEntry



def ds2df(ds: CrystalDataset) -> pd.DataFrame:
    data = dict()
    for attr in ['id', 'structure', 'energy', 'metadata']:
        data[attr] = [getattr(e, attr, None) for e in ds]
    df = pd.DataFrame(data)
    # Identify and drop all-None columns
    all_none_cols = [col for col in df.columns if df[col].isna().all()]
    df.drop(columns=all_none_cols, inplace=True)
    # Report dropped columns
    if all_none_cols:
        print(f"Dropped columns with only None values: {', '.join(all_none_cols)}")
    return df

def df2ds(df: pd.DataFrame, message=None):
    entries = [CrystalEntry(**row) for row in df.to_dict(orient='records')]
    return CrystalDataset(entries, message=message)


def cell_pos_atomtypes_from_pmg_structure(structure):

    # 1. 3×3 lattice matrix (row-major, Å)
    cell = structure.lattice.matrix  # numpy.ndarray, shape (3, 3)

    # 2. Cartesian coordinates (Å) of all sites
    positions = structure.cart_coords  # numpy.ndarray, shape (N, 3)

    # 3. Chemical symbols for each site
    atom_types = [site.specie.symbol for site in structure]  # list[str] length N

    # 4. Periodic-boundary flags along a, b, c
    pbc = [True, True, True]  # pymatgen structures are fully periodic

    return cell, positions, atom_types

