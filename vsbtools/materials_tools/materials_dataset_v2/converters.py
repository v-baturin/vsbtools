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

