from dataclasses import dataclass
from pathlib import Path
import pandas as pd
from pymatgen.core import Structure


@dataclass
class ALIGNN_output:
    path: str | Path

    def query(self, elements) -> pd.DataFrame:
        df = pd.DataFrame(columns=['id', 'formula', 'e_total', 'natoms', 'structure', 'source'])
        for current_csv in self.path.rglob('*.csv'):
            df_temp = pd.read_csv(current_csv, dtype={'Energy': float})
            current_path = current_csv.parent
            df_temp["id"] = df_temp["File"].apply(lambda x: f"{current_path.relative_to(self.path)}/{x}")
            df_temp["structure"] = df_temp["File"].apply(lambda x: Structure.from_file(current_path / x))
            df_temp["formula"] = df_temp["structure"].apply(lambda x: x.composition.formula)
            df_temp["natoms"] = df_temp["structure"].apply(lambda x: len(x))
            df_temp["e_total"] = df_temp["Energy"] * df_temp["natoms"]
            df = pd.concat([df, df_temp], ignore_index=True)
        return df
