import os
import shutil
from typing import Dict, Callable
from pymatgen.core import Composition

from ....genutils.pareto_tools import pareto_subdataframe_indices
from ..crystal_dataset import CrystalDataset
from .diffusion_analysis_scripts.guidance_stats import callables_from_ds
from ..analysis import (
    phase_diagram_tools as pdt,
    symmetry_tools as st,
    summary
)


def build_guidance_summary_table(name_ds_dict: Dict[str, CrystalDataset],
                                 raw_stage: str = 'parse_raw',
                                 ref_stage: str = 'poll_db',
                                 target_stages: list | str | None = 'deduplicate_all',
                                 callables: Dict[str, Callable] | None = None,
                                 max_pareto_front: int | None = None):
    """
    For given dict builds tables containing the values of target property
    """
    # gl_name = 'placeholder'

    if callables is None:
        callables, _targets, _ = callables_from_ds(name_ds_dict[raw_stage])
        target_stages = [target_stages] if isinstance(target_stages, str) else target_stages

    pd_tk = pdt.PhaseDiagramTools(name_ds_dict[ref_stage])
    callables["e_hull/at"] = pd_tk.height_above_hull_pa

    stk = st.SymmetryToolkit(a_sym_prec=1e-3)
    callables["symmetry"] = stk.sym_group_symbol

    for stage in target_stages:
        summary_df = summary.collect_summary_df(name_ds_dict[stage], native_columns=("id", "composition", "energy"),
                                                callables=callables)
        summary.print_pretty_df(summary_df, name_ds_dict[stage].base_path / 'table.txt', sort_by='e_hull/at')

        summary_df.to_csv(name_ds_dict[stage].base_path / "summary.csv")

        losses = [call_name for call_name in callables.keys() if 'loss' in call_name]

        losses_prefixes = ['']
        if len(losses) > 1:
            losses_prefixes = [loss.replace('loss_', '') + '_' for loss in losses]

        if max_pareto_front is not None:
            for j, gl_name in enumerate(losses):
                max_el_db = summary_df[summary_df['composition'].apply(lambda x: len(Composition(x))) == len(name_ds_dict[stage].elements)].sort_values('e_hull/at')
                fronts_idx, rank = pareto_subdataframe_indices(max_el_db, ['e_hull/at', gl_name],  max_pareto_front)
                for i, idx_pf in enumerate(fronts_idx):
                    max_el_db.iloc[idx_pf].to_csv(name_ds_dict[stage].base_path / f"{losses_prefixes[j]}pf_{i+1}.csv")
                    summary.print_pretty_df(max_el_db.iloc[idx_pf], name_ds_dict[stage].base_path / f"{losses_prefixes[j]}pf_{i+1}_table.txt",
                                            sort_by='e_hull/at')
                    pf_dir = name_ds_dict[stage].base_path / f"{losses_prefixes[j]}pf_{i+1}"
                    os.makedirs(pf_dir, exist_ok=True)
                    for k, idx in enumerate(idx_pf):
                        shutil.copy2(name_ds_dict[stage].base_path / "POSCARS" / f"{max_el_db.iloc[idx].id}POSCAR",
                                     pf_dir / f"{k+1}_{max_el_db.iloc[idx].id}POSCAR")



    return callables


if __name__ == "__main__":
    from pathlib import Path
    from vsbtools.materials_tools.materials_dataset.analysis.scenario_pipeline import Scenario
    from vsbtools.materials_tools.materials_dataset.scripts.diffusion_analysis_scripts.guidance_stats import get_guidance_generation_dirs
    from vsbtools.materials_tools.materials_dataset.io.yaml_csv_poscars import read

    repo = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/"
                "MG_postprocess_pipelines/PROCESSED/Ca-Cu-P-Si/"
                "Ca-Cu-P-Si__guidance_environment_mode_huber_Cu-P_4-2.2_Cu-Cu_0-2.9__diffusion_loss_weight_1-1-False__algo_None",
                )

    ds_dict = dict()

    for manifest_yaml in repo.rglob("manifest.yaml"):
        ds = read(manifest_yaml)
        ds_dict[ds.metadata["pipeline_stage"]] = ds
        print(ds_dict.keys())
