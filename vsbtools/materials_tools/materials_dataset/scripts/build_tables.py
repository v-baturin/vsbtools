from typing import Dict, Callable
from pymatgen.core import Element
from ..crystal_dataset import CrystalDataset
from .diffusion_analysis_scripts.mattergen_bridge import get_entry_fn
from ..analysis import (
    scenario_pipeline as sp,
    generation_postprocess as gp,
    phase_diagram_tools as pdt,
    symmetry_tools as st,
    atomic_environment_analysis as at_env,
    summary
)

guidance_vs_target_properties = {"environment": "compute_mean_coordination",
                                 "dominant_environment": "compute_target_share",
                                 "volume_pa": "volume_pa",
                                 "energy": "energy"}


def build_energy_vs_property_table(name_ds_dict: Dict[str, CrystalDataset],
                                   raw_stage: str = 'parse_raw',
                                   ref_stage: str = 'poll_db',
                                   target_stages: list | str | None = 'deduplicate_all',
                                   callables: Dict[str, Callable] | None = None):
    """
    For given dict builds tables containing the values of target property
    """

    if callables is None:
        callables = dict()
        target_stages = [target_stages] if isinstance(target_stages, str) else target_stages

        guidance_descr = name_ds_dict[raw_stage].metadata['batch_metadata']['guidance']
        guidance_names = list(guidance_descr.keys()) if isinstance(guidance_descr, dict) else [guidance_descr]
        assert len(guidance_names) == 1, "Only one guidance function per generation is supported"
        assert guidance_names[0] is not None, f"Invalid guidance name, check {name_ds_dict[raw_stage].base_path}"
        print(name_ds_dict[raw_stage].base_path)
        guidance_name = guidance_names[0]
        fn_name = guidance_vs_target_properties[guidance_name]
        if guidance_name in ('environment', 'dominant_environment'):
            for bond in name_ds_dict[raw_stage].metadata['batch_metadata']['guidance'][guidance_name].keys():
                if '-' not in bond:
                    continue
                params = dict(zip(('type_A', 'type_B'), [Element(e).Z for e in bond.split('-')]))
                if guidance_name == "dominant_environment":
                    params['target'] = name_ds_dict[raw_stage].metadata['batch_metadata']['guidance'][guidance_name][bond][0]
                fn = get_target_value_fn(fn_name, **params)
                callables[bond] = fn
        elif guidance_name == 'volume_pa':
            callables[guidance_name] = get_entry_fn(fn_name, **{})
            callables[guidance_name] = get_target_value_fn(fn_name, **{})

    pd_tk = pdt.PhaseDiagramTools(name_ds_dict[ref_stage])
    callables["e_hull/at"] = pd_tk.height_above_hull_pa

    stk = st.SymmetryToolkit(a_sym_prec=1e-3)
    callables["symmetry"] = stk.sym_group_symbol

    for stage in target_stages:
        summary_df = summary.collect_summary_df(name_ds_dict[stage], native_columns=("id", "composition", "energy"),
                                                callables=callables)
        summary.print_pretty_df(summary_df, name_ds_dict[stage].base_path / 'table.txt', sort_by='e_hull/at')

        summary_df.to_csv(name_ds_dict[stage].base_path / "summary.csv")

    return callables


if __name__ == "__main__":
    from pathlib import Path
    from vsbtools.materials_tools.materials_dataset.analysis.scenario_pipeline import Scenario
    from vsbtools.materials_tools.materials_dataset.scripts.diffusion_analysis_scripts.tools_for_histograms import get_guidance_generation_dirs
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
