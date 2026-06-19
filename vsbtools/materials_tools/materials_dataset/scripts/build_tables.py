import logging
import os
import shutil
from pathlib import Path
from typing import Dict, Callable
from pymatgen.core import Composition

from ....genutils.pareto_tools import pareto_subdataframe_indices
from ..crystal_dataset import CrystalDataset
from ..io.yaml_csv_poscars import read
from .diffusion_analysis_scripts.guidance_stats import callables_from_ds
from ..analysis import (
    phase_diagram_tools as pdt,
    symmetry_tools as st,
    summary
)


LOG = logging.getLogger(__name__)


def build_guidance_summary_table(name_ds_dict: Dict[str, CrystalDataset],
                                 raw_stage: str = 'parse_raw',
                                 target_stages: list | str | None = 'deduplicate_all',
                                 ref_stages: dict | None = None,
                                 auto_ref_stages: bool = False,
                                 callables: Dict[str, Callable] | None = None,
                                 max_pareto_front: int | None = None):
    """
    For given dict builds tables containing the values of target property
    """
    # gl_name = 'placeholder'

    target_stages = [target_stages] if isinstance(target_stages, str) else target_stages

    if callables is None:
        callables, _targets, _ = callables_from_ds(name_ds_dict[raw_stage])

    stk = st.SymmetryToolkit(a_sym_prec=1e-3)
    callables["symmetry"] = stk.sym_group_symbol

    for stage in target_stages:

        assert (ref_stages is not None and not auto_ref_stages) or (auto_ref_stages and ref_stages is None),\
            "Either ref_stages or auto_ref_stages should be specified (not both)"
        if ref_stages is None:
            ref_stage = 'poll_db_grace' if 'grace' in stage else 'poll_db'
        else:
            ref_stage = ref_stages[stage]
        pd_tk = pdt.PhaseDiagramTools(name_ds_dict[ref_stage])
        callables["e_hull/at"] = pd_tk.height_above_hull_pa

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

def stage_datasets_from_repo(repo_path: str | Path) -> Dict[str, CrystalDataset]:
    """Load all dataset manifests under a processed generation repo by stage name."""
    ds_dict: Dict[str, CrystalDataset] = {}
    for manifest_yaml in Path(repo_path).rglob("manifest.yaml"):
        dataset = read(manifest_yaml)
        ds_dict[dataset.metadata["pipeline_stage"]] = dataset
    return ds_dict


def build_guidance_summary_for_repo(
    repo_path: str | Path,
    *,
    raw_stage: str = "parse_raw",
    target_stages: list | str | None = "deduplicate_all",
    ref_stages: dict | None = None,
    auto_ref_stages: bool = False,
    callables: Dict[str, Callable] | None = None,
    max_pareto_front: int | None = None,
) -> Dict[str, Callable]:
    """Build summary tables and Pareto-front exports for one generation repo."""
    summary_callables = None if callables is None else dict(callables)
    return build_guidance_summary_table(
        stage_datasets_from_repo(repo_path),
        raw_stage=raw_stage,
        target_stages=target_stages,
        ref_stages=ref_stages,
        auto_ref_stages=auto_ref_stages,
        callables=summary_callables,
        max_pareto_front=max_pareto_front,
    )


def build_guidance_summary_for_processed_system(
    processed_system_repo: str | Path,
    *,
    raw_stage: str = "parse_raw",
    target_stages: list | str | None = "deduplicate_all",
    ref_stages: dict | None = None,
    auto_ref_stages: bool = False,
    callables: Dict[str, Callable] | None = None,
    max_pareto_front: int | None = None,
    non_guided_marker: str = "guidance_None",
) -> Dict[str, Callable]:
    """Build comparable summary artifacts for all generation repos in one system.

    Guided repos are summarized first. Their inferred callables are merged and
    then reused for the non-guided repo, which keeps the summary columns and
    Pareto-front exports aligned across the system.
    """
    processed_system_repo = Path(processed_system_repo)
    if not processed_system_repo.is_dir():
        raise NotADirectoryError(processed_system_repo)

    system_callables: Dict[str, Callable] = dict(callables or {})
    non_guided_repo: Path | None = None

    for generation_repo in sorted(processed_system_repo.iterdir()):
        if not generation_repo.is_dir():
            continue
        if non_guided_marker in generation_repo.as_posix():
            non_guided_repo = generation_repo
            continue

        repo_callables = build_guidance_summary_for_repo(
            generation_repo,
            raw_stage=raw_stage,
            target_stages=target_stages,
            ref_stages=ref_stages,
            auto_ref_stages=auto_ref_stages,
            callables=None if callables is None else system_callables,
            max_pareto_front=max_pareto_front,
        )
        system_callables.update(repo_callables)

    if non_guided_repo is None:
        LOG.warning(
            "Non-guided generation repo not found for processed system %s",
            processed_system_repo,
        )
        return system_callables

    return build_guidance_summary_for_repo(
        non_guided_repo,
        raw_stage=raw_stage,
        target_stages=target_stages,
        ref_stages=ref_stages,
        auto_ref_stages=auto_ref_stages,
        callables=system_callables,
        max_pareto_front=max_pareto_front,
    )


if __name__ == "__main__":
    from pathlib import Path
    from vsbtools.materials_tools.materials_dataset.analysis.scenario_pipeline import Scenario
    from vsbtools.materials_tools.materials_dataset.scripts.diffusion_analysis_scripts.guidance_stats import get_guidance_generation_dirs

    repo = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/"
                "MG_postprocess_pipelines/PROCESSED/Ca-Cu-P-Si/"
                "Ca-Cu-P-Si__guidance_environment_mode_huber_Cu-P_4-2.2_Cu-Cu_0-2.9__diffusion_loss_weight_1-1-False__algo_None",
                )

    ds_dict = dict()

    for manifest_yaml in repo.rglob("manifest.yaml"):
        ds = read(manifest_yaml)
        ds_dict[ds.metadata["pipeline_stage"]] = ds
        print(ds_dict.keys())
