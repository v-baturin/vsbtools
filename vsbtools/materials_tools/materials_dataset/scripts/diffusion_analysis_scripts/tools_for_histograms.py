import warnings
from typing import List, Dict
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import yaml
import torch
from .....genutils.misc import is_subtree
from ...crystal_dataset import CrystalDataset
from ...io.yaml_csv_poscars import read, load_yaml_recursively
from ....ext_software_io.mattergen_tools.parsers import input_parameters_to_dict
from ...scripts.diffusion_analysis_scripts.mattergen_bridge import (mattergen_cell_frac_types_fn_collection,
                                                                    mattergen_chemgraph_fn_collection,
                                                                    structure_to_tensors, entry2chemgraph)
from ...analysis.pipeline_legacy import LEGACY_INDEX_TO_NAME, LEGACY_NAME_TO_INDEX


def graph_name_from_ds(ds: CrystalDataset):
    if ds.metadata["pipeline_stage"] in [0, 'parse_raw']:
        param_dict_guided = input_parameters_to_dict(raw=ds.metadata["batch_metadata"])
        # print(param_dict_guided)
        if 'guidance' in param_dict_guided and param_dict_guided['guidance'] is not None:
            return ", ".join([f"{param}={val}" for param, val in zip(["$\\kappa$", "$\\gamma$", "norm"], param_dict_guided.get("diffusion_loss_weight", ['', '', '']))])
        else:
            return "non-guided"
    elif ds.metadata["pipeline_stage"] in [2, "poll_db"]:
        return "reference"
    return None


def get_guidance_generation_dirs(processed_repos_root: Path, system: str, guidance_sub_dict: Dict, include_non_guided: bool = True):
    """
    Get a list of directories containing postprocessing pipelines corresponding to a mattergen-generated structures
    with guidance corresponding to the guidance_sub_dict argument.

    Args:
        processed_repos_root: pathlib.Path -- root of processed repositories
        system: str -- dash('-') -separated elements (in any order), e.g. 'Ni-B-H'
        guidance_sub_dict: Dict[str, str|Dict] -- sub-dictionary looked upon in the dataset_info["metadata"]["batch_metadata"]
    Kwargs:
        include_non_guided: bool -- whether the non-guided generation is included
    sub_dict_examples:
        environment (mean coordination number): {'guidance': {'environment': {'mode': 'huber', 'Si-O': 6}}}
        volume_pa: {'guidance': {'volume_pa': 6.8}}
    """
    normalized_system = '-'.join(sorted(system.split('-')))
    search_dir = processed_repos_root / normalized_system
    gen_paths = []
    for gen_path in search_dir.glob(f"{normalized_system}*"):
        for stage_yaml in gen_path.rglob("manifest.yaml"):
            dataset_info = load_yaml_recursively(stage_yaml)
            if dataset_info["metadata"]["pipeline_stage"] in [0, 'parse_raw']:
                break
        else:
            continue
        if (dataset_info["metadata"]["batch_metadata"]["guidance"] == 'None' and include_non_guided) or \
            is_subtree(dataset_info["metadata"]["batch_metadata"], guidance_sub_dict):
            gen_paths.append(gen_path)
        continue
    return gen_paths


def get_volume_pa_gen_dirs(processed_repos_root: Path, system: str, guidance_name: str, target: float | int | None = None):
    guidance_sub_dict = {'guidance': {guidance_name: target}}
    return get_guidance_generation_dirs(processed_repos_root, system, guidance_sub_dict=guidance_sub_dict)


def get_environment_gen_dirs(processed_repos_root: Path, system: str, guidance_name: str,
                             bond: str = None, target: float | int | None = None):
    guidance_sub_dict = {'guidance': {guidance_name: {bond: target}}}
    return get_guidance_generation_dirs(processed_repos_root, system, guidance_sub_dict=guidance_sub_dict)


def collect_stage_dataset_dict(gen_dirs, stage, ref_stage, add_guid_descr=False):
    ds_dict = dict()
    ds = None
    for stage_yml in gen_dirs[0].rglob("manifest.yaml"):
        stage_desc = load_yaml_recursively(stage_yml)
        if stage_desc["metadata"]["pipeline_stage"] in (ref_stage, LEGACY_NAME_TO_INDEX.get(ref_stage, None)):
            ds_dict["reference"] = read(stage_yml)
    for gen_dir in gen_dirs:
        name = None
        for stage_yml in gen_dir.rglob("manifest.yaml"):
            stage_desc = load_yaml_recursively(stage_yml)
            if stage_desc["metadata"]["pipeline_stage"] in (stage, LEGACY_NAME_TO_INDEX.get(stage, None)):
                ds = read(stage_yml)
            elif stage_desc["metadata"]["pipeline_stage"] in (0, 'parse_raw'):
                bmd = stage_desc["metadata"]["batch_metadata"]
                if 'guidance' in bmd and bmd['guidance'] not in (None, 'None'):
                    dlw = bmd.get("diffusion_loss_weight", ['', '', ''])
                    name = ", ".join([f"{param}={val}" for param, val in zip(["$\\kappa$", "$\\gamma$", "norm"], dlw)])
                    if add_guid_descr:
                        name += (", " + "_".join(f"{k}_{v}" for k, v in bmd['guidance'].items()))
                else:
                    name = 'Non-guided'
        if name is not None:
            ds_dict[name] = ds
    return ds_dict


def get_entry_fn(fn_name, **params):
    fn = lambda x: None
    if fn_name in ("environment", "average_cn"):
        fn_params = {type_X: at_sym for type_X, at_sym in
                                  zip(['type_A', 'type_B'], params.pop("guidance_target_bond").split('-'))}
        def fn(entry):
            return at_env.analyze_average_environment_in_entry(entry, **{**fn_params, **fn_params}).cpu().item()
    elif fn_name == "dominant_environment":
        fn_params = {type_X: at_sym for type_X, at_sym in
                                  zip(['type_A', 'type_B'], params["guidance_target_bond"].split('-'))}
        fn_params_cp = fn_params.copy()
        fn_params_cp["target"] = params["guidance_target_cn"]
        def fn(entry):
            return at_env.analyze_target_coordination_share_in_entry(entry, **fn_params_cp).cpu().item()

    elif fn_name in mattergen_chemgraph_fn_collection:
        def fn(entry):
            x = entry2chemgraph(entry)
            return mattergen_chemgraph_fn_collection[fn_name](x, t=None, **params).cpu().detach().numpy()[0]  # x is a batch normally, so to have a result of entry we need [0]
    elif fn_name in mattergen_cell_frac_types_fn_collection:
        def fn(entry):
            cell, frac, types = structure_to_tensors(entry.structure)
            return mattergen_cell_frac_types_fn_collection[fn_name](cell, frac, types,
                                                                    num_atoms=torch.tensor([len(types)],
                                                                                           device = cell.device),
                                                                    **params).cpu().detach().numpy()
    return fn


def calculate_values(ds_dict: dict, callable_name, fn=None, callable_params=None):
    values_dict = dict()
    if callable_params is not None and fn is None:
        callables = {callable_name: get_entry_fn(callable_name, **callable_params)}
    elif fn is not None and callable_params is None:
        callables = {callable_name: fn}
    else:
        raise RuntimeError("Provide either callable or callable params (not both)")
    for name, ds in ds_dict.items():
        summary_df = summary.collect_summary_df(ds, native_columns=("id", "composition", "energy", "metadata.duplicates"), callables=callables)
        summary_df = summary_df.loc[
            summary_df['composition'].apply(len).eq(len(ds.elements))].copy()  # We take only system having all elements
        # summary.print_pretty_df(summary_df, f"{ds.base_path.name}._table.txt")
        values = pd.to_numeric(summary_df[callable_name], errors='coerce').to_numpy().astype(float)
        values_dict[name] = values[~np.isnan(values)]
    return values_dict


def values_2_histo_data(values, name=None, bin_centers=None, bins = None, integer_bins = True, n_bins=10, norm=True, **kwargs):

    if bin_centers is None:
        if bins is not None:
            warnings.warn("bin_centers will be calculated automatically, provided bins will be discarded")
            bins = None
        if integer_bins:
            bin_centers = np.arange(np.floor(values.min()), np.ceil(values.max()) + 1)
        else:
            bin_centers = np.linspace(values.min(), values.max(), n_bins)

    if bins is None:
        half_shift = 0.5 * (bin_centers[1] - bin_centers[0])
        bins = np.hstack((bin_centers - half_shift, [bin_centers[-1] + half_shift]))
        if bins[0] > values.min():
            warnings.warn(f"Dataset {name}: values outside leftmost bin edge ({values[values < bins[0]]})")
        if bins[-1] < values.max():
            warnings.warn(f"Dataset {name}: values outside rightmost bin edge ({values[values > bins[-1]]})")

    # cap by max_bincenter but keep an overflow bin
    counts, _ = np.histogram(values, bins=bins)
    if norm:
        counts = counts / counts.sum()
    return bin_centers, counts


def histo_data_collection(ds_dict, callable_name, callable_params=None, fn=None, auto_adjust_bins=False, n_bins=None,
                          **kwargs):
    """
    Returns list of dictionaries: [{'label': str, 'bin_centers': iterable[floats], 'counts': iterable[floats]}]
    """
    histo_collection = []
    values_dict = calculate_values(ds_dict, callable_name, callable_params=callable_params)
    if auto_adjust_bins:
        extremities = np.unique(np.concatenate([np.asarray([v.min(), v.max()]) for v in values_dict.values()]))
        bin_centers = np.linspace(extremities.min(), extremities.max(), n_bins)
        shift = (bin_centers[1] - bin_centers[0])
        bins = np.hstack((bin_centers - shift, [bin_centers[-1] + shift]))
        kwargs['bins'] = bins
        kwargs['bin_centers'] = bin_centers
    for name, values in values_dict.items():
        if len(values) > 0:
            hist_data = dict(zip(("bin_centers", "counts"), values_2_histo_data(values, name, **kwargs)))
            hist_data['label'] = name
        else:
            hist_data = {'label': name, "bin_centers": None, "counts": None}
        histo_collection.append(hist_data)
    return histo_collection



def plot_multihistogram(multidata, target=None, title='', max_bincenter=None,
                        show_gaussian=False, gaussian_fill_alpha=0.2, gaussian_edge_lw=2.5, **kwargs):
    """
    Plot side-by-side histograms for multidata:
    multidata::  list of dictionaries: [{'label': str, 'bin_centers': iterable[floats], 'counts': iterable[floats]}]
    If show_gaussian is True, overlay a filled
    Gaussian (PDF) per dataset with the same surface (=1 after normalization),
    mean, and variance. Gaussian fills are excluded from the legend; edges are emphasized.
    """
    keep = [i for i, d in enumerate(multidata) if d["bin_centers"] is not None]
    kept_multidata = [multidata[kpt] for kpt in keep]
    max_bincenter = max_bincenter or np.inf
    example_data = multidata[keep[0]]['bin_centers']
    delta = example_data[1] - example_data[0]

    cmap = plt.get_cmap('tab20c')
    n_data = len(keep)
    width = 0.8 * delta / max(n_data, 1)
    offsets = np.linspace(-0.4 * delta + width/2, 0.4 * delta - width/2, n_data)

    plt.figure(figsize=(12, 6))

    all_bincenters = np.unique(np.concatenate([np.asarray(d['bin_centers']) for d in kept_multidata]))

    if max_bincenter is not None:
        all_bincenters = all_bincenters[all_bincenters <= max_bincenter]

    x_min, x_max = all_bincenters.min() - 0.5 * delta, all_bincenters.max() + 0.5 * delta
    x_line = np.linspace(x_min, x_max, 800)

    last_tick_with_plus = False

    for i, dataset in enumerate(multidata):
        if i not in keep:
            plt.bar([], [], color="k", label=f"{dataset['label']}: No data")

    for i, dataset in enumerate(kept_multidata):
        bin_centers = np.asarray(dataset['bin_centers'], dtype=float)
        counts_raw  = np.asarray(dataset['counts'], dtype=float)

        if max_bincenter is not None and max_bincenter <= bin_centers.max():
            mask = bin_centers <= max_bincenter
            bin_centers = bin_centers[mask]
            extra_counts_raw = counts_raw[~mask]
            counts_raw = counts_raw[mask]
            counts_raw[-1] += extra_counts_raw.sum()
            last_tick_with_plus = True


        total = counts_raw.sum()
        counts = counts_raw / total if total > 0 else counts_raw

        mu = np.sum(bin_centers * counts)

        color = cmap(i / max(n_data, 1))
        # Bars: keep default-ish edges (focus is on Gaussian edges)
        plt.bar(bin_centers + offsets[i], counts, width=width, label=f"{dataset['label']}, $\\mu=${mu:.2f}",
                align='center', edgecolor='black', linewidth=0.5, color=color, alpha=0.8, zorder=2)

        if show_gaussian and counts.sum() > 0:
            mu = np.sum(bin_centers * counts)
            var = np.sum((bin_centers - mu)**2 * counts)
            sigma = np.sqrt(max(var, 0.0))

            if sigma > 0:
                y = (1.0 / (sigma * np.sqrt(2*np.pi))) * np.exp(-0.5 * ((x_line - mu)/sigma)**2)
                # Filled area (no legend), then a thick outline for visibility
                plt.fill_between(x_line, y, 0, facecolor=color, edgecolor='none',
                                 alpha=gaussian_fill_alpha, zorder=1, label='_nolegend_')
                plt.plot(x_line, y, color=color, linewidth=gaussian_edge_lw,
                         alpha=1.0, zorder=3, label='_nolegend_')
            else:
                # Degenerate distribution: translucent slab + thick central line
                eps = 0.15
                plt.axvspan(mu - eps, mu + eps, facecolor=color, edgecolor='none',
                            alpha=gaussian_fill_alpha, zorder=1, label='_nolegend_')
                plt.axvline(mu, color=color, linewidth=gaussian_edge_lw,
                            alpha=1.0, zorder=3, label='_nolegend_')

    if target is not None and all_bincenters.min() <= target <= all_bincenters.max():
        plt.axvspan(target - 0.5 * delta, target + 0.5 * delta, color='orange', alpha=0.3, zorder=0, label=f'Target={target}')
    tick_labels = [f"{x:.2f}" for x in all_bincenters]
    if last_tick_with_plus:
        tick_labels[-1] += "+"
    plt.xticks(all_bincenters, tick_labels)
    if title:
        plt.title(title)
    plt.legend(fontsize='small', loc='upper left')
    plt.tight_layout()

if __name__ == "__main__":
    procesed_dir = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/PROCESSED")
    system = "Si-O"
    dirs = get_environment_gen_dirs(processed_repos_root=procesed_dir, system=system, bond='Si-O', target=4)



