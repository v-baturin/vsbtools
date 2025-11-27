from typing import List, Dict
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import yaml
import torch
from ...crystal_dataset import CrystalDataset
from ...io.yaml_csv_poscars import read
from ...analysis import (
    atomic_environment_analysis as at_env,
    summary
)
from ....ext_software_io.mattergen_tools.parsers import input_parameters_to_dict
from ...scripts.diffusion_analysis_scripts.mattergen_bridge import (mattergen_cell_frac_types_fn_collection,
                                                                    mattergen_chemgraph_fn_collection,
                                                                    structure_to_tensors, entry2chemgraph)


def load_yaml_with_batch_data(yaml_fname):
    with open(yaml_fname) as ym_fid:
        ym_dict = yaml.safe_load(ym_fid.read())
    if 'batch_metadata' in ym_dict['metadata']:
        ym_dict['metadata']['batch_metadata'] = yaml.safe_load(ym_dict['metadata']['batch_metadata'])
    return ym_dict

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

def get_average_cn_gen_dirs(processed_repos_root: Path, system: str, guidance_name: str, bond: str, target: float | int | None = None):
    normalized_system = '-'.join(sorted(system.split('-')))
    search_dir = processed_repos_root / normalized_system
    gen_paths = []
    for gen_path in search_dir.glob(f"{normalized_system}*"):
        for stage_yaml in gen_path.rglob("manifest.yaml"):
            dataset_info = load_yaml_with_batch_data(stage_yaml)
            if dataset_info["metadata"]["pipeline_stage"] in [0, 'parse_raw']:
                break
        else:
            continue
        if (dataset_info["metadata"]["batch_metadata"]["guidance"] == 'None' or
                (set(dataset_info["metadata"]["batch_metadata"]["guidance"].keys()) == {guidance_name,} and
                bond in dataset_info["metadata"]["batch_metadata"]["guidance"][guidance_name].keys() and
                (not target or dataset_info["metadata"]["batch_metadata"]["guidance"][guidance_name][bond] == target))):
            gen_paths.append(gen_path)
        continue
    return gen_paths

def collect_stage_dataset_dict(gen_dirs, stage, ref_stage):
    ds_dict = dict()  #
    for stage_yml in gen_dirs[0].rglob("manifest.yaml"):
        stage_desc = load_yaml_with_batch_data(stage_yml)
        if stage_desc["metadata"]["pipeline_stage"] == ref_stage:
            ds_dict["reference"] = read(stage_yml)
    for gen_dir in gen_dirs:
        for stage_yml in gen_dir.rglob("manifest.yaml"):
            stage_desc = load_yaml_with_batch_data(stage_yml)
            if stage_desc["metadata"]["pipeline_stage"] == stage:
                ds = read(stage_yml)
            elif stage_desc["metadata"]["pipeline_stage"] in (0, 'parse_raw'):
                bmd = stage_desc["metadata"]["batch_metadata"]
                if 'guidance' in bmd and bmd['guidance'] not in (None, 'None'):
                    dlw = bmd.get("diffusion_loss_weight", ['', '', ''])
                    name = ", ".join([f"{param}={val}" for param, val in zip(["$\\kappa$", "$\\gamma$", "norm"], dlw)])
                else:
                    name = 'Non-guided'
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
            return mattergen_chemgraph_fn_collection[fn_name](x, **params).cpu().item()
    elif fn_name in mattergen_cell_frac_types_fn_collection:
        def fn(entry):
            cell, frac, types = structure_to_tensors(entry.structure)
            return mattergen_cell_frac_types_fn_collection[fn_name](cell, frac, types,
                                                                    num_atoms=torch.tensor([len(types)],
                                                                                           device = cell.device),
                                                                    **params).cpu().item()
    return fn


def calculate_values(ds_dict: dict, callable_name, callable=None, callable_params=None):
    values_dict = dict()
    if callable_params is not None and callable is None:
        callables = {callable_name: get_entry_fn(callable_name, **callable_params)}
    elif callable is not None and callable_params is None:
        callables = {callable_name: callable}
    else:
        raise RuntimeError("Provide either callable or callable params (not both)")
    for name, ds in ds_dict.items():
        summary_df = summary.collect_summary_df(ds, native_columns=("id", "composition", "energy"), callables=callables)
        df = summary_df.loc[
            summary_df['composition'].apply(len).eq(len(ds.elements))]  # We take only system having all elements
        values = pd.to_numeric(df[callable_name], errors='coerce').to_numpy()
        values_dict[name] = values[~np.isnan(values)]
    return values_dict


def values_2_histo_data(values, max_bincenter=10, norm=True, **kwargs):

    # lo = np.floor(values.min()) - 0.5
    # hi = np.ceil(values.max()) + 0.5
    bin_centers = np.arange(np.floor(values.min()), np.ceil(values.max()) + 1)
    bins = np.hstack((bin_centers - 0.5, [bin_centers[-1] + 0.5]))

    # cap by max_bincenter but keep an overflow bin
    mask = bin_centers <= max_bincenter
    bins = np.hstack((bin_centers[mask] - 0.5, [bins[-1]]))
    bin_centers = bin_centers[mask]
    counts, _ = np.histogram(values, bins=bins)
    if norm:
        counts = counts / counts.sum()
    return bin_centers, counts

def histo_data_collection(ds_dict, callable_name, callable_params=None, callable=None, **kwargs):
    """
    Returns list of dictionaries: [{'name': str, 'bin_centers': iterable[floats], }]
    """
    histo_collection = []
    histo_collection_dict = dict()
    values_dict = calculate_values(ds_dict, callable_name, callable_params=callable_params, callable=callable)
    for name, values in values_dict.items():
        hist_data = dict(zip(("bin_centers", "counts"), values_2_histo_data(values, **kwargs)))
        hist_data['label'] = name
        histo_collection.append(hist_data)
    return histo_collection



def plot_multihistogram(multidata, target=None, title='', max_bincenter=10,
                        show_gaussian=False, gaussian_fill_alpha=0.2, gaussian_edge_lw=2.5):
    """
    Plot side-by-side histograms. If show_gaussian is True, overlay a filled
    Gaussian (PDF) per dataset with the same surface (=1 after normalization),
    mean, and variance. Gaussian fills are excluded from the legend; edges are emphasized.
    """
    cmap = plt.get_cmap('tab20c')
    n_data = len(multidata)
    width = 0.8 / max(n_data, 1)
    offsets = np.linspace(-0.4 + width/2, 0.4 - width/2, n_data)

    plt.figure(figsize=(12, 6))

    all_bincenters = np.unique(np.concatenate([np.asarray(d['bin_centers']) for d in multidata]))
    x_min, x_max = all_bincenters.min() - 0.5, all_bincenters.max() + 0.5
    x_line = np.linspace(x_min, x_max, 800)

    for i, dataset in enumerate(multidata):
        bin_centers = np.asarray(dataset['bin_centers'], dtype=float)
        counts_raw  = np.asarray(dataset['counts'], dtype=float)

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

    if target is not None and np.isin(target, all_bincenters):
        plt.axvspan(target - 0.5, target + 0.5, color='orange', alpha=0.3, zorder=0, label='Target')

    plt.xticks(all_bincenters, [str(x) if x < max_bincenter else f"{max_bincenter}+" for x in all_bincenters])
    if title:
        plt.title(title)
    plt.legend(fontsize='small', loc='upper left')
    plt.tight_layout()


# def ds_dict_from_path(repo_path, system, guidance_target_bond: str, guidance_target_cn: int, guidance_type = 'environment', batch_treatments: list | None = None):
#     if batch_treatments is None:
#         batch_treatments = ['all', 'deduplicated', 'stable', 'stable_deduplicated']
#     datasets_dict = get_environment_datasets_dict(repo_path, system, guidance_type, guidance_target_bond, guidance_target_cn, batch_treatments)
#     labels = ['reference', 'non-guided', *[k for k in datasets_dict['all'].keys() if k not in ['reference', 'non-guided']]]
#     return labels, datasets_dict

# def plot_multihist_for_treatment(ds_dict, target_bond, target_cn, title='', max_bin_center=10, guidance_type = 'environment', show_gaussian=True):
#     multidata = []
#     for lbl in labels:
#         bin_centers, counts = ds_2_histo_data(ds_dict[batch_treatment][lbl], guidance_type, {"guidance_target_bond": target_bond}, max_bincenter=max_bin_center)
#         # print(lbl, bin_centers, counts, sum(counts))
#         multidata.append({'label': lbl, 'bin_centers': bin_centers, 'counts': counts})
#         # print(f"{lbl}: max_bincenter = {bin_centers.max()}")
#
#
#     plot_multihistogram(multidata, target=target_cn, max_bincenter=max_bin_center, show_gaussian=show_gaussian)
#     plt.title(title)
#     plt.show()

# _GEN_FOLDER_FILTERS = {"average_cn": get_average_cn_folders, }

if __name__ == "__main__":
    procesed_dir = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/PROCESSED")
    system = "Si-O"
    dirs = get_average_cn_gen_dirs(processed_repos_root=procesed_dir, system=system, bond='Si-O', target=4)



