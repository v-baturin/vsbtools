from typing import List, Dict
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ...crystal_dataset import CrystalDataset
from ...io.yaml_csv_poscars import read
from ...analysis import (
    atomic_environment_analysis as at_env,
    summary
)
from ....ext_software_io.mattergen_tools.parsers import input_parameters_to_dict


def graph_name_from_ds(ds: CrystalDataset):
    if ds.metadata["pipeline_stage"] == 'parse_raw':
        param_dict_guided = input_parameters_to_dict(raw=ds.metadata["batch_metadata"])
        # print(param_dict_guided)
        if 'guidance' in param_dict_guided and param_dict_guided['guidance'] is not None:
            return ", ".join([f"{param}={val}" for param, val in zip(["$\\kappa$", "$\\gamma$", "norm"], param_dict_guided.get("diffusion_loss_weight", ['','','']))])
        else:
            return "non-guided"
    elif ds.metadata["pipeline_stage"] == 2:
        return "reference"
    return None


def get_datasets_dict(processed_repos_root, system, guidance_type, guidance_target_bond, guidance_target_cn, batch_treatments):
    """
    processed_repos_root: path
    system: str
    guidance_type: str
    guidance_target_bond: str
    guidance_target_cn: int
    batch_treatments: List[str]
    """
    normalized_system = '-'.join(sorted(system.split('-')))
    system_search_path = processed_repos_root / normalized_system
    datasets_dict = {bt: dict() for bt in batch_treatments}
    datasets_dict = dict()
    for generation_path in system_search_path.glob(f"{normalized_system}*"): # go over all path and populate dictionary of datasets to plot histograms
        for manifest in generation_path.rglob("manifest.yaml"):
            ds = read(manifest)
            datasets_dict[ds.metadata["pipeline_stage"]] = ds



        raw = read(next(generation_path.glob('0_*')) / "manifest.yaml")
        # print(raw.base_path)
        guidance_label = graph_name_from_ds(raw)
        input_param_dict = input_parameters_to_dict(raw=raw.metadata["batch_metadata"])
        if input_param_dict['guidance'] is not None and \
           (guidance_type not in input_param_dict['guidance'] or \
            guidance_target_bond not in input_param_dict['guidance'][guidance_type] or \
            guidance_target_cn != input_param_dict['guidance'][guidance_type][guidance_target_bond]):
            # print(f"requested: {{{guidance_type}: {{{guidance_target_bond}: {guidance_target_cn}}}}}")
            # print(f"got:{input_param_dict['guidance']}")
            # print("Skipping\n")
            continue
        if generation_path.name.endswith('_all'):
            sym = read(next(generation_path.glob('1_*')) / "manifest.yaml") #these numbers correspond to pipeline in vsbtools.
            ref = read(next(generation_path.glob('2_*')) / "manifest.yaml")
            dedup = read(next(generation_path.glob('6_*')) / "manifest.yaml")
            datasets_dict['all'][guidance_label] = sym
            # print(f"Got {len(sym)} entries into datasets_dict['all'][{guidance_label}] from {sym.base_path}")
            datasets_dict['deduplicated'][guidance_label] = dedup
            # print(f"Got {len(dedup)} entries into datasets_dict['deduplicated'][{guidance_label}] from {dedup.base_path}")
            if 'reference' not in datasets_dict['all']:
                datasets_dict['all']['reference'] = ref
                # print(f"Got {len(ref)} entries into datasets_dict['all']['reference'] from {ref.base_path}")
                datasets_dict['deduplicated']['reference'] = ref
                # print(f"Got {len(ref)} entries into datasets_dict['deduplicated']['reference'] from {ref.base_path}")
        elif generation_path.name.endswith('_stable'):
            ref_stable = read(next(generation_path.glob('0_*')) / "manifest.yaml")
            stable = read(next(generation_path.glob('5_*')) / "manifest.yaml")
            stable_dedup = read(next(generation_path.glob('6_*')) / "manifest.yaml")
            guidance_label = graph_name_from_ds(raw)
            datasets_dict['stable'][guidance_label] = stable
            # print(f"Got {len(stable)} entries into datasets_dict['stable'][{guidance_label}] from {stable.base_path}")
            datasets_dict['stable_deduplicated'][guidance_label] = stable_dedup
            # print(f"Got {len(stable_dedup)} entries into datasets_dict['stable_dedup'][{guidance_label}] from {stable_dedup.base_path}")
            if 'reference' not in datasets_dict['stable']:
                datasets_dict['stable']['reference'] = ref_stable
                # print(f"Got {len(ref_stable)} entries into datasets_dict['stable']['reference'] from {ref_stable.base_path}")
                datasets_dict['stable_deduplicated']['reference'] = ref_stable
                # print(f"Got {len(ref_stable)} entries into datasets_dict['stable_deduplicated']['reference'] from {ref_stable.base_path}")
    return datasets_dict


# guidance_fns = {"environment": at_env.analyze_average_environment_in_entry, "dominant_environment": ''}

def get_guidance_fn(guidance_name, **params):
    labels = dict()

    if guidance_name == "environment":

        guidance_params_av_env = {type_X: at_sym for type_X, at_sym in
                                  zip(['type_A', 'type_B'], params["guidance_target_bond"].split('-'))}

        def average_env_fn(entry):
            return at_env.analyze_average_environment_in_entry(entry, **guidance_params_av_env).cpu().item()

        return average_env_fn
    elif guidance_name == "dominant_environment":
        guidance_params_av_env = {type_X: at_sym for type_X, at_sym in
                                  zip(['type_A', 'type_B'], params["guidance_target_bond"].split('-'))}
        guidance_params_target_share = guidance_params_av_env.copy()
        guidance_params_target_share["target"] = params["guidance_target_cn"]

        def average_target_share_fn(entry):
            return at_env.analyze_target_coordination_share_in_entry(entry, **guidance_params_target_share).cpu().item()

        return average_target_share_fn


def ds_2_histo_data(ds, guidance_name: str, guidance_params: dict, max_bincenter=10, norm=True):
    print(f"guidance_name = {guidance_name}, guidance_params={guidance_params}")
    callables = {guidance_name: get_guidance_fn(guidance_name, **guidance_params)}
    summary_df = summary.collect_summary_df(ds, native_columns=("id", "composition", "energy"), callables=callables)
    df = summary_df.loc[
        summary_df['composition'].apply(len).eq(len(ds.elements))]  # We take only system having all elements
    values = pd.to_numeric(df[guidance_name], errors='coerce').to_numpy()
    values = values[~np.isnan(values)]

    lo = np.floor(values.min()) - 0.5
    hi = np.ceil(values.max()) + 0.5
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


def ds_dict_from_path(repo_path, system, guidance_target_bond: str, guidance_target_cn: int, guidance_type = 'environment', batch_treatments: list | None = None):
    if batch_treatments is None:
        batch_treatments = ['all', 'deduplicated', 'stable', 'stable_deduplicated']
    datasets_dict = get_datasets_dict(repo_path, system, guidance_type, guidance_target_bond, guidance_target_cn, batch_treatments)
    labels = ['reference', 'non-guided', *[k for k in datasets_dict['all'].keys() if k not in ['reference', 'non-guided']]]
    return labels, datasets_dict

def plot_multihist_for_treatment(ds_dict, labels, target_bond, target_cn, batch_treatment = 'deduplicated', title='', max_bin_center=10, guidance_type = 'environment', show_gaussian=True):
    multidata = []
    for lbl in labels:
        bin_centers, counts = ds_2_histo_data(ds_dict[batch_treatment][lbl], guidance_type, {"guidance_target_bond": target_bond}, max_bincenter=max_bin_center)
        # print(lbl, bin_centers, counts, sum(counts))
        multidata.append({'label': lbl, 'bin_centers': bin_centers, 'counts': counts})
        # print(f"{lbl}: max_bincenter = {bin_centers.max()}")

    plot_multihistogram(multidata, target=target_cn, max_bincenter=max_bin_center, show_gaussian=show_gaussian)
    plt.title(title)
    plt.show()


if __name__ == "__main__":
    pass


