import warnings
from typing import List, Dict
from pathlib import Path
import re
import numpy as np

from .....genutils.misc import is_subtree
from ...io.yaml_csv_poscars import read, load_yaml_recursively
from ...scripts.diffusion_analysis_scripts.mattergen_bridge import get_target_value_fn
from ...analysis.pipeline_legacy import LEGACY_INDEX_TO_NAME, LEGACY_NAME_TO_INDEX

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from ....visualisation_utils.formatting import cm2inch

plt.rcParams['xtick.major.pad'] = '0'
plt.rcParams['ytick.major.pad'] = '0.'
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7

DEFAULT_PLOT_PRIORITY = ['reference', 'Non-guided']


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
                    name = ", ".join([f"{param}={val}" for param, val in zip(["$k$", "$g$", "norm"], dlw)])
                    if add_guid_descr:
                        name += (", " + "_".join(f"{k}_{v}" for k, v in bmd['guidance'].items()))
                else:
                    name = 'Non-guided'
        if name is not None:
            ds_dict[name] = ds
    return ds_dict



def calculate_values(ds_dict: dict, callable_name, fn=None, callable_params=None):
    values_dict = dict()
    if callable_params is not None and fn is None:
        predicate =  get_target_value_fn(callable_name, **callable_params)
    elif fn is not None and callable_params is None:
        predicate = fn
    else:
        raise RuntimeError("Provide either callable or callable params (not both)")
    for name, ds in ds_dict.items():
        ds_all_elements = [entry for entry in ds if len(entry.composition) == len(ds.elements)] # only max number of elements
        values_dict[name] = np.array([predicate(entry) for entry in ds_all_elements])
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
                        show_gaussian=False, gaussian_fill_alpha=0.2, gaussian_edge_lw=2.5,
                        apply_standard_order=True,
                        custom_priority=None,
                        **kwargs):
    """
    Plot side-by-side histograms for multidata:
    multidata:: list of dicts:
      {'label': str, 'bin_centers': iterable[float] | None, 'counts': iterable[float] | None}

    Missing datasets (bin_centers or counts is None) are still shown in the legend
    using proxy handles with the intended color.

    Color grouping:
      Labels matching any of special_label_patterns use special_cmap; others use other_cmap.
      (Configured via kwargs; see below.)
    """

    # --- optional config via kwargs (doesn't break existing calls) ---
    special_label_patterns = kwargs.pop("special_label_patterns", ("Non-guided", "reference"))
    special_cmap_name      = kwargs.pop("special_cmap", "copper")
    other_cmap_name        = kwargs.pop("other_cmap", "ocean")
    pattern_case_sensitive = kwargs.pop("pattern_case_sensitive", False)
    missing_suffix         = kwargs.pop("missing_suffix", ": No data")
    bar_kwargs             = kwargs.pop("bar_kwargs", {})       # extra kwargs passed to ax.bar
    legend_kwargs          = kwargs.pop("legend_kwargs", {})    # extra kwargs passed to ax.legend

    # --- ordering ---
    if apply_standard_order:
        priority = DEFAULT_PLOT_PRIORITY if custom_priority is None else custom_priority
        order = {lab: i for i, lab in enumerate(priority)}

        def sort_key(d):
            lab = d.get("label")
            if lab in order:
                return 0, order[lab]
            return 1, 0

        multidata = sorted(multidata, key=sort_key)

    n_all = len(multidata)

    # --- classify labels into color groups ---
    flags = 0 if pattern_case_sensitive else re.IGNORECASE
    pat = "|".join(re.escape(p) for p in special_label_patterns if p)
    special_re = re.compile(pat, flags=flags) if pat else None

    def is_special(label: str) -> bool:
        if special_re is None:
            return False
        return special_re.search(label or "") is not None

    special_idx = [i for i, d in enumerate(multidata) if is_special(d.get("label", ""))]
    special_set = set(special_idx)
    other_idx = [i for i in range(n_all) if i not in special_set]

    NCOLOR = kwargs.pop("ncolor", 4)

    cmap_special = plt.get_cmap(special_cmap_name)
    cmap_other = plt.get_cmap(other_cmap_name)

    def discrete_palette(cmap, n=NCOLOR):
        # sample midpoints: 0.05, 0.15, ..., 0.95 (for n=10) -> avoids extremes
        return [cmap((n - k - 0.5) / n) for k in range(n)]

    pal_special = discrete_palette(cmap_special, NCOLOR)
    pal_other = discrete_palette(cmap_other, NCOLOR)

    color_list = [None] * len(multidata)

    # sequential assignment within each group (wrap around if group > NCOLOR)
    for k, i in enumerate(special_idx):
        color_list[i] = pal_special[k % NCOLOR]

    for k, i in enumerate(other_idx):
        color_list[i] = pal_other[k % NCOLOR]

    # --- identify non-empty datasets ---
    keep = []
    for i, d in enumerate(multidata):
        bc = d.get("bin_centers")
        ct = d.get("counts")
        if bc is None or ct is None:
            continue
        bc = np.asarray(bc, dtype=float)
        ct = np.asarray(ct, dtype=float)
        if bc.size == 0 or ct.size == 0:
            continue
        keep.append(i)

    fig, ax = plt.subplots(figsize=cm2inch((8.1, 10)))

    legend_handles, legend_labels = [], []

    # If nothing to plot: still show legend with correct colors
    if len(keep) == 0:
        for i, d in enumerate(multidata):
            label = d.get("label", f"#{i}")
            color = color_list[i]
            legend_handles.append(Patch(facecolor=color, edgecolor='black', alpha=0.8))
            legend_labels.append(f"{label}{missing_suffix}" if missing_suffix else label)
        if title:
            ax.set_title(title)
        if not legend_kwargs.get("loc", False):
            legend_kwargs["loc"] = "upper left"
        ax.legend(legend_handles, legend_labels, fontsize='small', **legend_kwargs)
        fig.tight_layout()
        return fig, ax

    max_bincenter_val = np.inf if max_bincenter is None else float(max_bincenter)

    # Geometry from the first non-empty dataset
    example_bc = np.asarray(multidata[keep[0]]["bin_centers"], dtype=float)
    if example_bc.size >= 2:
        delta = float(example_bc[1] - example_bc[0])
        if delta == 0:
            delta = 1.0
    else:
        delta = 1.0

    kept_pos = {idx: pos for pos, idx in enumerate(keep)}
    n_kept = len(keep)

    width = 0.8 * delta / max(n_kept, 1)
    offsets = np.linspace(-0.4 * delta + width/2, 0.4 * delta - width/2, n_kept)

    # Union of bin centers (for axis range/ticks/gaussian x grid)
    all_bincenters = np.unique(np.concatenate([
        np.asarray(multidata[i]["bin_centers"], dtype=float) for i in keep
    ]))
    all_bincenters = all_bincenters[np.isfinite(all_bincenters)]
    all_bincenters.sort()
    all_bincenters = all_bincenters[all_bincenters <= max_bincenter_val]

    if all_bincenters.size == 0:
        # Still show legend; nothing sensible to place on x-axis.
        for i, d in enumerate(multidata):
            label = d.get("label", f"#{i}")
            color = color_list[i]
            legend_handles.append(Patch(facecolor=color, edgecolor='black', alpha=0.8))
            legend_labels.append(f"{label}{missing_suffix}" if missing_suffix else label)
        if title:
            ax.set_title(title)
        ax.legend(legend_handles, legend_labels, fontsize='small', loc='upper left', **legend_kwargs)
        fig.tight_layout()
        return fig, ax

    x_min, x_max = all_bincenters.min() - 0.5 * delta, all_bincenters.max() + 0.5 * delta
    x_line = np.linspace(x_min, x_max, 800)

    last_tick_with_plus = False

    # --- plot in the same order as multidata so legend order is stable ---
    for i, dataset in enumerate(multidata):
        label = dataset.get("label", f"#{i}")
        color = color_list[i]

        if i not in kept_pos:
            # proxy legend entry for missing dataset
            legend_handles.append(Patch(facecolor=color, edgecolor='black', alpha=0.8))
            legend_labels.append(f"{label}{missing_suffix}" if missing_suffix else label)
            continue

        pos = kept_pos[i]

        bin_centers = np.asarray(dataset["bin_centers"], dtype=float)
        counts_raw  = np.asarray(dataset["counts"], dtype=float)

        # Trim & accumulate tail into last bin if max_bincenter is active
        if np.isfinite(max_bincenter_val) and bin_centers.size > 0 and max_bincenter_val <= bin_centers.max():
            mask = bin_centers <= max_bincenter_val
            if np.any(mask):
                extra = counts_raw[~mask]
                bin_centers = bin_centers[mask]
                counts_raw = counts_raw[mask]
                if extra.size > 0:
                    counts_raw[-1] += extra.sum()
                    last_tick_with_plus = True
            # if nothing passes mask, treat as effectively missing-on-plot but keep legend entry
            else:
                legend_handles.append(Patch(facecolor=color, edgecolor='black', alpha=0.8))
                legend_labels.append(f"{label}{missing_suffix}" if missing_suffix else label)
                continue

        total = counts_raw.sum()
        counts = counts_raw / total if total > 0 else counts_raw

        mu = float(np.sum(bin_centers * counts)) if counts.sum() > 0 else float("nan")

        bars = ax.bar(
            bin_centers + offsets[pos], counts, width=width,
            align="center",
            edgecolor="black", linewidth=0.5,
            color=color, alpha=0.8, zorder=2,
            **bar_kwargs
        )

        legend_handles.append(bars[0])
        if np.isfinite(mu):
            legend_labels.append(f"{label}, $\\mu$={mu:.2f}")
        else:
            legend_labels.append(label)

        # Gaussian overlay (excluded from legend)
        if show_gaussian and counts.sum() > 0:
            var = float(np.sum((bin_centers - mu) ** 2 * counts))
            sigma = float(np.sqrt(max(var, 0.0)))

            if sigma > 0:
                y = (1.0 / (sigma * np.sqrt(2*np.pi))) * np.exp(-0.5 * ((x_line - mu) / sigma) ** 2)
                ax.fill_between(
                    x_line, y, 0,
                    facecolor=color, edgecolor="none",
                    alpha=gaussian_fill_alpha, zorder=1, label="_nolegend_"
                )
                ax.plot(
                    x_line, y,
                    color=color, linewidth=gaussian_edge_lw,
                    alpha=1.0, zorder=3, label="_nolegend_"
                )
            else:
                eps = 0.15
                ax.axvspan(mu - eps, mu + eps, facecolor=color, edgecolor="none",
                           alpha=gaussian_fill_alpha, zorder=1, label="_nolegend_")
                ax.axvline(mu, color=color, linewidth=gaussian_edge_lw,
                           alpha=1.0, zorder=3, label="_nolegend_")

    # Target span (add to legend explicitly because we use manual handles)
    if target is not None and all_bincenters.min() <= target <= all_bincenters.max():
        span = ax.axvspan(target - 0.5 * delta, target + 0.5 * delta,
                          color="#4287f5", alpha=0.3, zorder=0)
        legend_handles.append(span)
        legend_labels.append(f"Target={target}")

# ------ setting ticks
    xtick_labels = [f"{x:.2f}".rstrip('0').rstrip('.') for x in all_bincenters] if any([x != int (x) for x in all_bincenters]) \
        else ([f"{int(x)}" for x in all_bincenters])
    if last_tick_with_plus and xtick_labels:
        xtick_labels[-1] += "+"
    current_yticks = ax.get_yticks()
    current_yticks = current_yticks[(current_yticks<=1.) | (current_yticks == max(current_yticks))]
    ytick_labels = [f"{y * 100:.0f}" for y in current_yticks]
    print(current_yticks)
    print(ytick_labels)
    ytick_labels[np.where(current_yticks == current_yticks.max())[0][0]] = '%'
    ax.set_xticks(all_bincenters)
    ax.set_xticklabels(xtick_labels)
    ax.set_yticks(current_yticks)
    ax.set_yticklabels(ytick_labels)

# ------ setting title
    if title:
        ax.set_title(title)
    if not legend_kwargs.get("loc", False):
        legend_kwargs["loc"] = "upper left"
    ax.legend(legend_handles, legend_labels, fontsize='small', **legend_kwargs)
    fig.tight_layout()
    return fig, ax

if __name__ == "__main__":
    procesed_dir = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/PROCESSED")
    system = "Si-O"
    dirs = get_environment_gen_dirs(processed_repos_root=procesed_dir, system=system, bond='Si-O', target=4)



