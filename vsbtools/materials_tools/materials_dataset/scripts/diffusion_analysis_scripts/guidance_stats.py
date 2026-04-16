import warnings
from typing import List, Dict, Callable, Any
from pathlib import Path
import re
import numpy as np
from pymatgen.core import Element

from .....genutils.misc import is_subtree
from ....ext_software_io.mattergen_tools.parsers import fname_friendly_serialize
from ...io.yaml_csv_poscars import read, load_yaml_recursively
from ...scripts.diffusion_analysis_scripts.mattergen_bridge import get_target_value_fn, get_loss_fn, clear_globals
from ...analysis.pipeline_legacy import LEGACY_INDEX_TO_NAME, LEGACY_NAME_TO_INDEX

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from ....visualisation_utils.formatting import cm2inch
from .pvalue_utils import get_p_value, get_two_proportion_z_test

plt.rcParams['xtick.major.pad'] = 4.
plt.rcParams['ytick.major.pad'] = 0.
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7

DEFAULT_PLOT_PRIORITY = ['reference', 'Non-guided']
NOT_IMPLEMENTED_LOSSES = ["energy"]
guidance_vs_target_properties = {"environment": "compute_mean_coordination",
                                 "dominant_environment": "compute_target_share",
                                 "volume_pa": "volume_pa",
                                 "energy": "energy"}


def callables_from_ds(ds,
                      include_losses: bool = True,
                      clear_global_state: bool = True) -> tuple[Dict[str, Callable], Dict[str, Any], str]:
    """
    Build guidance-related callables and corresponding target values from dataset metadata
    (usually parse_raw stage dataset).
    """
    if clear_global_state:
        clear_globals()

    callables = dict()
    targets = dict()
    guidance_descr = ds.metadata['batch_metadata']['guidance']
    guidance_names = list(guidance_descr.keys()) if isinstance(guidance_descr, dict) else [guidance_descr]
    assert len(guidance_names) == 1, "Only one guidance function per generation is supported"
    assert guidance_names[0] is not None, f"Invalid guidance name, check {ds.base_path}"

    guidance_name = guidance_names[0]
    fn_name = guidance_vs_target_properties[guidance_name]
    if guidance_name in ('environment', 'dominant_environment'):
        for bond, target in ds.metadata['batch_metadata']['guidance'][guidance_name].items():
            if '-' not in bond:
                continue
            callable_name = bond
            params = dict(zip(('type_A', 'type_B'), [Element(e).Z for e in bond.split('-')]))
            if isinstance(target, list) and len(target) == 2:
                params['r_cut'] = target[1]
                callable_name += f"_{target[1]}"
            if guidance_name == "dominant_environment":
                params['target'] = target[0]
                targets[callable_name] = 1
            else:
                targets[callable_name] = target[0] if isinstance(target, list) else target
            callables[callable_name] = get_target_value_fn(fn_name, **params)
    elif guidance_name == 'volume_pa':
        callables[guidance_name] = get_target_value_fn(fn_name, **{})
        targets[guidance_name] = ds.metadata['batch_metadata']['guidance'][guidance_name]

    if include_losses and guidance_name not in NOT_IMPLEMENTED_LOSSES:
        target = ds.metadata['batch_metadata']['guidance'][guidance_name]
        gl_name = f'loss_{guidance_name}_{fname_friendly_serialize(target, target.keys()) if isinstance(target, dict) else target}'
        callables[gl_name] = get_loss_fn(guidance_name, target=target)
        targets[gl_name] = 0

    return callables, targets, guidance_name


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
            if stage_desc["metadata"]["pipeline_stage"] in (0, 'parse_raw'):
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

def calculate_values(ds_dict: dict, callable_name=None, callable_params=None, fn=None, filter_max_el=True, **kwargs):
    values_dict = dict()
    if callable_params is None:
        callable_params = dict()
    if fn is None and callable_name is not None:
        predicate =  get_target_value_fn(callable_name, **callable_params)
    elif fn is not None and callable_params is None:
        predicate = fn
    else:
        raise RuntimeError("Provide either callable or function name+params (but not both)")
    for name, ds in ds_dict.items():
        if filter_max_el:
            values_dict[name] = np.array([predicate(entry) for entry in ds if len(entry.composition) == len(ds.elements)])
        else:
            values_dict[name] = np.array([predicate(entry) for entry in ds])
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


def histo_data_collection(ds_dict, callable_name, callable_params=None, auto_adjust_bins=False, n_bins=None,
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
                        group_width=0.8,
                        shift_xticks=True,
                        print_mu=False,
                        simplified_legend=False,
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

    special_label_patterns = kwargs.pop("special_label_patterns", ("Non-guided", "reference"))
    special_cmap_name = kwargs.pop("special_cmap", "berlin_r")
    other_cmap_name = kwargs.pop("other_cmap", "berlin")
    pattern_case_sensitive = kwargs.pop("pattern_case_sensitive", False)
    missing_suffix = kwargs.pop("missing_suffix", ": No data")
    bar_kwargs = kwargs.pop("bar_kwargs", {})
    legend_kwargs = kwargs.pop("legend_kwargs", {})
    ncolor = kwargs.pop("ncolor", 4)

    def sort_key_with_priority(item):
        priority = DEFAULT_PLOT_PRIORITY if custom_priority is None else custom_priority
        order = {lab: i for i, lab in enumerate(priority)}
        label = item.get("label")
        if label in order:
            return 0, order[label]
        return 1, 0

    if apply_standard_order:
        multidata = sorted(multidata, key=sort_key_with_priority)

    flags = 0 if pattern_case_sensitive else re.IGNORECASE
    pattern = "|".join(re.escape(p) for p in special_label_patterns if p)
    special_re = re.compile(pattern, flags=flags) if pattern else None

    def is_special(label):
        return special_re is not None and special_re.search(label or "") is not None

    special_indices = [i for i, d in enumerate(multidata) if is_special(d.get("label", ""))]
    special_set = set(special_indices)
    other_indices = [i for i in range(len(multidata)) if i not in special_set]

    def discrete_palette(cmap, n=ncolor):
        return [cmap(k / n) for k in range(n)]

    pal_special = discrete_palette(plt.get_cmap(special_cmap_name), ncolor)
    pal_other = discrete_palette(plt.get_cmap(other_cmap_name), ncolor)
    colors = [None] * len(multidata)
    for k, i in enumerate(special_indices):
        colors[i] = pal_special[k % ncolor]
    for k, i in enumerate(other_indices):
        colors[i] = pal_other[k % ncolor]

    valid_indices = []
    for i, d in enumerate(multidata):
        bc = d.get("bin_centers")
        ct = d.get("counts")
        if bc is None or ct is None:
            continue
        if np.asarray(bc, dtype=float).size == 0 or np.asarray(ct, dtype=float).size == 0:
            continue
        valid_indices.append(i)

    fig, ax = plt.subplots(figsize=cm2inch((8.1, 10)))
    legend_handles, legend_labels = [], []

    def add_missing_dataset_legend(i):
        label = multidata[i].get("label", f"#{i}")
        legend_handles.append(Patch(facecolor=colors[i], edgecolor='black', alpha=0.8))
        legend_labels.append(f"{label}{missing_suffix}" if missing_suffix else label)

    def finalize_and_return():
        if title:
            ax.set_title(title)
        if not legend_kwargs.get("loc", False):
            legend_kwargs["loc"] = "upper left"
        if simplified_legend:
            guid_idx = None
            if len(other_indices) == 1:
                for i, l in enumerate(legend_labels):
                    if '$=' in l:
                        guid_idx = i
            if guid_idx is not None:
                legend_labels[guid_idx] = 'Guided'

        ax.legend(legend_handles, legend_labels, fontsize='small', **legend_kwargs)
        fig.tight_layout()
        return fig, ax

    if not valid_indices:
        for i, d in enumerate(multidata):
            add_missing_dataset_legend(i)
        return finalize_and_return()

    max_center = np.inf if max_bincenter is None else float(max_bincenter)
    first_bc = np.asarray(multidata[valid_indices[0]]["bin_centers"], dtype=float)
    delta = float(first_bc[1] - first_bc[0]) if first_bc.size >= 2 else 1.0
    delta = 1.0 if delta == 0 else delta

    index_to_position = {idx: pos for pos, idx in enumerate(valid_indices)}
    n_valid = len(valid_indices)
    width = group_width * delta / max(n_valid, 1)
    offsets = np.linspace(-0.5 * group_width * delta + width / 2, 0.5 * group_width * delta - width / 2, n_valid)

    all_bincenters = np.unique(np.concatenate([
        np.asarray(multidata[i]["bin_centers"], dtype=float) for i in valid_indices
    ]))
    all_bincenters = all_bincenters[np.isfinite(all_bincenters)]
    all_bincenters.sort()
    all_bincenters = all_bincenters[all_bincenters <= max_center]

    if all_bincenters.size == 0:
        for i in range(len(multidata)):
            add_missing_dataset_legend(i)
        if "loc" not in legend_kwargs:
            legend_kwargs["loc"] = "upper left"
        return finalize_and_return()

    x_min, x_max = all_bincenters.min() - 0.5 * delta, all_bincenters.max() + 0.5 * delta
    x_line = np.linspace(x_min, x_max, 800)
    last_tick_with_plus = False

    for i, dataset in enumerate(multidata):
        label = dataset.get("label", f"#{i}")
        color = colors[i]

        if i not in index_to_position:
            add_missing_dataset_legend(i)
            continue

        pos = index_to_position[i]
        bin_centers = np.asarray(dataset["bin_centers"], dtype=float)
        counts_raw = np.asarray(dataset["counts"], dtype=float)

        if np.isfinite(max_center) and bin_centers.size > 0 and max_center <= bin_centers.max():
            mask = bin_centers <= max_center
            if np.any(mask):
                extra = counts_raw[~mask]
                bin_centers = bin_centers[mask]
                counts_raw = counts_raw[mask]
                if extra.size > 0:
                    counts_raw[-1] += extra.sum()
                    last_tick_with_plus = True
            else:
                add_missing_dataset_legend(i)
                continue

        total = counts_raw.sum()
        counts = counts_raw / total if total > 0 else counts_raw
        mu = float(np.sum(bin_centers * counts)) if counts.sum() > 0 else float("nan")

        bars = ax.bar(
            bin_centers + offsets[pos], counts, width=width,
            align="center",
            edgecolor="black", linewidth=0.5,
            color=color, alpha=0.8, zorder=2, label="_nolegend_",
            **bar_kwargs
        )
        legend_handles.append(bars[0])
        if print_mu and np.isfinite(mu):
            legend_labels.append(f"{label}, $\\mu$={mu:.2f}")
        else:
            legend_labels.append(label)

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

    if target is not None and all_bincenters.min() <= target <= all_bincenters.max():
        span = ax.axvspan(target - 0.5 * delta, target + 0.5 * delta,
                          color="#4287f5", alpha=0.3, zorder=0)
        legend_handles.append(span)
        legend_labels.append(f"Target={target}")

    are_integer_ticks = np.all(np.isclose(all_bincenters, np.round(all_bincenters)))
    xtick_labels = [f"{int(x)}" for x in all_bincenters] if are_integer_ticks else \
        [f"{x:.2f}".rstrip('0').rstrip('.') for x in all_bincenters]
    if last_tick_with_plus and xtick_labels:
        xtick_labels[-1] += "+"

    current_yticks = ax.get_yticks()
    current_yticks = current_yticks[(current_yticks <= 1.0) | (current_yticks == max(current_yticks))]
    ytick_labels = [f"{y * 100:.0f}" for y in current_yticks]
    max_tick_idx = int(np.argmax(current_yticks))
    ytick_labels[max_tick_idx] = '%'

    ax.set_xticks(all_bincenters)
    ax.set_xticklabels(xtick_labels)
    ax.set_yticks(current_yticks)
    ax.set_yticklabels(ytick_labels)

    if shift_xticks and len(all_bincenters) > 1:
        xt = ax.get_xticks()
        mid = 0.5 * (xt[:-1] + xt[1:])
        ax.set_xticks(xt)
        ax.tick_params(axis='x', which='major', length=0)
        ax.set_xticks(mid, minor=True)
        ax.tick_params(axis='x', which='minor', length=6)

    return finalize_and_return()


def plot_multi_kde(values_dict, target=None, title='', max_value=None,
                   bandwidth=None, grid_size=800,
                   kde_fill_alpha=0.2, kde_edge_lw=2.5,
                   apply_standard_order=True,
                   custom_priority=None,
                   kernel=None,
                   print_mu=False,
                   simplified_legend=False,
                   **kwargs):
    """
    Plot KDE-smoothed curves (sum of Gaussian kernels over values) for each item in values_dict.

    values_dict: dict[str, iterable[float]]
      - each key is the dataset label
      - each value is raw scalar samples
    """
    special_label_patterns = kwargs.pop("special_label_patterns", ("Non-guided", "reference"))
    special_cmap_name = kwargs.pop("special_cmap", "berlin_r")
    other_cmap_name = kwargs.pop("other_cmap", "berlin")
    pattern_case_sensitive = kwargs.pop("pattern_case_sensitive", False)
    missing_suffix = kwargs.pop("missing_suffix", ": No data")
    legend_kwargs = kwargs.pop("legend_kwargs", {})
    ncolor = kwargs.pop("ncolor", 4)
    legacy_max = kwargs.pop("max_bincenter", None)
    if legacy_max is not None:
        warnings.warn("`max_bincenter` is deprecated for KDE; use `max_value`.", DeprecationWarning)
        if max_value is None:
            max_value = legacy_max

    def gaussian_kernel(u):
        return np.exp(-0.5 * u * u) / np.sqrt(2.0 * np.pi)

    kernel_fn = gaussian_kernel if kernel is None else kernel

    multidata = [{"label": label, "values": values} for label, values in values_dict.items()]

    if apply_standard_order:
        priority = DEFAULT_PLOT_PRIORITY if custom_priority is None else custom_priority
        order = {lab: i for i, lab in enumerate(priority)}

        def sort_key(d):
            lab = d.get("label")
            if lab in order:
                return 0, order[lab]
            return 1, 0

        multidata = sorted(multidata, key=sort_key)

    flags = 0 if pattern_case_sensitive else re.IGNORECASE
    pat = "|".join(re.escape(p) for p in special_label_patterns if p)
    special_re = re.compile(pat, flags=flags) if pat else None

    def is_special(label: str) -> bool:
        if special_re is None:
            return False
        return special_re.search(label or "") is not None

    special_idx = [i for i, d in enumerate(multidata) if is_special(d.get("label", ""))]
    special_set = set(special_idx)
    other_idx = [i for i in range(len(multidata)) if i not in special_set]

    def discrete_palette(cmap, n=ncolor):
        return [cmap(k / n) for k in range(n)]

    pal_special = discrete_palette(plt.get_cmap(special_cmap_name), ncolor)
    pal_other = discrete_palette(plt.get_cmap(other_cmap_name), ncolor)

    color_list = [None] * len(multidata)
    for k, i in enumerate(special_idx):
        color_list[i] = pal_special[k % ncolor]
    for k, i in enumerate(other_idx):
        color_list[i] = pal_other[k % ncolor]

    processed = []
    max_val = np.inf if max_value is None else float(max_value)
    for i, d in enumerate(multidata):
        vals = np.asarray(d.get("values", []), dtype=float)
        vals = vals[np.isfinite(vals)]
        if np.isfinite(max_val):
            vals = vals[vals <= max_val]
        processed.append(vals)

    keep = [i for i, vals in enumerate(processed) if vals.size > 0]

    fig, ax = plt.subplots(figsize=cm2inch((8.1, 10)))
    legend_handles, legend_labels = [], []

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

    all_values = np.concatenate([processed[i] for i in keep])
    x_min = float(all_values.min())
    x_max = float(all_values.max())
    if target is not None:
        x_min = min(x_min, float(target))
        x_max = max(x_max, float(target))
    if x_max == x_min:
        x_min -= 0.5
        x_max += 0.5
    x_pad = 0.1 * (x_max - x_min)
    x_line = np.linspace(x_min - x_pad, x_max + x_pad, int(grid_size))

    for i, dataset in enumerate(multidata):
        label = dataset.get("label", f"#{i}")
        color = color_list[i]
        vals = processed[i]

        if vals.size == 0:
            legend_handles.append(Patch(facecolor=color, edgecolor='black', alpha=0.8))
            legend_labels.append(f"{label}{missing_suffix}" if missing_suffix else label)
            continue

        n = vals.size
        sigma = float(np.std(vals, ddof=1)) if n > 1 else 0.0
        q75, q25 = np.percentile(vals, [75, 25])
        iqr = float(q75 - q25)
        scale = sigma
        if iqr > 0:
            scale = min(sigma, iqr / 1.34) if sigma > 0 else iqr / 1.34
        if scale <= 0:
            data_span = float(np.max(vals) - np.min(vals))
            scale = data_span / 10.0 if data_span > 0 else 1.0
        bw = float(bandwidth) if bandwidth is not None else 0.9 * scale * (n ** (-1.0 / 5.0))
        bw = max(bw, 1e-8)

        z = (x_line[:, None] - vals[None, :]) / bw
        kernel_vals = np.asarray(kernel_fn(z), dtype=float)
        if kernel_vals.shape != z.shape:
            raise ValueError("`kernel` must support array input and return an array of the same shape.")
        density = kernel_vals.sum(axis=1) / (n * bw)

        ax.fill_between(
            x_line, density, 0,
            facecolor=color, edgecolor="none",
            alpha=kde_fill_alpha, zorder=1, label="_nolegend_"
        )
        line, = ax.plot(
            x_line, density,
            color=color, linewidth=kde_edge_lw,
            alpha=1.0, zorder=2, label="_nolegend_"
        )

        legend_handles.append(line)
        if print_mu:
            mu = float(np.mean(vals))
            legend_labels.append(f"{label}, $\\mu$={mu:.2f}")
        else:
            legend_labels.append(label)

    if target is not None:
        ymax = ax.get_ylim()[1]
        target_line = ax.vlines(target, ymin=0.0, ymax=ymax, color="#424242", alpha=0.8, linewidth=2.0, zorder=0)
        legend_handles.append(Line2D([0], [0], color=target_line.get_colors()[0], linewidth=target_line.get_linewidths()[0]))
        legend_labels.append(f"Target={target}")

    if title:
        ax.set_title(title)
    if not legend_kwargs.get("loc", False):
        legend_kwargs["loc"] = "upper left"
    if simplified_legend:
        guided_count = 0
        for i in range(len(legend_labels)):
            if '$=' in legend_labels[i]:
                legend_labels[i] = 'Guided' if guided_count == 0 else f'Guided({guided_count})'
                guided_count += 1
    ax.legend(legend_handles, legend_labels, fontsize='small', **legend_kwargs)
    fig.tight_layout()
    return fig, ax

if __name__ == "__main__":
    procesed_dir = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/PROCESSED")
    system = "Si-O"
    dirs = get_environment_gen_dirs(processed_repos_root=procesed_dir, system=system, bond='Si-O', target=4)
