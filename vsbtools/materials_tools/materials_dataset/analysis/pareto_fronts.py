from __future__ import annotations

import re
import warnings
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator

try:
    from ....genutils.misc import serialize_structure
except ImportError:
    from genutils.misc import serialize_structure
from ...visualisation_utils.formatting import cm2inch

try:
    import pandas as pd
except ModuleNotFoundError as _pandas_error:
    pd = None
    _PANDAS_ERROR = _pandas_error

try:
    from pymatgen.core import Composition as _Composition
except ModuleNotFoundError as _pymatgen_error:
    _Composition = None
    _PYMATGEN_ERROR = _pymatgen_error

__all__ = ["format_article_axes", "plot_pareto", "plot_pareto_fronts"]

DEFAULT_AXIS_LABELS = {
    "loss": r"$\ell$",
    "guidance_loss": r"$\ell$",
    "e_hull/at": r"$e_\mathrm{ch}/$eV",
}
TICK_LABEL_SIZE = 10
MIN_VERTICAL_TICKS = 4
MAX_VERTICAL_TICKS = 6


def format_article_axes(ax) -> None:
    ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_SIZE)

    ymin, ymax = ax.get_ylim()
    if not (np.isfinite(ymin) and np.isfinite(ymax)):
        return

    if ymax <= ymin:
        delta = max(abs(ymin) * 0.01, 1e-6)
        ymax = ymin + delta

    data_ymin = ax.dataLim.y0 if np.isfinite(ax.dataLim.y0) else ymin
    data_ymax = ax.dataLim.y1 if np.isfinite(ax.dataLim.y1) else ymax

    plot_ymax = max(ymax, data_ymax)
    near_zero_threshold = 0.06 * max(plot_ymax, 1e-12)
    data_is_near_zero = (data_ymin >= 0.0) and (data_ymin <= near_zero_threshold)

    if data_is_near_zero:
        plot_ymin = -0.06 * plot_ymax
        include_zero_tick = True
    else:
        data_span = max(data_ymax - data_ymin, 1e-6)
        plot_ymin = data_ymin - 0.06 * data_span
        include_zero_tick = False

    ax.set_ylim(plot_ymin, plot_ymax)

    ticks = None
    for nbins in range(MAX_VERTICAL_TICKS - 1, 1, -1):
        locator = MaxNLocator(
            nbins=nbins,
            min_n_ticks=MIN_VERTICAL_TICKS,
            steps=[1, 2, 2.5, 5, 10],
            prune=None,
        )
        candidate_ticks = locator.tick_values(plot_ymin, plot_ymax)
        eps = 1e-12 * max(1.0, abs(plot_ymin), abs(plot_ymax))
        candidate_ticks = candidate_ticks[
            (candidate_ticks >= plot_ymin - eps) & (candidate_ticks <= plot_ymax + eps)
        ]
        if candidate_ticks.size == 0:
            continue

        if include_zero_tick and not np.any(np.isclose(candidate_ticks, 0.0, atol=eps)):
            candidate_ticks = np.sort(np.append(candidate_ticks, 0.0))

        if MIN_VERTICAL_TICKS <= candidate_ticks.size <= MAX_VERTICAL_TICKS:
            ticks = candidate_ticks
            break

    if ticks is None:
        ticks = np.linspace(plot_ymin, plot_ymax, MIN_VERTICAL_TICKS)
        if include_zero_tick and not np.any(np.isclose(ticks, 0.0, atol=1e-12)):
            ticks = np.sort(np.append(ticks, 0.0))

    ax.set_yticks(ticks)

    y_span = abs(plot_ymax - plot_ymin)
    decimals = 3 if y_span <= 0.2 else 2
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda y, _: f"{y:.{decimals}f}".rstrip("0").rstrip("."))
    )


def _require_pandas():
    if pd is None:
        raise RuntimeError("Pareto front plotting requires pandas. Install the `materials_tools` extra.") from _PANDAS_ERROR
    return pd


def _read_dataset(manifest_path: Path):
    try:
        from ..io.yaml_csv_poscars import read
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Dataset manifest loading requires pymatgen. Install the `materials_tools` extra."
        ) from exc
    return read(manifest_path)


def _axis_label_map(prefix: str) -> dict[str, str]:
    labels = dict(DEFAULT_AXIS_LABELS)
    stripped_prefix = prefix.strip("_")
    if stripped_prefix:
        labels[f"loss_{stripped_prefix}"] = r"$\ell$"
    return labels


def _resolve_column(df: pd.DataFrame, requested: str, csv_path: Path) -> str:
    if requested in df.columns:
        return requested

    matches = [column for column in df.columns if requested in column]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ValueError(f"{csv_path}: column {requested!r} was not found.")
    raise ValueError(f"{csv_path}: column {requested!r} matched multiple columns: {matches}")


def _front_index(csv_path: Path) -> int:
    match = re.search(r"pf_(\d+)\.csv$", csv_path.name)
    if not match:
        raise ValueError(f"No trailing Pareto front index found in {csv_path.name!r}.")
    return int(match.group(1))


def _raw_stage_dataset(stage_path: Path):
    raw_stage_path = next(stage_path.parent.glob("0_*"), None)
    if raw_stage_path is None:
        return None
    manifest_path = raw_stage_path / "manifest.yaml"
    if not manifest_path.exists():
        return None
    return _read_dataset(manifest_path)


def _filter_front_df(
    df: pd.DataFrame,
    col1: str,
    col2: str,
    *,
    n_elements: int | None,
    trim_col1: float | None,
    trim_col2: float | None,
    max_loss: float | None,
    max_ehull: float | None,
) -> pd.DataFrame:
    if n_elements is not None and "composition" in df.columns:
        if _Composition is None:
            raise RuntimeError(
                "Composition filtering requires pymatgen. Install the `materials_tools` extra."
            ) from _PYMATGEN_ERROR
        df = df[df["composition"].apply(lambda value: len(_Composition(value)) == n_elements)]

    for trim_value, column in ((trim_col1, col1), (trim_col2, col2)):
        if trim_value is not None:
            df = df[df[column] <= trim_value]

    for column in (col1, col2):
        column_lower = column.lower()
        if max_loss is not None and "loss" in column_lower:
            df = df[df[column] <= max_loss]
        if max_ehull is not None and ("hull" in column_lower or column == "e_hull/at"):
            df = df[df[column] <= max_ehull]
    return df


def _default_title(stage_path: Path) -> str | None:
    raw_dataset = _raw_stage_dataset(stage_path)
    if raw_dataset is None:
        return None
    return serialize_structure(raw_dataset.metadata["batch_metadata"].get("guidance"))


def plot_pareto_fronts(
    stage_path: Path | str,
    col1: str = "loss",
    col2: str = "e_hull/at",
    trim_col1: float | None = None,
    trim_col2: float | None = None,
    n_fronts: int = 3,
    id_col: str = "id",
    show_all_background: bool = False,
    ax=None,
    max_loss: float | None = 2.5,
    max_ehull: float | None = 0.5,
    prefix: str = "",
    title: str | None = None,
    show_title: bool = False,
    article_axes: bool = False,
    **kwargs: Any,
):
    """
    Plot Pareto front CSVs from a postprocessing stage directory.

    Points whose IDs contain letters are drawn with ``+`` markers and treated
    as reference points; all other IDs are drawn as generated structures.
    """
    pandas = _require_pandas()
    stage_path = Path(stage_path)
    if show_all_background:
        warnings.warn("show_all_background is not implemented for Pareto-only CSV inputs.")

    marker_size = kwargs.pop("marker_size", 25)
    marker_alpha = kwargs.pop("marker_alpha", 0.9)
    line_width = kwargs.pop("line_width", 1.0)
    cmap = plt.get_cmap(kwargs.pop("cmap", "tab10"))
    legend = kwargs.pop("legend", True)
    grid = kwargs.pop("grid", True)

    if kwargs:
        warnings.warn(f"Unused plot_pareto_fronts kwargs: {sorted(kwargs)}")

    raw_dataset = _raw_stage_dataset(stage_path)
    n_elements = len(raw_dataset.elements) if raw_dataset is not None else None
    axis_labels = _axis_label_map(prefix)

    if ax is None:
        fig, ax = plt.subplots(figsize=cm2inch((8.1, 8.0)))
    else:
        fig = ax.figure

    fronts = {}
    resolved_col1 = None
    resolved_col2 = None
    for csv_path in sorted(stage_path.glob(f"{prefix}pf_*.csv")):
        df = pandas.read_csv(csv_path)
        x_col = _resolve_column(df, col1, csv_path)
        y_col = _resolve_column(df, col2, csv_path)

        if resolved_col1 is None:
            resolved_col1 = x_col
            resolved_col2 = y_col
        elif (x_col, y_col) != (resolved_col1, resolved_col2):
            raise ValueError(
                f"{csv_path}: resolved columns {(x_col, y_col)} differ from "
                f"previous Pareto CSV columns {(resolved_col1, resolved_col2)}."
            )

        for required_col in (x_col, y_col, id_col):
            if required_col not in df.columns:
                raise ValueError(f"{csv_path}: column {required_col!r} not found in CSV.")

        fronts[_front_index(csv_path)] = _filter_front_df(
            df,
            x_col,
            y_col,
            n_elements=n_elements,
            trim_col1=trim_col1,
            trim_col2=trim_col2,
            max_loss=max_loss,
            max_ehull=max_ehull,
        )

    if not fronts:
        warnings.warn(f"No Pareto front CSV files found under {stage_path}")
        return ax

    for i, front_i in enumerate(sorted(fronts)):
        if front_i > n_fronts:
            continue

        sub = fronts[front_i]
        if sub.empty:
            continue

        sub_sorted = sub.sort_values(resolved_col1)
        color = cmap(i % cmap.N)

        ax.plot(
            sub_sorted[resolved_col1].to_numpy(),
            sub_sorted[resolved_col2].to_numpy(),
            "-",
            linewidth=line_width,
            color=color,
            label=f"PF{front_i}",
        )

        ids_str = sub_sorted[id_col].astype(str)
        has_letters = ids_str.str.contains(r"[a-z]", case=False, regex=True)

        ref_points = sub_sorted[has_letters]
        generated_points = sub_sorted[~has_letters]

        if not ref_points.empty:
            ax.scatter(
                ref_points[resolved_col1],
                ref_points[resolved_col2],
                s=marker_size,
                alpha=marker_alpha,
                marker="+",
                color=color,
                label=None,
            )
        if not generated_points.empty:
            ax.scatter(
                generated_points[resolved_col1],
                generated_points[resolved_col2],
                s=marker_size,
                alpha=marker_alpha,
                marker="o",
                color=color,
                label=None,
            )

    x_label = axis_labels.get(resolved_col1, resolved_col1)
    y_label = axis_labels.get(resolved_col2, resolved_col2)
    if "volume" in stage_path.parent.name:
        if "loss" in resolved_col1:
            x_label += r"$/\mathrm{\AA}^3$"
        if "loss" in resolved_col2:
            y_label += r"$/\mathrm{\AA}^3$"

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if show_title:
        ax.set_title(title if title is not None else _default_title(stage_path))

    if legend:
        handles, _labels = ax.get_legend_handles_labels()
        handles += [
            Line2D([0], [0], color="k", linestyle="None", marker="+", label="Ref."),
            Line2D([0], [0], color="k", linestyle="None", marker="o", label="Gen."),
        ]
        ax.legend(handles=handles)

    if grid:
        ax.grid(True)
    if article_axes:
        format_article_axes(ax)
    ax.set_box_aspect(1)
    fig.tight_layout()
    return ax


plot_pareto = plot_pareto_fronts
