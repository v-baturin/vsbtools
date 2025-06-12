from __future__ import annotations
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.entries.computed_entries import ComputedEntry
from typing import Sequence, Iterable


def _ternary_xy(fracs: Sequence[float]) -> tuple[float, float]:
    """
    Convert barycentric fractions (bottom, left, right) → 2-D xy
    for an equilateral triangle with:
        bottom element  at (0, 0)
        right  element  at (1, 0)
        left   element  at (0.5, √3/2)
    """
    b, l, r = np.asarray(fracs, float) / np.sum(fracs)
    x = r + 0.5 * l
    y = (np.sqrt(3) / 2) * l
    return x, y


def plot_ternary_pd(
    pd: PhaseDiagram,
    ordering: Sequence[str],         # [bottom, left, right]
    *,
    show_unstable: float | None = None,
    unstable_alpha: float = 0.25,
    stable_kwargs: dict | None = None,
    unstable_kwargs: dict | None = None,
):
    """
    Draw a ternary phase diagram and (optionally) fade-in unstable entries.

    Parameters
    ----------
    pd : PhaseDiagram
        Any ternary PhaseDiagram instance.
    ordering : 3-element list/tuple/Sequence[str]
        Chemical symbols in the order *bottom – left – right*.
    show_unstable : float | None
        Hull distance (eV/atom) threshold for plotting unstable phases.
        ``None`` ⇒ do not plot unstable entries.
    unstable_alpha : float
        Transparency for the faded unstable markers (0 … 1).
    stable_kwargs / unstable_kwargs : dict
        Extra keyword arguments forwarded to ``plt.scatter``.
    """
    stable_kwargs = stable_kwargs or dict(s=60, edgecolors="k", zorder=3)
    unstable_kwargs = unstable_kwargs or dict(
        s=20, marker="s", c="grey", zorder=2, alpha=unstable_alpha
    )

    # let PDPlotter do the heavy lifting for tie-lines etc.
    plotter = PDPlotter(pd, show_unstable=False, backend="matplotlib")
    ax: plt.Axes = plotter.get_plot(ordering=ordering, fill=True)

    # --- Stable entries (circles) -------------------------------------------------
    for e in pd.stable_entries:
        x, y = _ternary_xy([e.composition.get_atomic_fraction(el) for el in ordering])
        ax.scatter(x, y, **stable_kwargs)
        ax.text(x, y, e.name, ha="center", va="center", fontsize=8)

    # --- Unstable entries (optional, squares, faded) ------------------------------
    if show_unstable is not None and show_unstable > 0:
        for e in pd.unstable_entries:
            if pd.get_e_above_hull(e) < show_unstable:
                x, y = _ternary_xy(
                    [e.composition.get_atomic_fraction(el) for el in ordering]
                )
                ax.scatter(x, y, **unstable_kwargs)

    # cosmetic axes clean-up
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()
    return ax
