"""
Structure sanity checks for ASE and pymatgen.

Provides:
    check_structure_sanity_ase(atoms, ...)
    check_structure_sanity_pmg(structure, ...)
"""

from __future__ import annotations

import numpy as np
from functools import lru_cache
from typing import Any, Dict, Tuple

# ---- Optional imports for each framework ------------------------------------

try:
    from ase import Atoms  # type: ignore
    from ase.data import covalent_radii as _ase_covalent_radii  # type: ignore
except ImportError:  # ASE not installed; functions will raise if used
    Atoms = Any
    _ase_covalent_radii = None

try:
    from pymatgen.core import Structure  # type: ignore
    from pymatgen.core.periodic_table import Element  # type: ignore
except ImportError:  # pymatgen not installed; functions will raise if used
    Structure = Any
    Element = None

# ---- Constants --------------------------------------------------------------

FOUR_THIRDS_PI = 4.1887902047863905  # 4/3 * pi
_EPS_VOL = 1e-6
_EPS_DET = 1e-8
_EPS_PARALLEL = 1e-10


# ---- Shared helper logic ----------------------------------------------------

def _plane_spacing_info(
    cell: np.ndarray, min_plane_spacing: float
) -> Tuple[bool, float, str]:
    """
    Compute min lattice-plane spacing and a qualitative cell degeneracy flag.
    """
    cell = np.asarray(cell, dtype=float)
    if cell.shape != (3, 3):
        raise ValueError(f"cell must be (3, 3), got {cell.shape}")

    a, b, c = cell
    det = float(np.linalg.det(cell))

    if abs(det) <= _EPS_DET:
        return False, 0.0, "Lattice vectors are (nearly) linearly dependent."

    bc = np.cross(b, c)
    ca = np.cross(c, a)
    ab = np.cross(a, b)

    n_bc = np.linalg.norm(bc)
    n_ca = np.linalg.norm(ca)
    n_ab = np.linalg.norm(ab)

    if min(n_bc, n_ca, n_ab) <= _EPS_PARALLEL:
        return False, 0.0, "At least one pair of lattice vectors is nearly parallel."

    d_a = abs(det) / n_bc
    d_b = abs(det) / n_ca
    d_c = abs(det) / n_ab
    min_d = float(min(d_a, d_b, d_c))

    if min_d >= min_plane_spacing:
        return True, min_d, "OK"
    else:
        return (
            False,
            min_d,
            f"Minimum plane spacing ({min_d:.3f} Å) < threshold ({min_plane_spacing:.3f} Å).",
        )


def _evaluate_structure(
    volume: float,
    sum_atomic_volume: float,
    cell: np.ndarray,
    min_ratio: float,
    max_ratio: float,
    min_plane_spacing: float,
) -> Tuple[bool, Dict[str, Any]]:
    """
    Shared evaluation logic: uses volume, ΣV_atom, and cell.
    """
    volume = float(volume)

    if volume <= _EPS_VOL:
        info = {
            "volume": volume,
            "sum_atomic_volume": float(sum_atomic_volume),
            "volume_ratio": np.inf,
            "packing_fraction": 0.0,
            "min_plane_spacing": 0.0,
            "ok_ratio": False,
            "ok_cell": False,
            "reason": "Cell volume is zero or numerically negligible.",
        }
        return False, info

    if sum_atomic_volume <= 0.0:
        info = {
            "volume": volume,
            "sum_atomic_volume": float(sum_atomic_volume),
            "volume_ratio": np.inf,
            "packing_fraction": 0.0,
            "min_plane_spacing": 0.0,
            "ok_ratio": False,
            "ok_cell": False,
            "reason": "Non-positive sum of covalent atomic volumes.",
        }
        return False, info

    volume_ratio = volume / sum_atomic_volume  # V_cell / ΣV_atom
    packing_fraction = sum_atomic_volume / volume  # ΣV_atom / V_cell
    ok_ratio = (min_ratio <= volume_ratio <= max_ratio)

    ok_cell, min_d, reason_cell = _plane_spacing_info(cell, min_plane_spacing)

    is_ok = ok_ratio and ok_cell

    if is_ok:
        reason = "OK"
    elif not ok_ratio and not ok_cell:
        reason = (
            f"Unreasonable volume_ratio={volume_ratio:.2f} "
            f"and degenerate/near-degenerate cell: {reason_cell}"
        )
    elif not ok_ratio:
        reason = (
            f"Unreasonable volume_ratio={volume_ratio:.2f}, "
            f"expected in [{min_ratio}, {max_ratio}]."
        )
    else:
        reason = reason_cell

    info = {
        "volume": volume,
        "sum_atomic_volume": float(sum_atomic_volume),
        "volume_ratio": float(volume_ratio),
        "packing_fraction": float(packing_fraction),
        "min_plane_spacing": float(min_d),
        "ok_ratio": bool(ok_ratio),
        "ok_cell": bool(ok_cell),
        "reason": reason,
    }

    return is_ok, info


# ---- ASE front-end ----------------------------------------------------------

def check_density_sanity_ase(
    atoms: "Atoms",
    min_ratio: float = 1.0,
    max_ratio: float = 20.0,
    min_plane_spacing: float = 0.3,
) -> Tuple[bool, Dict[str, Any]]:
    """
    Sanity check for ASE Atoms.

    Uses:
        - cell volume
        - sum of covalent-volume spheres (ASE covalent radii)
        - minimum lattice-plane spacing

    Returns (is_ok, info_dict).
    """
    if _ase_covalent_radii is None:
        raise ImportError(
            "ASE is not available or ase.data.covalent_radii could not be imported."
        )

    vol = float(atoms.get_volume())
    Z = atoms.get_atomic_numbers()
    radii = _ase_covalent_radii[Z]  # vectorized lookup
    atomic_volumes = FOUR_THIRDS_PI * (radii ** 3)
    sum_atomic_volume = float(atomic_volumes.sum())

    cell = np.array(atoms.get_cell(complete=True), dtype=float)

    return _evaluate_structure(
        volume=vol,
        sum_atomic_volume=sum_atomic_volume,
        cell=cell,
        min_ratio=min_ratio,
        max_ratio=max_ratio,
        min_plane_spacing=min_plane_spacing,
    )


# ---- pymatgen front-end -----------------------------------------------------

@lru_cache(maxsize=None)
def _pmg_radius_from_Z(Z: int) -> float:
    """Covalent radius lookup for pymatgen Element, with caching."""
    if Element is None:
        raise ImportError("pymatgen is not available.")
    el = Element.from_Z(int(Z))
    r = el.covalent_radius
    if r is None:
        r = el.atomic_radius or 1.0  # conservative fallback
    return float(r)


def check_density_sanity_pmg(
    structure: "Structure",
    min_ratio: float = 1.0,
    max_ratio: float = 20.0,
    min_plane_spacing: float = 0.3,
) -> Tuple[bool, Dict[str, Any]]:
    """
    Sanity check for pymatgen Structure.

    Uses:
        - cell volume
        - sum of covalent-volume spheres (pymatgen Element radii)
        - minimum lattice-plane spacing

    Returns (is_ok, info_dict).
    """
    if Element is None:
        raise ImportError(
            "pymatgen is not available or pymatgen.core.periodic_table.Element "
            "could not be imported."
        )

    vol = float(structure.volume)
    # For partially occupied sites, Structure.atomic_numbers uses the "major"
    # species; that is sufficient for a crude sanity check.
    Z_array = np.asarray(structure.atomic_numbers, dtype=int)
    radii = np.array([_pmg_radius_from_Z(int(Z)) for Z in Z_array], dtype=float)
    atomic_volumes = FOUR_THIRDS_PI * (radii ** 3)
    sum_atomic_volume = float(atomic_volumes.sum())

    cell = np.array(structure.lattice.matrix, dtype=float)

    return _evaluate_structure(
        volume=vol,
        sum_atomic_volume=sum_atomic_volume,
        cell=cell,
        min_ratio=min_ratio,
        max_ratio=max_ratio,
        min_plane_spacing=min_plane_spacing,
    )
