
from __future__ import annotations

from typing import Callable, Iterable
import numpy as np
import torch
from ase.data import covalent_radii, atomic_numbers as _Z

# Device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Standard, fully differentiable distance kernels g(d)
_STANDARD_KERNELS: dict[str, Callable] = {
    "gaussian": lambda dist, sigma=1.0, **kwargs: torch.exp(- (dist / sigma).pow(2)),
    "sigmoid":  lambda dist, alpha=8.0, r_cut=None, **kwargs: torch.sigmoid(alpha * (r_cut - dist)),
}

# ------------------------- Core helper (shared) -------------------------

def _soft_neighbor_counts_per_A(
    cell: torch.Tensor,
    frac: torch.Tensor,
    types: torch.Tensor,
    type_A: int | str,
    type_B: int | str,
    *,
    kernel: str | Callable = "sigmoid",
    sigma: float = 1.0,
    r_cut: float | None = None,
    alpha: float = 8.0,
    remove_self: bool = True,
    images_27: bool = True,
    **kwargs,
) -> torch.Tensor:
    """
    Compute per-A soft neighbor counts C_i = sum_{j in A_B} g(d_ij).
    Periodic boundary conditions handled via 27-image expansion by default.

    Args:
      cell   : (3,3) Cartesian cell tensor (Å).
      frac   : (N,3) fractional coordinates (0..1).
      types  : (N,) atomic numbers (int).
      type_A : A species (int Z or symbol).
      type_B : B species (int Z or symbol).
      kernel : "sigmoid", "gaussian", or callable g(dist)->[0,1].
      sigma, r_cut, alpha : kernel parameters.
      remove_self : if True and A==B, subtract g(0) once per A atom.
      images_27  : if True use 27 images; otherwise assumes minimum-image not needed.
    Returns:
      counts : (n_A,) tensor of soft counts for each A atom.
    """

    ker = kernel if hasattr(kernel, "__call__") else _STANDARD_KERNELS[kernel]

    ZA = _Z[type_A] if isinstance(type_A, str) else int(type_A)
    ZB = _Z[type_B] if isinstance(type_B, str) else int(type_B)

    device = frac.device
    types = types.to(torch.int64)

    mask_A = (types == ZA)
    mask_B = (types == ZB)
    idx_A = mask_A.nonzero(as_tuple=True)[0]
    idx_B = mask_B.nonzero(as_tuple=True)[0]

    # Handle degenerate cases with a 0.0 tensor that preserves graph
    if idx_A.numel() == 0 or idx_B.numel() == 0:
        return cell.sum()[None] * 0.0

    if r_cut is None:
        # A relaxed, chemistry-informed default
        r_cut = float(covalent_radii[ZA] + covalent_radii[ZB] + 0.5)

    # Build image shifts
    if images_27:
        rng = torch.arange(-1, 2, device=device)
        shifts = torch.stack(torch.meshgrid(rng, rng, rng, indexing="ij"), dim=-1).reshape(-1, 3)  # (27,3)
    else:
        shifts = torch.zeros((1, 3), device=device)

    frac_A = frac[idx_A]                   # (n_A,3)
    frac_B = frac[idx_B]                   # (n_B,3)
    n_A = frac_A.shape[0]

    # Expand B over images
    frac_B_images = (frac_B.unsqueeze(1) + shifts.unsqueeze(0)).reshape(-1, 3)  # (n_B*27, 3)

    # Pairwise A -> (B, images) distances under PBC
    dfrac = frac_A.unsqueeze(1) - frac_B_images.unsqueeze(0)   # (n_A, n_B*27, 3)
    dcart = torch.matmul(dfrac, cell)                           # (n_A, n_B*27, 3)
    dist  = dcart.norm(dim=-1)                                  # (n_A, n_B*27)

    G = ker(dist, sigma=sigma, r_cut=r_cut, alpha=alpha, **kwargs)  # (n_A, n_B*27)
    counts = G.sum(dim=1)                                       # (n_A,)

    # Remove self contributions when A==B: subtract g(0) once per A atom.
    if remove_self and ZA == ZB:
        g0 = ker(torch.zeros(1, device=device), sigma=sigma, r_cut=r_cut, alpha=alpha, **kwargs)[0]
        counts = counts - g0

    return counts

# ------------------------- Public metrics -------------------------

def compute_average_coordination(
    cell: torch.Tensor,
    frac: torch.Tensor,
    types: torch.Tensor,
    type_A: int | str,
    type_B: int | str,
    kernel: str | Callable = "sigmoid",
    sigma: float = 1.0,
    r_cut: float | None = None,
    alpha: float = 8.0,
    **kwargs,
) -> torch.Tensor:
    r"""
    Old metric (kept for backward compatibility):
        E(A,B) = (1/|A|) * sum_{i in A_A} sum_{j in A_B} g(d_{ij})

    Returns a scalar torch.Tensor with gradients to cell and positions.
    """
    C = _soft_neighbor_counts_per_A(
        cell, frac, types, type_A, type_B,
        kernel=kernel, sigma=sigma, r_cut=r_cut, alpha=alpha, **kwargs
    )  # (n_A,)
    if C.numel() == 0:
        return cell.sum() * 0.0
    return C.mean()

def compute_target_share(
    cell: torch.Tensor,
    frac: torch.Tensor,
    types: torch.Tensor,
    type_A: int | str,
    type_B: int | str,
    *,
    target: float,
    tau: float = 0.5,
    kernel: str | Callable = "sigmoid",
    sigma: float = 1.0,
    r_cut: float | None = None,
    alpha: float = 8.0,
    **kwargs,
) -> torch.Tensor:
    r"""
    New, "resolved" metric:
        (1/|A|) * sum_i  h( sum_j g(d_{ij}) - target ),
    with h(x) = exp( - (x/tau)^2 ).

    Returns a scalar torch.Tensor in (0,1] ~ fraction of A atoms
    whose soft B-neighbor count is near `target`.
    """
    C = _soft_neighbor_counts_per_A(
        cell, frac, types, type_A, type_B,
        kernel=kernel, sigma=sigma, r_cut=r_cut, alpha=alpha, **kwargs
    )  # (n_A,)
    if C.numel() == 0:
        return cell.sum() * 0.0
    H = torch.exp(-((C - float(target)) / float(tau)).pow(2))
    return H.mean()

def compute_multiple_target_share(
    cell: torch.Tensor,
    frac: torch.Tensor,
    types: torch.Tensor,
    type_A: int | str,
    type_B: int | str,
    targets: Iterable[float],
    *,
    tau: float = 0.5,
    kernel: str | Callable = "sigmoid",
    sigma: float = 1.0,
    r_cut: float | None = None,
    alpha: float = 8.0,
    **kwargs,
) -> torch.Tensor:
    """
    Vectorized "histogram-like" shares for multiple targets.
    Returns a (K,) tensor with shares[k] ~ fraction near targets[k].
    """
    C = _soft_neighbor_counts_per_A(
        cell, frac, types, type_A, type_B,
        kernel=kernel, sigma=sigma, r_cut=r_cut, alpha=alpha, **kwargs
    )  # (n_A,)
    if C.numel() == 0:
        return cell.sum()[None] * 0.0
    T = torch.as_tensor(list(targets), dtype=C.dtype, device=C.device).view(1, -1)  # (1,K)
    H = torch.exp(-((C.view(-1, 1) - T) / float(tau)).pow(2))                        # (n_A,K)
    return H.mean(dim=0)                                                             # (K,)

# ------------------------- ASE convenience wrappers -------------------------

def _to_torch_from_atoms(atoms):
    cell = torch.tensor(np.asarray(atoms.cell), dtype=torch.float64, device=device)
    frac = torch.tensor(atoms.get_scaled_positions(), dtype=torch.float64, device=device)
    types = torch.tensor(atoms.get_atomic_numbers(), dtype=torch.int64, device=device)
    return cell, frac, types

def compute_average_coordination_atoms(atoms, type_A, type_B, kernel="sigmoid", sigma=1.0, r_cut=None, alpha=8.0):
    cell, frac, types = _to_torch_from_atoms(atoms)
    return compute_average_coordination(cell, frac, types, type_A, type_B,
                                        kernel=kernel, sigma=sigma, r_cut=r_cut, alpha=alpha)

def compute_target_share_atoms(atoms, type_A, type_B, *, target, tau=0.5,
                               kernel="sigmoid", sigma=1.0, r_cut=None, alpha=8.0):
    cell, frac, types = _to_torch_from_atoms(atoms)
    return compute_target_share(cell, frac, types, type_A, type_B,
                                target=target, tau=tau, kernel=kernel,
                                sigma=sigma, r_cut=r_cut, alpha=alpha)

def compute_multi_target_share_atoms(atoms, type_A, type_B, targets, *, tau=0.5,
                                           kernel="sigmoid", sigma=1.0, r_cut=None, alpha=8.0):
    cell, frac, types = _to_torch_from_atoms(atoms)
    return compute_multiple_target_share(cell, frac, types, type_A, type_B, targets,
                                         tau=tau, kernel=kernel, sigma=sigma, r_cut=r_cut, alpha=alpha)