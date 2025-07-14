import numpy as np
from ase.io import read
import torch
from ase.data import covalent_radii, chemical_symbols, atomic_numbers
from itertools import product

# SIGMOID = lambda dist, cutoff, delta: 1.0 / (1.0 + np.exp((dist - cutoff) / delta))
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def _compute_species_pair_single(
        cell: torch.Tensor,
        frac: torch.Tensor,
        types: torch.Tensor,
        type_A: int,
        type_B: int,
        kernel: str = "gaussian",
        sigma: float = 1.0,
        r_cut: float | None = None,
        alpha: float = 8.0
) -> torch.Tensor:
    """
    Compute a differentiable species‐pair value f[A,B], where
        f[A,B] = (1 / |A_A|) * sum_{i in A_A} sum_{j in A_B} g(d_ij),
    under PBC.  g(d) is either:
      - Gaussian:      exp[−(d/sigma)^2]
      - Soft‐cutoff:   sigmoid(alpha * (r_cut − d))

    Args:
      cell      (3×3) Tensor: rows are lattice vectors.
      frac      (N×3) Tensor: fractional coords.
      types     list of N ints: atomic numbers.
      type_A    int: atomic number of species A (the one we want to compute the environment of).
      type_B    int: atomic number of species B (the one we are looking in the environment of species A).
      kernel    "gaussian" or "sigmoid"
      sigma     width for Gaussian (Å)
      r_cut     cutoff for sigmoid (Å)
      alpha     sharpness for sigmoid

    Returns:
      f_AB (scalar Tensor), with gradients flowing to cell and positions.
    """
    type_A = atomic_numbers[type_A] if isinstance(type_A, str) else type_A
    type_B = atomic_numbers[type_B] if isinstance(type_B, str) else type_B
    device = frac.device
    mask_A = (types == type_A)
    mask_B = (types == type_B)
    idx_A = mask_A.nonzero(as_tuple=True)[0]
    idx_B = mask_B.nonzero(as_tuple=True)[0]

    if idx_A.numel() == 0 or idx_B.numel() == 0:
        return cell.sum() * 0.0  # No A or B atoms, return 0

    # Prepare cutoffs if needed
    if r_cut is None:
        r_cut = covalent_radii[type_A] + covalent_radii[type_B] + 1.# (N, N)

    # PBC images
    shifts = torch.stack(torch.meshgrid(
        torch.arange(-1, 2, device=device),
        torch.arange(-1, 2, device=device),
        torch.arange(-1, 2, device=device),
        indexing='ij'
    ), dim=-1).reshape(-1, 3)  # (27, 3)

    # Get frac for A and B
    frac_A = frac[idx_A]  # (n_A, 3)
    frac_B = frac[idx_B]  # (n_B, 3)

    # Expand B atoms to all images
    frac_B_images = frac_B.unsqueeze(1) + shifts.unsqueeze(0)  # (n_B, 27, 3)
    frac_B_images = frac_B_images.reshape(-1, 3)  # (n_B*27, 3)

    # Compute all distances from each A to all B images
    d = frac_A.unsqueeze(1) - frac_B_images.unsqueeze(0)  # (n_A, n_B*27, 3)
    dc = torch.matmul(d, cell)  # (n_A, n_B*27, 3)
    dist = dc.norm(dim=-1)  # (n_A, n_B*27)

    # Kernel
    if kernel == "gaussian":
        G = torch.exp(- (dist / sigma).pow(2))
    elif kernel == "sigmoid":
        G = torch.sigmoid(alpha * (r_cut - dist))
    else:
        raise ValueError("kernel must be 'gaussian' or 'sigmoid'")

    # Sum over all pairs
    n_A = mask_A.sum()
    if n_A < 1e-8:
        # to keep the graph differentiable, return 0
        return cell.sum() * 0.0
    f_AB = G.sum() / n_A - int(type_A == type_B)  # Subtract self-interaction
    return f_AB

def compute_species_pair(atoms, type_A, type_B, kernel="gaussian", sigma=1.0, r_cut=None, alpha=8.0):
    """
    Compute the species pair value for a given structure.

    Args:
        atoms: ASE Atoms object.
        type_A: Atomic number of species A.
        type_B: Atomic number of species B.
        kernel: Kernel type ('gaussian' or 'sigmoid').
        sigma: Width for Gaussian kernel.
        r_cut: Cutoff for sigmoid kernel.
        alpha: Sharpness for sigmoid kernel.

    Returns:
        f_AB (scalar Tensor).
    """

    cell = torch.tensor(np.array(atoms.cell), dtype=torch.float64, device=device)
    frac = torch.tensor(atoms.get_scaled_positions(), dtype=torch.float64, device=device)
    atomic_numbers = torch.tensor(atoms.get_atomic_numbers(), dtype=torch.int64, device=device)

    return _compute_species_pair_single(cell, frac, atomic_numbers, type_A, type_B, kernel, sigma, r_cut, alpha)

# Usage example:
# counts = count_neighbors('POSCAR', delta=0.15)
# print(counts)


INTER_ATOMIC_CUTOFF = {}

if __name__ == '__main__':
    atoms = read('/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/STRUCTURE_GENERATION/Guidance/Mo5SiB2_POSCAR')
    num_atoms = torch.tensor([len(atoms)], dtype=torch.int64, device=device)
    print(compute_species_pair(atoms, 42, 42, kernel="sigmoid", alpha=4800.0))

