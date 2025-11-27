import sys
import torch
import numpy as np
from pymatgen.core import Structure
sys.path.append("/home/vsbat/my_git_projects/mattergenbis_vb")
from mattergen.common.data.chemgraph import ChemGraph
from mattergen.diffusion.diffusion_loss import volume, volume_pa, compute_mean_coordination, compute_target_share

mattergen_fn_collection = {"volume": volume, "volume_pa": volume_pa,
                        compute_mean_coordination: "compute_mean_coordination",
                        compute_target_share: "compute_target_share"}

def structure_to_tensors(struct):
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # 3×3 lattice matrix (a, b, c as rows)
    cell = torch.tensor(
        np.array(struct.lattice.matrix),
        dtype=torch.float32,
        device=device,
    )

    # fractional coordinates
    pos = torch.tensor(
        np.array(struct.frac_coords),
        dtype=torch.float32,
        requires_grad=True,
        device=device,
    )

    # atomic numbers
    atomic_numbers = torch.tensor(
        [site.specie.number for site in struct],
        dtype=torch.int64,
        device=device,
    )

    return cell, pos, atomic_numbers


def structure_to_single_chemgraph(struct):
    cell, pos, atomic_numbers = structure_to_tensors(struct)
    return ChemGraph(
        cell=cell,
        atomic_numbers=atomic_numbers,
        pos=pos,
        num_atoms=torch.tensor([len(atomic_numbers)], device=cell.device),
    )


def entry2tensors(entry):
    return structure_to_tensors(entry.structure)


def entry2chemgraph(entry):
    return structure_to_single_chemgraph(entry.structure)
