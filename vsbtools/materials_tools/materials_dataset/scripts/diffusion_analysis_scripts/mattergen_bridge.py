import socket
import sys
from pathlib import Path
import copy
import torch
import numpy as np
from pymatgen.core import Structure
host = socket.gethostname()
MATTERGEN_PYTHON_PATHS = {'nina': "/home/vsbat/work/mattergen/mattergenbis_vb",
                          "taurus": "/home/vsbat/my_git_projects/mattergenbis_vb",
                          "serpens": "/home/vsbat/git_packages/mattergenbis_vb"}

if host not in MATTERGEN_PYTHON_PATHS:
    mgen_path = Path(input("Enter full path containing mattergen package"))
else:
    mgen_path = Path(MATTERGEN_PYTHON_PATHS[host])

sys.path.append(mgen_path.as_posix())
from mattergen.common.data.chemgraph import ChemGraph
from mattergen.diffusion.diffusion_loss import (volume, volume_pa, compute_mean_coordination, compute_target_share,
                                                _soft_neighbor_counts_per_A_single, LOSS_REGISTRY, clear_globals)

mattergen_chemgraph_fn_collection = {"volume": volume, "volume_pa": volume_pa}

mattergen_cell_frac_types_fn_collection = {"compute_mean_coordination": compute_mean_coordination,
                                           "compute_target_share": compute_target_share,
                                           "_soft_neighbor_counts_per_A_single": _soft_neighbor_counts_per_A_single}

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

def get_target_value_fn(fn_name, **params):
    fn = lambda x: None
    if fn_name in mattergen_chemgraph_fn_collection:
        def fn(entry):
            x = entry2chemgraph(entry)
            return mattergen_chemgraph_fn_collection[fn_name](x, t=None, **params).cpu().detach().numpy()[0]  # x is a batch normally, so to have a result of entry we need [0]
    elif fn_name in mattergen_cell_frac_types_fn_collection:
        def fn(entry):
            cell, frac, types = structure_to_tensors(entry.structure)
            return mattergen_cell_frac_types_fn_collection[fn_name](cell, frac, types,
                                                                    num_atoms=torch.tensor([len(types)],
                                                                                           device = cell.device),
                                                                    **params).cpu().detach().numpy().item()
    return fn

def get_loss_fn(fn_name, **params):
    assert fn_name in LOSS_REGISTRY, f"Loss function: {fn_name} not implemented"
    def fn(entry):
        assert 'target' in params, "{target: target_dict} is required for loss calculation"
        clear_globals()
        x = entry2chemgraph(entry)
        call_params = copy.deepcopy(params)
        return LOSS_REGISTRY[fn_name](x, t=None, **call_params).cpu().detach().numpy().item()
    return fn
