import copy
import os
import socket
import sys
from pathlib import Path

import numpy as np
import torch

MATTERGEN_PYTHON_PATHS = {
    "nina": "/home/vsbat/work/mattergen/mattergenbis_vb",
    "taurus": "/home/vsbat/my_git_projects/mattergenbis_vb",
    "serpens": "/home/vsbat/git_packages/mattergenbis_vb",
}

_MATTERGEN_IMPORT_ERROR = None

host = socket.gethostname()
mgen_path = os.environ.get("MATTERGEN_PYTHON_PATH") or MATTERGEN_PYTHON_PATHS.get(host)
if mgen_path:
    sys.path.append(Path(mgen_path).as_posix())

try:
    from mattergen.common.data.chemgraph import ChemGraph
    from mattergen.diffusion.diffusion_loss import (
        LOSS_REGISTRY,
        _soft_neighbor_counts_per_A_single,
        clear_globals as _clear_globals_impl,
        compute_mean_coordination,
        compute_target_share,
        volume,
        volume_pa,
    )
except ModuleNotFoundError as exc:
    ChemGraph = None
    LOSS_REGISTRY = {}
    _soft_neighbor_counts_per_A_single = None
    _clear_globals_impl = None
    compute_mean_coordination = None
    compute_target_share = None
    volume = None
    volume_pa = None
    _MATTERGEN_IMPORT_ERROR = exc

mattergen_chemgraph_fn_collection = {"volume": volume, "volume_pa": volume_pa}
mattergen_cell_frac_types_fn_collection = {
    "compute_mean_coordination": compute_mean_coordination,
    "compute_target_share": compute_target_share,
    "_soft_neighbor_counts_per_A_single": _soft_neighbor_counts_per_A_single,
}


def _require_mattergen():
    if ChemGraph is None:
        if mgen_path is None:
            raise ModuleNotFoundError(
                "mattergen is not configured. Set MATTERGEN_PYTHON_PATH or add hostname to MATTERGEN_PYTHON_PATHS."
            ) from _MATTERGEN_IMPORT_ERROR
        raise ModuleNotFoundError(
            f"mattergen import failed from path '{mgen_path}'."
        ) from _MATTERGEN_IMPORT_ERROR


def structure_to_tensors(struct):
    _require_mattergen()
    device = "cuda" if torch.cuda.is_available() else "cpu"

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
    _require_mattergen()
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


def clear_globals():
    _require_mattergen()
    _clear_globals_impl()


def get_target_value_fn(fn_name, **params):
    _require_mattergen()
    fn = lambda x: None
    if fn_name in mattergen_chemgraph_fn_collection:
        def fn(entry):
            x = entry2chemgraph(entry)
            return mattergen_chemgraph_fn_collection[fn_name](x, t=None, **params).cpu().detach().numpy()[0]
    elif fn_name in mattergen_cell_frac_types_fn_collection:
        def fn(entry):
            cell, frac, types = structure_to_tensors(entry.structure)
            return mattergen_cell_frac_types_fn_collection[fn_name](
                cell,
                frac,
                types,
                num_atoms=torch.tensor([len(types)], device=cell.device),
                **params,
            ).cpu().detach().numpy().item()
    return fn


def get_loss_fn(fn_name, **params):
    _require_mattergen()
    assert fn_name in LOSS_REGISTRY, f"Loss function: {fn_name} not implemented"

    def fn(entry):
        assert "target" in params, "{target: target_dict} is required for loss calculation"
        clear_globals()
        x = entry2chemgraph(entry)
        call_params = copy.deepcopy(params)
        return LOSS_REGISTRY[fn_name](x, t=None, **call_params).cpu().detach().numpy().item()

    return fn

