import copy
from pathlib import Path

import numpy as np
import torch
from ....external_paths import add_sys_path, import_from_path_validator, resolve_external_path

_MATTERGEN_IMPORT_ERROR = None
mgen_path = None

ChemGraph = None
LOSS_REGISTRY = {}
_soft_neighbor_counts_per_A_single = None
_clear_globals_impl = None
compute_mean_coordination = None
compute_target_share = None
volume = None
volume_pa = None
mattergen_chemgraph_fn_collection = {}
mattergen_cell_frac_types_fn_collection = {}


def _refresh_function_collections() -> None:
    global mattergen_chemgraph_fn_collection, mattergen_cell_frac_types_fn_collection
    mattergen_chemgraph_fn_collection = {"volume": volume, "volume_pa": volume_pa}
    mattergen_cell_frac_types_fn_collection = {
        "compute_mean_coordination": compute_mean_coordination,
        "compute_target_share": compute_target_share,
        "_soft_neighbor_counts_per_A_single": _soft_neighbor_counts_per_A_single,
    }


def _configure_mattergen_path(*, prompt: bool) -> Path | None:
    global mgen_path
    if mgen_path is not None:
        return Path(mgen_path)

    resolved_path = resolve_external_path(
        name="MatterGen source tree",
        config_key="mattergen_python_path",
        env_var="MATTERGEN_PYTHON_PATH",
        validator=import_from_path_validator(
            ("mattergen.common.data.chemgraph", "mattergen.diffusion.diffusion_loss")
        ),
        prompt=prompt,
        required=prompt,
        prompt_text="Enter path to MatterGen source tree: ",
    )
    if resolved_path is not None:
        add_sys_path(resolved_path, prepend=False)
        mgen_path = resolved_path
    return resolved_path


def _import_mattergen() -> None:
    global ChemGraph, LOSS_REGISTRY, _soft_neighbor_counts_per_A_single, _clear_globals_impl
    global compute_mean_coordination, compute_target_share, volume, volume_pa
    global _MATTERGEN_IMPORT_ERROR

    try:
        from mattergen.common.data.chemgraph import ChemGraph as _ChemGraph
        from mattergen.diffusion.diffusion_loss import (
            LOSS_REGISTRY as _LOSS_REGISTRY,
            _soft_neighbor_counts_per_A_single as _soft_neighbor_counts_per_A_single_impl,
            clear_globals as _clear_globals_impl_imported,
            compute_mean_coordination as _compute_mean_coordination,
            compute_target_share as _compute_target_share,
            volume as _volume,
            volume_pa as _volume_pa,
        )
    except ModuleNotFoundError as exc:
        _MATTERGEN_IMPORT_ERROR = exc
        return

    ChemGraph = _ChemGraph
    LOSS_REGISTRY = _LOSS_REGISTRY
    _soft_neighbor_counts_per_A_single = _soft_neighbor_counts_per_A_single_impl
    _clear_globals_impl = _clear_globals_impl_imported
    compute_mean_coordination = _compute_mean_coordination
    compute_target_share = _compute_target_share
    volume = _volume
    volume_pa = _volume_pa
    _MATTERGEN_IMPORT_ERROR = None
    _refresh_function_collections()


_configure_mattergen_path(prompt=False)
_import_mattergen()


def _require_mattergen():
    if ChemGraph is None:
        _configure_mattergen_path(prompt=True)
        _import_mattergen()
    if ChemGraph is None:
        if mgen_path is None:
            raise ModuleNotFoundError(
                "mattergen is not configured. Set MATTERGEN_PYTHON_PATH, configure external_paths.json, "
                "or enter the path when prompted."
            ) from _MATTERGEN_IMPORT_ERROR
        raise ModuleNotFoundError(
            f"mattergen import failed from path '{mgen_path}'."
        ) from _MATTERGEN_IMPORT_ERROR


def _device_from_force_gpu(force_gpu: int) -> torch.device:
    if not isinstance(force_gpu, int):
        raise TypeError(f"force_gpu must be int, got {type(force_gpu).__name__}")
    if torch.cuda.is_available():
        cuda_count = torch.cuda.device_count()
        if force_gpu < 0 or force_gpu >= cuda_count:
            raise ValueError(f"force_gpu={force_gpu} is out of range for {cuda_count} visible CUDA device(s)")
        return torch.device(f"cuda:{force_gpu}")
    return torch.device("cpu")


def structure_to_tensors(struct, force_gpu: int):
    _require_mattergen()
    device = _device_from_force_gpu(force_gpu)

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


def structure_to_single_chemgraph(struct, force_gpu: int):
    _require_mattergen()
    cell, pos, atomic_numbers = structure_to_tensors(struct, force_gpu=force_gpu)
    return ChemGraph(
        cell=cell,
        atomic_numbers=atomic_numbers,
        pos=pos,
        num_atoms=torch.tensor([len(atomic_numbers)], device=cell.device),
    )


def entry2tensors(entry, force_gpu: int):
    return structure_to_tensors(entry.structure, force_gpu=force_gpu)


def entry2chemgraph(entry, force_gpu: int):
    return structure_to_single_chemgraph(entry.structure, force_gpu=force_gpu)


def clear_globals():
    _require_mattergen()
    _clear_globals_impl()


def get_target_value_fn(fn_name, force_gpu: int = 0, **params):
    _require_mattergen()
    fn = lambda x: None
    if fn_name in mattergen_chemgraph_fn_collection:
        def fn(entry):
            x = entry2chemgraph(entry, force_gpu=force_gpu)
            return mattergen_chemgraph_fn_collection[fn_name](x, t=None, **params).cpu().detach().numpy()[0]
    elif fn_name in mattergen_cell_frac_types_fn_collection:
        def fn(entry):
            cell, frac, types = structure_to_tensors(entry.structure, force_gpu=force_gpu)
            return mattergen_cell_frac_types_fn_collection[fn_name](
                cell,
                frac,
                types,
                num_atoms=torch.tensor([len(types)], device=cell.device),
                **params,
            ).cpu().detach().numpy().item()
    return fn


def get_loss_fn(fn_name, force_gpu: int = 0, **params):
    _require_mattergen()
    assert fn_name in LOSS_REGISTRY, f"Loss function: {fn_name} not implemented"

    def fn(entry):
        assert "target" in params, "{target: target_dict} is required for loss calculation"
        clear_globals()
        x = entry2chemgraph(entry, force_gpu=force_gpu)
        call_params = copy.deepcopy(params)
        return LOSS_REGISTRY[fn_name](x, t=None, **call_params).cpu().detach().numpy().item()

    return fn
