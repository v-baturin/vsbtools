# scenario_pipeline.py

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, MutableMapping, Optional, Tuple
from collections import defaultdict, deque

from ..crystal_dataset import CrystalDataset
from ..io.structures_dataset_io import StructureDatasetIO
from .symmetry_tools import SymmetryToolkit
from ..scripts.poll_databases import poll_databases
from .similarity_tools import SimilarityTools
from ..io.uspex_bridge import USPEXBridge
from ..energy_estimation import nn_estimator, mattersim_bridge
from ..analysis.phase_diagram_tools import PhaseDiagramTools

from .pipeline_legacy import LEGACY_INDEX_TO_NAME, LEGACY_NAME_TO_INDEX

# --- global constants / toolkit registry ---------------------------

TOOLKIT_DICT: Dict[str, Any] = {
    "structure_parser": StructureDatasetIO,
    "symmetry":         SymmetryToolkit,
    "similarity":       SimilarityTools,
    "uspex":            USPEXBridge,
    "phase_diag":       PhaseDiagramTools,
    "estimator":        nn_estimator.NNEstimator,
}

KNOWN_MODELS = {"mattersim": mattersim_bridge}
MAX_EHULL_PA = 0.1


# ==================================================================
# 1. Context: toolkits, global settings, stage outputs
# ==================================================================

@dataclass
class Context:
    """Execution context for the pipeline."""
    toolkit_options: MutableMapping[str, Dict[str, Any]] = field(
        default_factory=lambda: defaultdict(dict)
    )
    toolkits: Dict[str, Any] = field(default_factory=dict)
    globals: Dict[str, Any] = field(default_factory=dict)
    outputs: Dict[str, CrystalDataset] = field(default_factory=dict)

    def get_tool(self, name: str):
        """Lazy initialization of a toolkit by name from TOOLKIT_DICT."""
        if name not in TOOLKIT_DICT:
            raise KeyError(f"Unknown toolkit '{name}'")

        if name not in self.toolkits:
            # special handling for some toolkits
            # Auto-fill USPEX elements from ctx.globals if missing
            if name == "uspex":
                opts = self.toolkit_options["uspex"]
                if "elements" not in opts:
                    elems = self.globals.get("elements")
                    if elems is None:
                        raise RuntimeError(
                            "USPEXBridge requires 'elements'; "
                            "set ctx.globals['elements'] before using 'uspex'"
                        )
                    opts["elements"] = elems

            if name == "similarity":
                ub = self.get_tool("uspex")
                self.toolkit_options["similarity"].setdefault("dist_fn", ub.fp_dist)
            elif name == "estimator":
                for label, module in KNOWN_MODELS.items():
                    TOOLKIT_DICT["estimator"].register_model(label, module)

            options = self.toolkit_options.get(name, {})
            self.toolkits[name] = TOOLKIT_DICT[name](**options)

        return self.toolkits[name]


# ==================================================================
# 2. Operations and registry
# ==================================================================

OperationFn = Callable[[Context, Dict[str, CrystalDataset], Dict[str, Any]],
                       CrystalDataset]

OP_REGISTRY: Dict[str, OperationFn] = {}


def op(name: str) -> Callable[[OperationFn], OperationFn]:
    """Register an operation in the global registry."""
    def _wrap(fn: OperationFn) -> OperationFn:
        if name in OP_REGISTRY:
            raise ValueError(f"Operation '{name}' already registered")
        OP_REGISTRY[name] = fn
        return fn
    return _wrap


# ==================================================================
# 3. Scenario model and topological ordering
# ==================================================================

@dataclass
class StageSpec:
    """Specification of a single stage from the scenario."""
    name: str
    op: str
    needs: List[str] = field(default_factory=list)
    params: Dict[str, Any] = field(default_factory=dict)
    description: Optional[str] = None
    legacy_index: Optional[int] = None


@dataclass
class Scenario:
    version: int
    globals: Dict[str, Any] = field(default_factory=dict)
    stages: Dict[str, StageSpec] = field(default_factory=dict)

    # Legacy mapping: numeric index <-> stage name
    legacy_index_to_name: Dict[int, str] = field(init=False, default_factory=dict)
    legacy_name_to_index: Dict[str, int] = field(init=False, default_factory=dict)

    @staticmethod
    def from_mapping(cfg: Dict[str, Any]) -> "Scenario":
        if "version" not in cfg:
            raise ValueError("Scenario must have a 'version' field")

        g = cfg.get("globals", {}) or {}
        raw_stages = cfg.get("stages", {}) or {}
        stages: Dict[str, StageSpec] = {}

        for name, sc in raw_stages.items():
            stages[name] = StageSpec(
                name=name,
                op=sc["op"],
                needs=list(sc.get("needs", []) or []),
                params=dict(sc.get("params", {}) or {}),
                description=sc.get("description"),
                # legacy_index=int(sc["legacy_index"]) if "legacy_index" in sc else None,
            )

        scenario = Scenario(version=int(cfg["version"]), globals=g, stages=stages)
        scenario._build_legacy_mappings()
        return scenario

    def _build_legacy_mappings(self) -> None:
        """Initialize legacy index/name mappings from legacy module and scenario."""
        # Base mapping from legacy module
        self.legacy_index_to_name.update(LEGACY_INDEX_TO_NAME)
        self.legacy_name_to_index.update(LEGACY_NAME_TO_INDEX)

        # Optional overrides from scenario (if you decide to use legacy_index in StageSpec)
        for name, spec in self.stages.items():
            legacy_idx = getattr(spec, "legacy_index", None)
            if legacy_idx is None:
                continue
            idx = int(legacy_idx)
            if idx in self.legacy_index_to_name and self.legacy_index_to_name[idx] != name:
                raise ValueError(
                    f"Legacy index {idx} is assigned to both "
                    f"'{self.legacy_index_to_name[idx]}' and '{name}'"
                )
            self.legacy_index_to_name[idx] = name
            self.legacy_name_to_index[name] = idx


    def resolve_stage_name_from_metadata(self, ds: CrystalDataset) -> str:
        """Resolve canonical stage name from dataset metadata.

        Supports:
        - legacy int in metadata["pipeline_stage"] (0, 1, 2, ...)
        - string stage name in metadata["pipeline_stage"]
        - stringified int in metadata["pipeline_stage"]
        """
        meta: Dict[str, Any] = getattr(ds, "metadata", None) or {}
        raw = meta.get("pipeline_stage")

        # Case 1: direct string name
        if isinstance(raw, str):
            if raw in self.stages:
                return raw
            # maybe it's a stringified integer; try legacy mapping
            try:
                idx = int(raw)
            except ValueError:
                idx = None
            if idx is not None and idx in self.legacy_index_to_name:
                return self.legacy_index_to_name[idx]
            raise KeyError(f"Cannot resolve stage name from string '{raw}'")

        # Case 2: direct int (legacy index)
        if isinstance(raw, int):
            if raw in self.legacy_index_to_name:
                return self.legacy_index_to_name[raw]
            raise KeyError(f"Unknown legacy stage index: {raw}")

        raise KeyError(f"Cannot resolve stage name from metadata: {raw!r}")


def load_scenario_file(path: Path) -> Scenario:
    """Load scenario from YAML or JSON."""
    text = path.read_text(encoding="utf-8")
    suffix = path.suffix.lower()

    if suffix in {".yml", ".yaml"}:
        import yaml  # PyYAML
        data = yaml.safe_load(text)
    elif suffix == ".json":
        data = json.loads(text)
    else:
        raise ValueError(f"Unsupported scenario format: '{suffix}'")

    if not isinstance(data, dict):
        raise TypeError("Scenario root must be a mapping (dict)")
    return Scenario.from_mapping(data)


def topo_order(stages: Dict[str, StageSpec]) -> List[str]:
    """Topological sort of stage DAG.

    Returns a list of stage names such that each stage appears after all its
    dependencies. Raises RuntimeError if a cycle is detected.
    """
    indeg: Dict[str, int] = {k: 0 for k in stages}
    graph: Dict[str, List[str]] = {k: [] for k in stages}

    for name, spec in stages.items():
        for parent in spec.needs:
            if parent not in stages:
                raise KeyError(f"Unknown parent '{parent}' for stage '{name}'")
            indeg[name] += 1
            graph[parent].append(name)

    q = deque([k for k, d in indeg.items() if d == 0])
    order: List[str] = []

    while q:
        u = q.popleft()
        order.append(u)
        for v in graph[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                q.append(v)

    if len(order) != len(stages):
        raise RuntimeError("Cycle detected in scenario stages")

    return order


# ==================================================================
# 4. Scenario executor
# ==================================================================

class ScenarioPipeline:
    def __init__(self, scenario: Scenario):
        self.scenario = scenario
        self.ctx = Context()

        # clear USPEX cache as in original PPPipeline.__post_init__
        USPEXBridge.uspex_entry_from_de.cache_clear()

        # merge globals.toolkit_options into context
        tk_opts = scenario.globals.get("toolkit_options", {}) or {}
        if not isinstance(self.ctx.toolkit_options, defaultdict):
            self.ctx.toolkit_options = defaultdict(dict, tk_opts)
        else:
            self.ctx.toolkit_options.update(tk_opts)

        # other globals
        self.ctx.globals = {
            k: v for k, v in scenario.globals.items()
            if k != "toolkit_options"
        }

    @classmethod
    def from_file(cls, path: Path) -> "ScenarioPipeline":
        return cls(load_scenario_file(path))

    def _stage_sequence(
        self,
        targets: Optional[Iterable[str]],
        skip: Optional[Iterable[str]],
    ) -> List[str]:
        order = topo_order(self.scenario.stages)
        all_names = set(order)

        if targets is None:
            # By default, run only stages that are not already in outputs
            targets_set = all_names - set(self.ctx.outputs.keys())
        else:
            targets_set = set(targets)
            unknown_targets = targets_set - all_names
            if unknown_targets:
                raise KeyError(f"Unknown target stages: {sorted(unknown_targets)}")

        skip_set = set(skip) if skip else set()
        return [name for name in order if name in targets_set and name not in skip_set]

    def run(
        self,
        targets: Optional[Iterable[str]] = None,
        skip: Optional[Iterable[str]] = None,
    ) -> Iterable[Tuple[str, CrystalDataset]]:
        """
        Execute the scenario.

        Yields pairs (stage_name, CrystalDataset).
        """
        for name in self._stage_sequence(targets, skip):
            spec = self.scenario.stages[name]
            inputs = {p: self.ctx.outputs[p] for p in spec.needs}

            if spec.op not in OP_REGISTRY:
                raise KeyError(f"Operation '{spec.op}' not registered")

            params = dict(spec.params)  # copy to avoid mutating original scenario
            ds = OP_REGISTRY[spec.op](self.ctx, inputs, params)

            # minimal provenance
            if ds.metadata is None:
                ds.metadata = {}
            else:
                ds.metadata = dict(ds.metadata)

            ds.metadata["pipeline_stage"] = name

            ds.parent_ids = [
                getattr(inp, "dataset_id", None)
                for inp in inputs.values()
                if hasattr(inp, "dataset_id")
            ]

            idx = self.scenario.legacy_name_to_index.get(name)
            if idx is not None:
                ds.metadata.setdefault("pipeline_stage_index", str(idx))

            self.ctx.outputs[name] = ds
            yield name, ds


# ==================================================================
# 5. Operations implementing the original PPPipeline stages
# ==================================================================

# --- parse_raw ----------------------------------------------------

@op("parse_raw")
def op_parse_raw(
    ctx: Context,
    inputs: Dict[str, CrystalDataset],
    params: Dict[str, Any],
) -> CrystalDataset:
    if inputs:
        raise AssertionError("parse_raw must have no parents")

    # elements can come from params or from globals
    elements = params.get("elements", ctx.globals.get("elements"))
    if elements is None:
        raise AssertionError("elements must be specified (in params or globals)")

    parser: StructureDatasetIO = ctx.get_tool("structure_parser")
    ds = parser.load_from_directory(elements=elements)

    if len(ds) == 0:
        raise AssertionError("parse_raw produced an empty dataset")

    if ds.metadata is None:
        ds.metadata = {}
    ds.metadata["message"] = f"loaded {len(ds)} structures"

    return ds


# --- symmetrize_raw ----------------------------------------------

@op("symmetrize_raw")
def op_symmetrize_raw(
    ctx: Context,
    inputs: Dict[str, CrystalDataset],
    params: Dict[str, Any],
) -> CrystalDataset:
    if len(inputs) != 1:
        raise AssertionError("symmetrize_raw expects exactly one parent")
    parent = next(iter(inputs.values()))

    symtool: SymmetryToolkit = ctx.get_tool("symmetry")
    return symtool.get_symmetrized_dataset(parent)


# --- poll_db ------------------------------------------------------

@op("poll_db")
def op_poll_db(
    ctx: Context,
    inputs: Dict[str, CrystalDataset],
    params: Dict[str, Any],
) -> CrystalDataset:
    if len(inputs) != 1:
        raise AssertionError("poll_db expects exactly one parent (parse_raw)")

    parse_ds = next(iter(inputs.values()))
    elements = getattr(parse_ds, "elements", None)
    if elements is None:
        raise AssertionError("parse_raw dataset must expose 'elements'")

    # default max_ehull as in the original code
    max_ehull = params.setdefault("max_ehull", MAX_EHULL_PA)

    # optional additional estimation
    estimate_energies = bool(params.pop("estimate_energies", False))

    # remaining parameters are passed directly to poll_databases
    db = poll_databases(elements, **params)

    if not estimate_energies:
        return db

    estimator = ctx.get_tool("estimator")
    estimated = estimator.estimate_dataset_energies(db)

    # parent_ids is empty as in the original code
    estimated.parent_ids = []

    msg0 = (db.metadata or {}).get("message", "").strip()
    msg1 = (estimated.metadata or {}).get("message", "").strip()
    estimated.metadata = dict(estimated.metadata or {})
    if msg0 and msg1:
        estimated.metadata["message"] = f"{msg0}. {msg1}"
    elif msg0:
        estimated.metadata["message"] = msg0
    elif msg1:
        estimated.metadata["message"] = msg1

    return estimated


# --- augment_raw_by_db -------------------------------------------

@op("augment_raw_by_db")
def op_augment_raw_by_db(
    ctx: Context,
    inputs: Dict[str, CrystalDataset],
    params: Dict[str, Any],
) -> CrystalDataset:
    if set(inputs.keys()) != {"symmetrize_raw", "poll_db"}:
        # allow arbitrary names but require exactly two parents
        if len(inputs) != 2:
            raise AssertionError("augment_raw_by_db expects symmetrize_raw and poll_db as parents")

    sym_ds = inputs.get("symmetrize_raw", list(inputs.values())[0])
    db_ds = inputs.get("poll_db", list(inputs.values())[1])

    similarity_tk: SimilarityTools = ctx.get_tool("similarity")
    augmentation = similarity_tk.get_unseen_in_ref(sym_ds, ref_ds=db_ds)

    print(augmentation.metadata.get("message", ""))

    # mark DB entries with reproduced=True
    reproduced = augmentation.metadata.get("reproduced", {})
    labeled_entries = []
    for e in db_ds:
        if str(e.id) in reproduced:
            new_metadata = {**(e.metadata or {}), "reproduced": True}
            labeled_entries.append(e.copy_with(metadata=new_metadata))
        else:
            labeled_entries.append(e)

    # reconstruct CrystalDataset parameters more carefully
    db_parameters = db_ds.__dict__.copy()
    meta = db_parameters.pop("metadata", {})
    base_path = db_parameters.pop("_base_path", None)
    if base_path is not None:
        db_parameters["base_path"] = base_path
    for k in list(db_parameters):
        if k.startswith("_"):
            db_parameters.pop(k, None)

    labeled_db = CrystalDataset(labeled_entries, **db_parameters)
    labeled_db.metadata = meta

    merged = labeled_db.merge(augmentation)
    merged.metadata = dict(merged.metadata or {})
    merged.metadata["message"] = "Merged DB data into generated set"
    merged.metadata["reproducibility"] = augmentation.metadata.get("reproducibility")

    merged.parent_ids = [
        getattr(sym_ds, "dataset_id", None),
        getattr(db_ds, "dataset_id", None),
    ]

    return merged


# --- estimate -----------------------------------------------------

@op("estimate")
def op_estimate(
    ctx: Context,
    inputs: Dict[str, CrystalDataset],
    params: Dict[str, Any],
) -> CrystalDataset:
    if len(inputs) != 1:
        raise AssertionError("estimate expects exactly one parent (augment_raw_by_db)")
    parent = next(iter(inputs.values()))

    estimator = ctx.get_tool("estimator")
    return estimator.estimate_dataset_energies(parent)


# --- filter_hull --------------------------------------------------

@op("filter_hull")
def op_filter_hull(
    ctx: Context,
    inputs: Dict[str, CrystalDataset],
    params: Dict[str, Any],
) -> CrystalDataset:
    if not inputs:
        raise AssertionError("filter_hull requires at least one parent")

    max_ehull = params.get("max_ehull", MAX_EHULL_PA)

    # dataset for phase-diagram construction
    base_stage = params.get("base_stage")
    if base_stage is not None:
        if base_stage not in inputs:
            raise KeyError(f"base_stage '{base_stage}' not among parents {list(inputs)}")
        base_ds = inputs[base_stage]
    else:
        base_ds = inputs.get("estimate", next(iter(inputs.values())))

    ctx.toolkit_options["phase_diag"] = {"dataset": base_ds}
    pd_tk: PhaseDiagramTools = ctx.get_tool("phase_diag")

    # dataset we actually filter
    select_from = inputs.get("estimate", base_ds)

    filtered = select_from.filter(
        lambda entry: pd_tk.height_above_hull_pa(entry) < max_ehull
    )
    filtered.metadata = dict(filtered.metadata or {})
    filtered.metadata["message"] = f"Selected entries with e_hull < {max_ehull:.3f}"

    return filtered


# --- deduplicate --------------------------------------------------

@op("deduplicate")
def op_deduplicate(
    ctx: Context,
    inputs: Dict[str, CrystalDataset],
    params: Dict[str, Any],
) -> CrystalDataset:
    if len(inputs) != 1:
        raise AssertionError("deduplicate expects exactly one parent (filter_hull)")
    parent = next(iter(inputs.values()))

    similarity_tk: SimilarityTools = ctx.get_tool("similarity")
    ds, _, _ = similarity_tk.deduplicate(
        parent,
        check_clusters_file=bool(params.get("check_clusters_file", True)),
        check_dist_matrix_file=bool(params.get("check_dist_matrix_file", True)),
    )

    ds.parent_ids = [getattr(parent, "dataset_id", None)]
    ds.metadata = dict(ds.metadata or {})
    ds.metadata["message"] = f"Deduplicated version of {getattr(parent, 'dataset_id', 'unknown')}"

    return ds


# --- postprocess_dft (placeholder) --------------------------------

@op("postprocess_dft")
def op_postprocess_dft(
    ctx: Context,
    inputs: Dict[str, CrystalDataset],
    params: Dict[str, Any],
) -> CrystalDataset:
    raise NotImplementedError("Implement DFT post-processing when needed")
