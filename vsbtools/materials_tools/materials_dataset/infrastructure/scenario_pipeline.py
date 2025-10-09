"""
Scenario-driven, DAG-based pipeline orchestrator for CrystalDataset postprocessing
------------------------------------------------------------------------------
Commercial-style, framework-free design that replaces a rigid, linear pipeline
with a declarative, scenario-driven Directed Acyclic Graph (DAG).

Key ideas
---------
- **Operations** are small, independent, side-effect-light functions.
- **Scenario** (JSON or YAML) declares nodes, their dependencies (parents),
  and parameters.
- **Registry** maps operation names to implementations, enabling plug-in ops.
- **Context** manages toolkits, options (dependency injection), and outputs.
- **Topological execution**: resolve dependencies (with cycle checks) and run.
- **Provenance**: attach minimal metadata (parents, stage name).
- **No external frameworks**: stdlib only; YAML is optional if PyYAML is present.

Scenario example (YAML)
-----------------------
# Save as scenario.yaml (or JSON with same structure)
version: 1
globals:
  source_path: ./data/raw
  elements: [Ca, B, H]
  root_source_name: "NA"
  toolkit_options:
    structure_parser:
      root: ./data/raw
      source_name: "NA"
    uspex:
      elements: [Ca, B, H]
    estimator: {}
    symmetry: {}
    similarity: {}
    phase_diag: {}

stages:
  parse_raw:
    op: parse_directory
    tool: structure_parser
    needs: []
    params:
      elements: [Ca, B, H]

  symmetrize_raw:
    op: symmetrize_dataset
    tool: symmetry
    needs: [parse_raw]

  poll_db:
    op: poll_databases
    tool: builtin
    needs: [parse_raw]
    params:
      max_ehull: 0.1
      estimate_energies: true
      elements_refdata_cache: ./cache/refdata

  augment_raw_by_db:
    op: augment_with_reference
    tool: similarity
    needs: [symmetrize_raw, poll_db]

  estimate:
    op: estimate_energies
    tool: estimator
    needs: [augment_raw_by_db]

  filter_hull:
    op: filter_by_hull
    tool: phase_diag
    needs: [estimate]
    params:
      max_ehull: 0.1
      # base_stage can override PD reference; default uses the sole parent
      # base_stage: estimate

  deduplicate:
    op: deduplicate
    tool: similarity
    needs: [filter_hull]
    params:
      check_clusters_file: true
      check_dist_matrix_file: true

Usage
-----
from pathlib import Path
from scenario_pipeline import ScenarioPipeline

pipe = ScenarioPipeline.from_file(Path("scenario.yaml"))
for name, ds in pipe.run(targets=None, skip=None):
    print("Produced:", name, len(ds))

Adding new operations
---------------------
@op("my_new_op")
def my_new_op(ctx, inputs, params):
    # inputs: dict[str, CrystalDataset], params: dict
    ...
    return output_dataset
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, MutableMapping, Optional, Set, Tuple
from collections import defaultdict, deque

# --- Project-specific imports (keep identical module paths to your project) ---
from ..crystal_dataset import CrystalDataset
from ..io.structures_dataset_io import StructureDatasetIO
from ..analysis.symmetry_tools import SymmetryToolkit
from ..scripts.poll_databases import poll_databases
from ..analysis.similarity_tools import SimilarityTools
from ..io.uspex_bridge import USPEXBridge
from ..energy_estimation import nn_estimator, mattersim_bridge
from ..analysis.phase_diagram_tools import PhaseDiagramTools

# ---- Toolkits registry (as in your current code) ----
TOOLKIT_DICT = {
    "structure_parser": StructureDatasetIO,
    "symmetry": SymmetryToolkit,
    "similarity": SimilarityTools,
    "uspex": USPEXBridge,
    "phase_diag": PhaseDiagramTools,
    "estimator": nn_estimator.NNEstimator,
}

KNOWN_MODELS = {"mattersim": mattersim_bridge}

MAX_EHULL_PA = 0.1

# ----------------------------------------------------------------------------
# Logging setup (adjust level/handlers in your app entrypoint if needed)
# ----------------------------------------------------------------------------
logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler()
    fmt = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(fmt)
    logger.addHandler(handler)
logger.setLevel(logging.INFO)

# ----------------------------------------------------------------------------
# Operation registry
# ----------------------------------------------------------------------------
OperationFn = Callable[["Context", Dict[str, CrystalDataset], Dict[str, Any]], CrystalDataset]

OP_REGISTRY: Dict[str, OperationFn] = {}

def op(name: str) -> Callable[[OperationFn], OperationFn]:
    """Decorator to register an operation under a stable name."""
    def _wrap(fn: OperationFn) -> OperationFn:
        if name in OP_REGISTRY:
            raise ValueError(f"Operation '{name}' already registered")
        OP_REGISTRY[name] = fn
        return fn
    return _wrap

# ----------------------------------------------------------------------------
# Context: toolkits, options, and produced artifacts
# ----------------------------------------------------------------------------
@dataclass
class Context:
    toolkit_options: MutableMapping[str, Dict[str, Any]] = field(default_factory=lambda: defaultdict(dict))
    toolkits: Dict[str, Any] = field(default_factory=dict)
    outputs: Dict[str, CrystalDataset] = field(default_factory=dict)

    def get_tool(self, toolkit_name: str):
        """Lazy-initialize toolkits and wire special cases."""
        if toolkit_name not in self.toolkits:
            if toolkit_name == "similarity":
                ub = self.get_tool("uspex")
                # Ensure similarity has an fp_dist function available from USPEX bridge
                self.toolkit_options["similarity"]["dist_fn"] = getattr(ub, "fp_dist", None)
            elif toolkit_name == "estimator":
                # Late-register known models for the estimator
                for label, module in KNOWN_MODELS.items():
                    TOOLKIT_DICT[toolkit_name].register_model(label, module)
            self.toolkits[toolkit_name] = TOOLKIT_DICT[toolkit_name](**self.toolkit_options.get(toolkit_name, {}))
        return self.toolkits[toolkit_name]

# ----------------------------------------------------------------------------
# Scenario model
# ----------------------------------------------------------------------------
@dataclass
class StageSpec:
    name: str
    op: str
    tool: Optional[str]
    needs: List[str] = field(default_factory=list)
    params: Dict[str, Any] = field(default_factory=dict)

@dataclass
class Scenario:
    version: int
    globals: Dict[str, Any] = field(default_factory=dict)
    stages: Dict[str, StageSpec] = field(default_factory=dict)

    @staticmethod
    def from_mapping(cfg: Dict[str, Any]) -> "Scenario":
        if "version" not in cfg:
            raise ValueError("Scenario must include a 'version' field")
        g = cfg.get("globals", {}) or {}
        stages_cfg = cfg.get("stages", {}) or {}
        stages: Dict[str, StageSpec] = {}
        for name, sc in stages_cfg.items():
            stages[name] = StageSpec(
                name=name,
                op=sc["op"],
                tool=sc.get("tool"),
                needs=list(sc.get("needs", []) or []),
                params=dict(sc.get("params", {}) or {}),
            )
        return Scenario(version=int(cfg["version"]), globals=g, stages=stages)

# ----------------------------------------------------------------------------
# Topological sorting with cycle detection (Kahn's algorithm)
# ----------------------------------------------------------------------------

def topo_order(stages: Dict[str, StageSpec]) -> List[str]:
    indeg = {k: 0 for k in stages}
    graph: Dict[str, List[str]] = {k: [] for k in stages}
    for k, spec in stages.items():
        for parent in spec.needs:
            if parent not in stages:
                raise KeyError(f"Unknown parent '{parent}' required by stage '{k}'")
            indeg[k] += 1
            graph[parent].append(k)
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

# ----------------------------------------------------------------------------
# Loader for YAML or JSON scenarios (YAML optional)
# ----------------------------------------------------------------------------

def load_scenario_file(path: Path) -> Scenario:
    text = Path(path).read_text()
    data: Dict[str, Any]
    try:
        import yaml  # type: ignore
        data = yaml.safe_load(text)
        if not isinstance(data, dict):
            raise TypeError("YAML scenario must map to a dict at the top level")
    except Exception:
        data = json.loads(text)
        if not isinstance(data, dict):
            raise TypeError("JSON scenario must map to an object at the top level")
    return Scenario.from_mapping(data)

# ----------------------------------------------------------------------------
# Pipeline orchestrator
# ----------------------------------------------------------------------------
class ScenarioPipeline:
    def __init__(self, scenario: Scenario):
        self.scenario = scenario
        self.ctx = Context()
        # Preload toolkit options from globals (if provided)
        tk_opts = scenario.globals.get("toolkit_options", {}) or {}
        if not isinstance(self.ctx.toolkit_options, defaultdict):
            self.ctx.toolkit_options = defaultdict(dict, tk_opts)
        else:
            self.ctx.toolkit_options.update(tk_opts)

    @classmethod
    def from_file(cls, path: Path) -> "ScenarioPipeline":
        return cls(load_scenario_file(path))

    def _stage_iter(self, targets: Optional[Iterable[str]], skip: Optional[Iterable[str]]) -> List[str]:
        order = topo_order(self.scenario.stages)
        targets_set = set(targets) if targets else set(order)
        unknown_targets = targets_set - set(order)
        if unknown_targets:
            raise KeyError(f"Unknown target stages: {sorted(unknown_targets)}")
        skip_set = set(skip) if skip else set()
        return [s for s in order if s in targets_set and s not in skip_set]

    def run(self, targets: Optional[Iterable[str]] = None, skip: Optional[Iterable[str]] = None) -> Iterable[Tuple[str, CrystalDataset]]:
        """Run the pipeline. Yields (stage_name, CrystalDataset)."""
        for name in self._stage_iter(targets, skip):
            spec = self.scenario.stages[name]
            logger.info("Running stage: %s (op=%s) ...", name, spec.op)
            inputs = {p: self.ctx.outputs[p] for p in spec.needs}
            if spec.op not in OP_REGISTRY:
                raise KeyError(f"Operation '{spec.op}' is not registered")
            out = OP_REGISTRY[spec.op](self.ctx, inputs, dict(spec.params))
            # Minimal provenance
            out.metadata = dict(out.metadata or {})
            out.metadata["pipeline_stage"] = name
            out.parent_ids = [ds.dataset_id for ds in inputs.values() if hasattr(ds, "dataset_id")]
            self.ctx.outputs[name] = out
            logger.info("Finished stage: %s", name)
            yield name, out

# ----------------------------------------------------------------------------
# Concrete operations (adapted from your current PPPipeline.run logic)
# ----------------------------------------------------------------------------

@op("parse_directory")
def op_parse_directory(ctx: Context, inputs: Dict[str, CrystalDataset], params: Dict[str, Any]) -> CrystalDataset:
    assert not inputs, "parse_directory takes no parents"
    elements: Optional[Set[str]] = set(params.get("elements", [])) if params.get("elements") else None
    io: StructureDatasetIO = ctx.get_tool("structure_parser")
    if elements is None:
        raise AssertionError("'elements' must be provided either in params or toolkit_options")
    ds = io.load_from_directory(elements=elements)
    assert len(ds) > 0, "Empty dataset parsed from directory"
    ds.metadata = dict(ds.metadata or {})
    ds.metadata["message"] = f"loaded {len(ds)} structures"
    return ds

@op("symmetrize_dataset")
def op_symmetrize(ctx: Context, inputs: Dict[str, CrystalDataset], params: Dict[str, Any]) -> CrystalDataset:
    assert set(inputs.keys()) == {"parse_raw"} or len(inputs) == 1, "symmetrize_dataset expects exactly one parent"
    parent = next(iter(inputs.values()))
    symtool: SymmetryToolkit = ctx.get_tool("symmetry")
    return symtool.get_symmetrized_dataset(parent)

@op("poll_databases")
def op_poll_db(ctx: Context, inputs: Dict[str, CrystalDataset], params: Dict[str, Any]) -> CrystalDataset:
    assert set(inputs.keys()) == {"parse_raw"} or len(inputs) == 1, "poll_databases expects parse_raw (or single parent)"
    parse_ds = next(iter(inputs.values()))
    elements = getattr(parse_ds, "elements", None)
    max_ehull = params.get("max_ehull", MAX_EHULL_PA)
    elements_refdata_cache = params.get("elements_refdata_cache")

    db = poll_databases(elements, cache_base_path=elements_refdata_cache, max_ehull=max_ehull)

    if params.get("estimate_energies"):
        estimator = ctx.get_tool("estimator")
        estimated = estimator.estimate_dataset_energies(db)
        # merge messages
        msg0 = (db.metadata or {}).get("message", "").strip()
        msg1 = (estimated.metadata or {}).get("message", "").strip()
        estimated.metadata = dict(estimated.metadata or {})
        estimated.metadata["message"] = (msg0 + (". " if msg0 and msg1 else "") + msg1)
        estimated.parent_ids = []
        return estimated
    return db

@op("augment_with_reference")
def op_augment(ctx: Context, inputs: Dict[str, CrystalDataset], params: Dict[str, Any]) -> CrystalDataset:
    # expects two parents: symmetrized set and polled db
    if len(inputs) != 2:
        raise AssertionError("augment_with_reference expects two parents: symmetrize_raw and poll_db")
    # Resolve by names if present
    sym_ds = inputs.get("symmetrize_raw") or next(iter(inputs.values()))
    db_ds = inputs.get("poll_db") or list(inputs.values())[1]

    similarity_tk: SimilarityTools = ctx.get_tool("similarity")
    augmentation = similarity_tk.get_unseen_in_ref(sym_ds, ref_ds=db_ds)

    # Label reproduced entries inside DB, preserve DB parameters
    labeled_entries = []
    for e in db_ds:
        if str(e.id) in augmentation.metadata.get("reproduced", {}):
            new_metadata = {**(e.metadata or {}), **{"reproduced": True}}
            labeled_entries.append(e.copy_with(**{"metadata": new_metadata}))
        else:
            labeled_entries.append(e)
    db_parameters = db_ds.__dict__.copy()
    meta = db_parameters.pop("metadata", {})
    db_parameters["base_path"] = db_parameters.pop("_base_path", None)
    to_drop = [k for k in db_parameters if k.startswith("_")]
    for k in to_drop:
        db_parameters.pop(k, None)
    labeled_db = CrystalDataset(labeled_entries, **db_parameters)
    labeled_db.metadata = meta

    merged = labeled_db.merge(augmentation)
    merged.metadata = dict(merged.metadata or {})
    merged.metadata["message"] = "Merged DB data into generated set"
    merged.parent_ids = [getattr(sym_ds, "dataset_id", None), getattr(db_ds, "dataset_id", None)]
    merged.metadata["reproducibility"] = augmentation.metadata.get("reproducibility")
    return merged

@op("estimate_energies")
def op_estimate(ctx: Context, inputs: Dict[str, CrystalDataset], params: Dict[str, Any]) -> CrystalDataset:
    if len(inputs) != 1:
        raise AssertionError("estimate_energies expects a single parent")
    parent = next(iter(inputs.values()))
    estimator = ctx.get_tool("estimator")
    return estimator.estimate_dataset_energies(parent)

@op("filter_by_hull")
def op_filter_hull(ctx: Context, inputs: Dict[str, CrystalDataset], params: Dict[str, Any]) -> CrystalDataset:
    if not inputs:
        raise AssertionError("filter_by_hull requires at least one parent")
    max_ehull = params.get("max_ehull", MAX_EHULL_PA)

    # Pick PD base dataset
    base_stage = params.get("base_stage")
    if base_stage:
        if base_stage not in inputs:
            raise KeyError(f"base_stage '{base_stage}' not among parents {list(inputs)}")
        base_ds = inputs[base_stage]
    else:
        base_ds = next(iter(inputs.values()))

    # Phase diagram tools need a dataset at construction time
    ctx.toolkit_options["phase_diag"] = {"dataset": base_ds}
    pd_tk: PhaseDiagramTools = ctx.get_tool("phase_diag")

    # We filter the (sole) child set we want to select from – by default use the only parent
    select_from = next(iter(inputs.values())) if len(inputs) == 1 else inputs.get("estimate", base_ds)

    filtered = select_from.filter(lambda entry: pd_tk.height_above_hull_pa(entry) < max_ehull)
    filtered.metadata = dict(filtered.metadata or {})
    filtered.metadata["message"] = f"Selected entries with e_hull < {max_ehull:.3f}"
    return filtered

@op("deduplicate")
def op_dedup(ctx: Context, inputs: Dict[str, CrystalDataset], params: Dict[str, Any]) -> CrystalDataset:
    if len(inputs) != 1:
        raise AssertionError("deduplicate expects a single parent")
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

# Placeholder for extension point, if you later add a dedicated DFT postprocess
@op("postprocess_dft")
def op_postprocess_dft(ctx: Context, inputs: Dict[str, CrystalDataset], params: Dict[str, Any]) -> CrystalDataset:
    raise NotImplementedError("Implement DFT postprocessing here as needed")
