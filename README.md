# vsbtools

`vsbtools` is a Python research toolbox centered on crystal-structure dataset
workflows used by the MatterGen generation postprocessing and guidance-analysis
pipeline.

## Package Map

| Module | Role |
| --- | --- |
| `vsbtools/materials_dataset` | Dataset model, staged processing, database comparison, ML energy estimation/relaxation, stability filtering, deduplication, and generation-quality analysis. |

Helper code needed by this workflow is kept inside `materials_dataset`, including
geometry checks, MatterGen parsing helpers, plotting helpers, cache/path
configuration, GRACE/MatterSim subprocess clients, and small internal utilities.

Common external names used by the workflow:

- MatterGen is a crystal generative model; `scout-matter` is the MatterGen fork
  used by the guidance workflow (<https://github.com/link-lab3629/scout-matter>).
- MatterSim and GRACE are ML interatomic-potential/energy-estimation backends;
  GRACE is accessed through `tensorpotential`.
- USPEX is an evolutionary crystal-structure-prediction code. `materials_dataset`
  retains a legacy-compatible structural fingerprint and optional bridge where
  needed for dataset deduplication.

## Crystal-Structure Dataset Workflows

`vsbtools/materials_dataset` processes crystal-structure
datasets from generation runs, simulations, database exports, or user-supplied
POSCAR, CIF, or extended XYZ (`extxyz`) collections.

The pipeline layer can clean structures, symmetrize them, compare against
reference databases, estimate or relax energies with ML models, filter by
stability, remove structural duplicates, and evaluate generation quality.

### Main Concepts

| File | Role |
| --- | --- |
| `crystal_entry.py` | `CrystalEntry`: one crystal structure plus id, optional energy, formula, and metadata. |
| `crystal_dataset.py` | `CrystalDataset`: read-oriented collection of entries with metadata, provenance, `merge()`, and `filter()`. |
| `analysis/scenario_pipeline.py` | Directed-acyclic-graph executor for cleanup, database polling, energy estimation, hull filtering, and deduplication stages. |
| `analysis/similarity_tools.py` | Structural deduplication workflows. |
| `analysis/phase_diagram_tools.py` | Phase diagram and energy-above-hull calculations. |
| `energy_estimation/` | Bridges to MatterSim and GRACE through `NNEstimator`. |
| `scripts/poll_databases.py` | Integrates reference structures from Alexandria, OQMD, Materials Project, and OPTIMADE providers. |

### Scenario Pipeline

`analysis/scenario_pipeline.py` runs YAML/JSON-defined workflows as a directed
acyclic graph of stages. Registered operations include raw parsing, density and
minimum-distance filtering, symmetrization, database polling, reference merging,
energy estimation, relaxation, hull filtering, and deduplication.

Reference database polling prefers OPTIMADE by default. Other supported sources
include local Alexandria JSON exports, local OQMD MySQL deployments, and
Materials Project through its API.

Database caches use the platform's per-user cache location and are created only
when data is written:

- Windows: `%LOCALAPPDATA%\vsbtools\Cache\DB_caches`
- macOS: `~/Library/Caches/vsbtools/DB_caches`
- Linux/Unix: `${XDG_CACHE_HOME:-~/.cache}/vsbtools/DB_caches`

Set `VSBTOOLS_CACHE_DIR` to override the `vsbtools` cache root on any platform.

### Reproducibility Notebook

A packaged reproducibility pipeline is provided for the MatterGen guidance
analysis workflow. It starts from the raw-generation archives in
`Examples/raw_generations`, postprocesses them into staged datasets, builds
summary/Pareto artifacts, plots descriptor distributions, and writes a run
manifest.

```text
vsbtools/materials_dataset/Examples/
├── README_reproducibility.md
├── setup_reproducibility_envs.sh
└── mg_generation_postprocessing_pipeline.ipynb
```

The setup script creates a contained workspace with three virtual environments:
`vsbtools`, `scout-matter`/MatterGen, and GRACE/`tensorpotential`.

```bash
bash vsbtools/materials_dataset/Examples/setup_reproducibility_envs.sh \
  --root ./vsbtools_reproducibility_env
```

See `Examples/README_reproducibility.md` for manual configuration, pinned
commit/tag setup, and launch instructions.

### Programmatic Use

```python
from pathlib import Path
from vsbtools.materials_dataset.analysis.scenario_pipeline import ScenarioPipeline
from vsbtools.materials_dataset.io import write

pipeline = ScenarioPipeline.from_file(Path("scenario.yaml"))
pipeline.ctx.globals["elements"] = ["Si", "O"]
pipeline.ctx.toolkit_options["structure_parser"]["root"] = Path("input_structures")

for stage_name, dataset in pipeline.run():
    write(dataset, enforce_base_path=Path("processed") / stage_name)
```

## Installation

Python 3.10 or newer is required.

```bash
python3 -m pip install -e .
python3 -m pip install -e ".[materials_dataset]"
```

Additional external requirements depend on the workflow:

| Workflow | External dependency |
| --- | --- |
| Optional legacy `USPEXBridge` similarity backend | Python USPEX package/installation importable as `USPEX`. |
| Database polling | Materials Project credentials and/or Alexandria/OQMD access. |
| ML energy estimation | MatterSim and/or GRACE environment. |
| Diffusion guidance analysis | MatterGen importability via `MATTERGEN_PYTHON_PATH` or host-specific configuration. |
| Packaged reproducibility notebook | Use `Examples/setup_reproducibility_envs.sh` to create contained `vsbtools`, `scout-matter`, and GRACE environments. |

## Running Tests

There is no single project test runner configured. The repository contains
`unittest`-style tests under `materials_dataset` subdirectories.

```bash
python3 -m unittest discover -s vsbtools/materials_dataset -t . -p "*_Test.py"
python3 -m unittest discover -s vsbtools/materials_dataset -t . -p "*_test.py"
```

Many tests require optional dependencies or local sample data.

## Development Notes

- Prefer the active `materials_dataset` API for new dataset-processing code.
- Many scripts assume an editable install with the repository root on `PYTHONPATH`.
- Generated artifacts such as `manifest.yaml`, `data.csv`, `POSCARS/`, cached database files, plots, and summary tables are expected outputs of normal workflows.
