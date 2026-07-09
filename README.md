# vsbtools

`vsbtools` is a Python research toolbox for atomistic materials workflows. It is
organized as several mostly independent subsystems rather than as one monolithic
application.

## Overview

The main workflow areas are:

- Crystal-structure datasets: staged pipelines for parsing generated or
  simulated structures, cleaning and symmetrizing them, comparing against
  reference databases, estimating or relaxing energies with ML models, filtering
  by stability, deduplicating structures, and analyzing generation quality. The
  MatterGen guidance-analysis tools and reproducibility notebook live inside
  this subsystem and build on the same staged-dataset representation.
- Gaussian calculation management: input generation, batch submission,
  monitoring, restart/error handling, and postprocessing of completed
  calculations from the Gaussian quantum-chemistry software package.
- Materials utilities: geometry checks, format conversion, USPEX-related helpers,
  stability-analysis routines, plotting helpers, and general filesystem,
  clustering, duplicate-analysis, and Pareto-front utilities.

## Package Map

| Module | Role |
| --- | --- |
| `vsbtools/materials_tools/materials_dataset` | Workflow layer for crystal-structure datasets: staged processing, database comparison, ML energy estimation/relaxation, stability filtering, deduplication, and generation-quality analysis. |
| `vsbtools/gaussian_calc_manager` | Independent Gaussian automation layer: input generation, submission, monitoring, restart/error recovery, and result postprocessing. |
| `vsbtools/materials_tools/visualisation_utils` | Plotting helpers for spectra, heatmaps, density-of-states and inverse-participation-ratio views, ternary phase diagrams, and figure formatting. |
| `vsbtools/materials_tools/tools_stability` | Stability-analysis utilities for binding, exchange, and fragmentation energy data. |
| `vsbtools/materials_tools/geom_utils` | Geometry checks, coordination metrics, primitive-cell extraction, structure sanity checks, and POSCAR modification helpers. |
| `vsbtools/materials_tools/ext_software_io` | Format conversion and parsers for POSCAR, XYZ, CIF, and Gaussian outputs. |
| `vsbtools/materials_tools/uspex_toolkit` | USPEX I/O, duplicate removal, template preparation, and fingerprint-distance integration. |
| `vsbtools/materials_tools/NN_energy_estimators` | Process/stream helpers for MatterSim and GRACE energy/relaxation calls. |
| `vsbtools/genutils` | General helpers for filesystem work, formatting, clustering, duplicate analysis, Pareto tools, and nested-object utilities. |
| `vsbtools/uspex_gather_stat` | Early/work-in-progress code for launching and tracking repeated USPEX calculations. |

Common external names used below:

- MatterGen is a crystal generative model; `scout-matter` is the MatterGen fork used by the guidance workflow (<https://github.com/link-lab3629/scout-matter>).
- MatterSim and GRACE are ML interatomic-potential/energy-estimation backends; GRACE is accessed through `tensorpotential`.
- USPEX is an evolutionary crystal-structure-prediction code used here only by optional legacy helpers (<https://uspex-team.org/>).

## Crystal-Structure Dataset Workflows

`vsbtools/materials_tools/materials_dataset` is the main materials-workflow
subsystem. It is designed to process arbitrary crystal-structure datasets from
generation runs, simulations, database exports, or user-supplied POSCAR, CIF,
or extended XYZ (`extxyz`) collections.

The subsystem separates the data model, disk I/O, and analysis pipelines. The
pipeline layer can clean structures, symmetrize them, compare against reference
databases, estimate or relax energies with ML models, filter by stability,
remove structural duplicates, and evaluate generation quality.

### Main Concepts

| File                              | Role                                                                                                                                                                                         |
|-----------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `crystal_entry.py`                | `CrystalEntry`: one crystal structure plus id, optional energy, formula, and metadata.                                                                                                       |
| `crystal_dataset.py`              | `CrystalDataset`: read-oriented collection of entries with metadata, provenance, `merge()`, and `filter()`.                                                                                  |
| `analysis/scenario_pipeline.py`   | Directed-acyclic-graph (DAG) executor for named dataset stages such as cleanup, database polling, energy estimation, hull filtering, and deduplication.                                      |
| `analysis/similarity_tools.py`    | Reference comparison and structural deduplication workflows.                                                                                                                                 |
| `analysis/phase_diagram_tools.py` | Phase diagram and energy-above-hull calculations.                                                                                                                                            |
| `energy_estimation/`              | Bridges to MatterSim and GRACE through `NNEstimator`.                                                                                                                                        |
| `scripts/poll_databases.py`       | Integrates reference structures from Alexandria, OQMD, Materials Project and OPTIMADE providers with optional hull filtering and deduplication. OPTIMADE is the preferred default interface. |

### Scenario Pipeline

`analysis/scenario_pipeline.py` runs YAML/JSON-defined workflows as a directed
acyclic graph (DAG) of stages. Each stage consumes one or more `CrystalDataset`
outputs and produces a new dataset with stage metadata and parent provenance.

Registered operations include:

| Operation | Purpose |
| --- | --- |
| `parse_raw` | Load input structures from a configured source directory. |
| `discard_bad_density` | Remove structures with pathological density. |
| `discard_close_atoms` | Remove structures with too-short interatomic distances. |
| `symmetrize` | Symmetrize structures through `SymmetryToolkit`. |
| `poll_db` | Fetch/cache reference structures from materials databases. |
| `merge_base_into_ref` | Compare generated/input structures against reference data and label reproduced entries. |
| `estimate` | Estimate energies using a registered ML model. |
| `relax` | Relax structures using a registered ML model. |
| `filter_hull` | Keep entries below a configured energy above hull. |
| `deduplicate` | Remove structural duplicates using the configured structural-distance backend. |

Reference database polling prefers OPTIMADE by default. OPTIMADE (Open
Databases Integration for Materials Design, <https://www.optimade.org/>) is a
common API specification for querying materials databases. `poll_databases()`
uses `pref_db="op"` and includes `optimade` in the default database list, so
OPTIMADE data is used as the reference source before additional local/provider
sources are merged in. Per-scenario configuration can still override `pref_db`
or pass provider-specific loader options. Other supported sources include the
Alexandria database (<https://alexandria.icams.rub.de/>) requiring local copies of json files, OQMD (Open Quantum
Materials Database, <https://oqmd.org/>) polling locally deployed database, and Materials Project
(<https://materialsproject.org/>) that leverages a web API.

Database caches use the platform's per-user cache location and are created only
when data is written:

- Windows: `%LOCALAPPDATA%\vsbtools\Cache\DB_caches`
- macOS: `~/Library/Caches/vsbtools/DB_caches`
- Linux/Unix: `${XDG_CACHE_HOME:-~/.cache}/vsbtools/DB_caches`

Set `VSBTOOLS_CACHE_DIR` to override the `vsbtools` cache root on any platform.

The raw-generation examples used by the notebook documentation are stored under:

```text
vsbtools/materials_tools/materials_dataset/Examples/raw_generations/
```

The rendered tutorial is kept in:

```text
vsbtools/materials_tools/materials_dataset/Doc/Crystal Dataset Use Cases.md
```

### Reproducibility Notebook

A packaged reproducibility pipeline is provided for the MatterGen guidance
analysis workflow. It starts from the raw-generation archives in
`Examples/raw_generations`, postprocesses them into staged datasets, builds
summary/Pareto artifacts, plots descriptor distributions, and writes a run
manifest.

```text
vsbtools/materials_tools/materials_dataset/Examples/
├── README_reproducibility.md
├── setup_reproducibility_envs.sh
└── mg_generation_postprocessing_pipeline.ipynb
```

The setup script creates a contained workspace with three virtual environments:
`vsbtools`, `scout-matter`/MatterGen, and GRACE/`tensorpotential`. Jupyter,
IPython, matplotlib, pip cache, and `vsbtools` external-path state are kept
under the chosen workspace root rather than under the user's global
configuration.

```bash
bash vsbtools/materials_tools/materials_dataset/Examples/setup_reproducibility_envs.sh \
  --root ./vsbtools_reproducibility_env
```

See `Examples/README_reproducibility.md` for manual configuration, pinned
commit/tag setup, and launch instructions.

### Programmatic Scenario Pipelines

The reproducibility notebook uses the same lower-level `ScenarioPipeline`
interface that can also be called directly from Python. A minimal pipeline run
looks like this:

```python
from pathlib import Path
from vsbtools.materials_tools.materials_dataset.analysis.scenario_pipeline import ScenarioPipeline
from vsbtools.materials_tools.materials_dataset.io import write

pipeline = ScenarioPipeline.from_file(Path("scenario.yaml"))
pipeline.ctx.globals["elements"] = ["Si", "O"]
pipeline.ctx.toolkit_options["structure_parser"]["root"] = Path("input_structures")

for stage_name, dataset in pipeline.run():
    write(dataset, enforce_base_path=Path("processed") / stage_name)
```

### Structural Similarity and Deduplication

Structural comparison starts from a representation of each crystal and a metric
for comparing those representations. The default structural-distance backend uses
the radial-distribution-function fingerprint introduced in the USPEX
structural-comparison workflow, implemented through DScribe descriptors
(<https://singroup.github.io/dscribe/>). Deduplication is the downstream use of
this similarity metric: structures that are close under the chosen criterion are
grouped, and a representative is kept from each group. If energies are available,
the lowest-energy structure in the group is used as the representative.

When a distance-based backend is used, `SimilarityTools.deduplicate()` can also
expose the pairwise distance matrix and the resulting clusters. By default,
these objects are kept in memory; pass `save_dist_matrix_file=True` and/or
`save_clusters_file=True` to persist `*_dist_matrix.pkl` or `*_clusters.pkl`
files. The `check_dist_matrix_file` and `check_clusters_file` flags only request
reuse of existing files.

The similarity layer is backend-based: the default DScribe implementation can be
replaced by another structural-comparison backend, including the optional
`USPEXBridge`. The same pattern is used elsewhere in `materials_dataset`; for
example, `NNEstimator` delegates energy estimation and relaxation to registered
model bridges such as MatterSim or GRACE.

### Diffusion Guidance Analysis

Diffusion-guidance tools live in:

```text
vsbtools/materials_tools/materials_dataset/scripts/diffusion_analysis_scripts/
```

| File | Purpose |
| --- | --- |
| `guidance_stats.py` | Locate guided/non-guided processed repositories, derive guidance callables from metadata, collect stage datasets, compute target-property values, build histogram/kernel-density-estimate data, plot distributions, and run p-value helpers. |
| `guided_generation_metrics.py` | Higher-level reports for stable/new/constraint-satisfying intersections, scout ratios, p-value summaries, and artifact export. |
| `mattergen_bridge.py` | Adapter from `CrystalEntry`/pymatgen structures to MatterGen target and loss functions. |
| `pvalue_utils.py` | Exact binomial and pooled two-proportion z-test utilities for comparing good-structure fractions. |

Example: compare mean-coordination distributions across processed runs.

```python
from pathlib import Path
from vsbtools.materials_tools.materials_dataset.scripts.diffusion_analysis_scripts.guidance_stats import (
    collect_stage_dataset_dict,
    get_mean_coordination_gen_dirs,
    histo_data_collection,
    plot_multihistogram,
)

processed_root = Path("/path/to/PROCESSED")

gen_dirs = get_mean_coordination_gen_dirs(
    processed_root,
    system="Si-O",
    bond="Si-O",
    target=6,
)

datasets = collect_stage_dataset_dict(
    gen_dirs,
    stage="deduplicate",
    ref_stage="poll_db",
    add_guid_descr=True,
)

histograms = histo_data_collection(
    datasets,
    callable_name="compute_mean_coordination",
    callable_params={"type_A": 14, "type_B": 8},
    auto_adjust_bins=True,
    n_bins=12,
    integer_bins=False,
)

fig, ax = plot_multihistogram(histograms, target=6, title="Si-O mean coordination")
fig.savefig("si_o_coordination.png", dpi=300)
```

Example: full guidance metric report.

```python
from vsbtools.materials_tools.materials_dataset.scripts.diffusion_analysis_scripts.guided_generation_metrics import (
    collect_guidance_metric_report,
    write_guidance_metric_artifacts,
)

report = collect_guidance_metric_report("/path/to/processed_repo")
written_paths = write_guidance_metric_artifacts(report)
```

### Input, Output, and Persistence Utilities

The read/write utilities are general `CrystalDataset` serialization helpers. They
can be used inside scenario pipelines, but they are not tied to any particular
analysis workflow. `CrystalDataset` objects can be saved as folders containing:

```text
manifest.yaml
data.csv
POSCARS/
```

Relevant files:

| File | Role |
| --- | --- |
| `io/structures_dataset_io.py` | Load CIF, POSCAR, multi-image POSCARS, or extended XYZ (`extxyz`) files from directories. |
| `io/yaml_csv_poscars.py` | Read/write `manifest.yaml`, `data.csv`, and `POSCARS/` dataset folders. |
| `io/preset_loaders.py` | Cached loaders for external databases and USPEX/CSV/POSCARS sources. |
| `io/sources/` | Source-specific parsers and clients. |


When writing a dataset, pass the target dataset folder with `enforce_base_path`.
When reading it back, pass the path to the folder's `manifest.yaml` file.

Minimal standalone usage:

```python
from vsbtools.materials_tools.materials_dataset.io import read, write

ds = read("processed/deduplicate/manifest.yaml")
filtered = ds.filter(lambda entry: entry.energy is not None)
write(filtered, enforce_base_path="processed/with_energy")
```

## Gaussian Calculation Manager

`vsbtools/gaussian_calc_manager` is independent from the crystal-dataset
pipeline. Here, Gaussian refers to the external Gaussian quantum-chemistry
software package. The manager handles many Gaussian calculations as a task
database with polling, submission, restart, and correction logic.

Main files:

| File | Role |
| --- | --- |
| `gau_recalc.py` | Main polling/submission loop. Creates or resumes a `GauCalcDB`, updates statuses, submits pending jobs, writes logs/statistics, and stops on `DONE`, `STOP`, or time limit. |
| `src/tasks_database.py` | Defines `GauTask` and `GauCalcDB`. Handles task folders, job scripts, queue checks, Gaussian output parsing, restart/failure correction, status statistics, and database dumps. |
| `src/gjf.py` | Parser/writer/editor for Gaussian `.gjf` files. Supports recursive corrector dictionaries for route keywords and numeric options. |
| `src/ext_software_io.py` | POSCARS/XYZ/Gaussian readers and Gaussian output parsing through `cclib`. |
| `db_postproc.py` | Postprocessing helpers for completed Gaussian databases: lowest-energy structures, POSCAR/XYZ export, text tables, and spectra plotting. |
| `machines.cjson` | Host/scheduler configuration: job template, submit command, job-id regex, and queue command. |
| `scenarios*.cjson` | Gaussian scenarios: initial `.gjf`, normal-termination flags, restart strategy, known error/failure matchers, and correction rules. |
| `gjf_templates/` | Gaussian input templates. |
| `job_templates/` | Local/cluster shell templates with placeholders such as `@INPUT`, `@OUTFILE`, `@JOBNAME`, and `@NPROCSHARED`. |

Task statuses:

| Status | Meaning |
| --- | --- |
| `P` | Pending. Input/job files should be copied and submitted. |
| `R` | Running. A scheduler/process job id is known. |
| `L` | Loaded from an existing calculation folder with an output file. |
| `D` | Done. Normal termination was detected. |
| `F` | Failed. Known or unexpected failure could not be recovered. |

Typical flow:

1. Prepare `machines.cjson`, a `scenarios*.cjson`, a Gaussian template, and a geometry source (`POSCARS` or XYZ batch).
2. Run `gau_recalc.py` from `vsbtools/gaussian_calc_manager`, because the script currently imports `src.*` as a local package.
3. The script creates `results`, `results_1`, etc. unless `--recalc_folder` is provided.
4. Each task gets its own folder with a `.gjf`, job shell script, output log, checkpoint file, and generated XYZ snapshots when parseable.
5. `GauCalcDB.update()` parses logs, refreshes geometries/energies, applies scenario correctors for known errors, and marks tasks `D` or `F`.
6. `GauCalcDB.submit_jobs()` submits pending tasks up to `--maxcalcs`.
7. `database.pkl`, `energies.txt`, `stats.txt`, and the main `log` are updated in the result folder.

Example:

```bash
cd vsbtools/gaussian_calc_manager
python3 gau_recalc.py \
  -c local \
  -g POSCARS \
  -m 2 \
  -s 5 \
  --maxiter 5 \
  --machines machines.cjson \
  --scenarios scenarios.cjson
```

Key options:

| Option | Meaning |
| --- | --- |
| `--machine` | Machine key in `machines.cjson`; auto-detected from hostname if omitted. |
| `--geoms_file` | Input geometry batch, default `POSCARS`. |
| `--recalc_folder` | Existing or new calculation root. |
| `--outfile_pattern` | Gaussian output filename pattern, default `log`. |
| `--maxcalcs` | Maximum concurrently running/submitted tasks. |
| `--maxiter` | Maximum correction/restart attempts per task. |
| `--min_mult` | Derive minimal spin multiplicity from total atomic number parity. |
| `--exec_time` | Stop polling after this many minutes. |

## Supporting Materials Modules

| Path | Purpose |
| --- | --- |
| `materials_tools/NN_energy_estimators` | Process/stream helpers for MatterSim and GRACE energy/relaxation calls. |
| `materials_tools/geom_utils` | Coordination metrics, primitive-cell extraction, POSCAR modifiers, density/min-distance checks, and symmetry utilities. |
| `materials_tools/ext_software_io` | POSCAR/XYZ/Gaussian/CIF conversion and Gaussian log parsing. |
| `materials_tools/uspex_toolkit` | USPEX I/O, duplicate removal, template preparation, and fingerprint distance integration. |
| `materials_tools/uspex_001_postproc` | Utilities for older USPEX 001 result folders. |
| `materials_tools/el_spectra_tools` | Density/orbital/k-point/localization helpers. |
| `materials_tools/tools_stability` | Binding, exchange, and fragmentation energy table helpers. |
| `materials_tools/visualisation_utils` | Plot formatting, heatmaps, density-of-states and inverse-participation-ratio plots, ternary phase diagrams, and rasterization helpers. |

## General Utilities and Work-In-Progress Modules

| Path | Purpose |
| --- | --- |
| `vsbtools/genutils` | Filesystem helpers, formatting helpers, clustering, duplicate analysis, Pareto tools, and nested-object utilities. |
| `vsbtools/uspex_gather_stat` | Early/work-in-progress code for launching and tracking repeated USPEX calculations. |

## Installation

Python 3.10 or newer is required.

```bash
python3 -m pip install -e .
```

Optional extras declared in `pyproject.toml`:

```bash
python3 -m pip install -e ".[gaussian_calc_manager]"
python3 -m pip install -e ".[materials_tools]"
python3 -m pip install -e ".[genutils]"
python3 -m pip install -e ".[uspex_gather_stat]"
```

Additional external requirements depend on the workflow:

| Workflow | External dependency |
| --- | --- |
| Gaussian manager | Gaussian executable, local/cluster scheduler, `cclib`. |
| Optional legacy `USPEXBridge` similarity backend | Python USPEX package/installation importable as `USPEX`. |
| Database polling | Materials Project credentials and/or Alexandria/OQMD access. |
| ML energy estimation | MatterSim and/or GRACE environment. |
| Diffusion guidance analysis | MatterGen importability via `MATTERGEN_PYTHON_PATH` or host-specific configuration. |
| Packaged reproducibility notebook | Use `Examples/setup_reproducibility_envs.sh` to create contained `vsbtools`, `scout-matter`, and GRACE environments. |
| Plotting/reporting | Common scientific stack such as `matplotlib`, `pandas`, and `PyYAML`. |

## Running Tests

There is no single project test runner configured. The repository contains
`unittest`-style tests under several `unittests` directories.

```bash
python3 -m unittest discover -s vsbtools -p "*_Test.py"
python3 -m unittest discover -s vsbtools -p "*_test.py"
```

Many tests require optional dependencies or local sample data.

## Development Notes

- Prefer the active `materials_dataset` API for new dataset-processing code.
- Treat `materials_dataset` and `gaussian_calc_manager` as independent projects inside the same toolbox.
- Many scripts assume they are run from their own directory or from an editable install with the repository root on `PYTHONPATH`.
- Generated artifacts such as `database.pkl`, `manifest.yaml`, `data.csv`, `POSCARS/`, Gaussian logs, and summary tables are expected outputs of normal workflows.
- The root package exposes `genutils`, `gaussian_calc_manager`, `materials_tools`, and `uspex_gather_stat` through `vsbtools.__all__`.
