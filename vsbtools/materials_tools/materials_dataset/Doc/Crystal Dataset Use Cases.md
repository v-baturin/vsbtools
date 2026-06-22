Examples assume they are run from the repository root `/home/vsbat/my_git_projects/vsbtools`. The raw-generation fixtures currently live under `unittests_datasets/raw_generations`; if a local checkout uses `unittest_databases/raw_generations`, change only `RAW_GENERATIONS` in the setup cell. Each first-level zip in `RAW_GENERATION_DIR` is one MatterGen generation: `Cu-Si-P_nonguided.zip` is single-run, while the guided archives contain multiple `run_*` directories. For examples, unzip one first-level archive into a temporary directory and then call the regular `StructureDatasetIO(...).load_from_directory(...)` on the extracted directory.

```python
from pathlib import Path

REPO_ROOT = Path('/home/vsbat/my_git_projects/vsbtools')
RAW_GENERATIONS = REPO_ROOT / 'vsbtools/materials_tools/materials_dataset/unittests_datasets/raw_generations'
RAW_GENERATION_DIR = RAW_GENERATIONS / 'Cu-Si-P'
GENERATION_ZIP = RAW_GENERATION_DIR / 'Cu-Si-P_guided_CuP3.zip'
GENERATION_NAME = GENERATION_ZIP.stem
ELEMENTS = {'Cu', 'Si', 'P'}
WORK = Path('tmp_crystal_dataset_examples')
WORK.mkdir(exist_ok=True)
```

# 1. Operations with a single Dataset

The central object is `CrystalDataset`: an iterable collection of `CrystalEntry` objects with shared metadata. The `raw_ds` example dataset is loaded once in section 1.1.1 and reused by later cells.

##     1.1 Load external data into a crystal dataset

Use parser classes for local files and preset loaders for external databases.

###         1.1.1 Raw generated data

MatterGen-style raw outputs are stored as first-level generation archives. This example unzips one archive into a temporary directory and then uses the normal directory loader on the extracted contents.

```python
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from vsbtools.materials_tools.materials_dataset.io.structures_dataset_io import StructureDatasetIO


@contextmanager
def temporary_unzip(zip_path: Path):
    zip_path = Path(zip_path)
    with TemporaryDirectory(prefix=f'{zip_path.stem}_') as tmpdir:
        extracted_root = Path(tmpdir) / zip_path.stem
        extracted_root.mkdir()
        with ZipFile(zip_path) as zip_file:
            zip_file.extractall(extracted_root)
        yield extracted_root


with temporary_unzip(GENERATION_ZIP) as generation_root:
    raw_ds = StructureDatasetIO(
        generation_root,
        id_prefix=f'{GENERATION_NAME}_',
        source_name='MatterGen',
    ).load_from_directory(elements=ELEMENTS)

len(raw_ds), GENERATION_ZIP.name, sorted(raw_ds.elements)
```

###         1.1.2 Poll Databases

Database polling is optional and may need API keys, network access, or a local cache.

```python
from vsbtools.materials_tools.materials_dataset.io.preset_loaders import load_from_materials_project

# Requires an MP API key/configured environment.
reference_ds = load_from_materials_project(
    elements={'Cu', 'Si', 'P'},
    message='Cu-Si-P reference phases from Materials Project',
)
print(len(reference_ds))
```

##     1.2 Clean the Dataset

Cleaning usually combines geometry checks, duplicate removal, and optional symmetry normalization.

###         1.2.1 Remove structures with too short bonds

Filter entries by the package geometry validator.

```python
from vsbtools.materials_tools.geom_utils.structure_checks import check_min_dist_pmg

clean_ds = raw_ds.filter(lambda entry: check_min_dist_pmg(entry.structure)[0])
print(f'{len(clean_ds)} / {len(raw_ds)} structures passed the minimum-distance check')
```

###         1.2.2 Remove Duplicates

The default in-package duplicate workflow uses `SimilarityTools` with a DScribe-backed fingerprint distance.

```python
from vsbtools.materials_tools.materials_dataset.analysis.similarity_tools import SimilarityTools
from vsbtools.materials_tools.materials_dataset.analysis.structural_distance.dscribe_bridge import DScribeBridge

bridge = DScribeBridge(elements={'Cu', 'Si', 'P'}, tol_FP=0.04)
similarity = SimilarityTools(bridge.fp_dist, bridge.tol_FP)

dedup_ds, clusters, best_idx = similarity.deduplicate(
    clean_ds,
    enforce_compositions_separation=True,
)
print(f'{len(dedup_ds)} unique structures in {len(clusters)} clusters')
```

###         1.2.3 Symmetrize the Entries

Symmetrization returns a new dataset with standardized structures.

```python
from vsbtools.materials_tools.materials_dataset.analysis.symmetry_tools import SymmetryToolkit

symmetry = SymmetryToolkit(a_sym_prec=1e-3, e_sym_prec=1e-3)
sym_ds = symmetry.get_symmetrized_dataset(clean_ds)
print(len(sym_ds))
```

##     1.3 Estimate the energy

Energy-dependent analysis expects each `CrystalEntry.energy` to be populated, either from a database, DFT, or an estimator.

###         1.3.1 Construct the convex hull and calculate heights*

Once energies are present, build a phase diagram from the reference dataset and evaluate generated entries against it.

\*special case of [[Crystal Dataset Use Cases#2.1 Heights of Dataset 1 above the convex hull from Dataset 2 (generated wrt reference)|relative convex hull heights]].

```python
from vsbtools.materials_tools.materials_dataset.analysis.phase_diagram_tools import PhaseDiagramTools

# `reference_ds` and `generated_ds` must contain energies in eV per structure.
pd_tools = PhaseDiagramTools(reference_ds)
heights_pa = [pd_tools.height_above_hull_pa(entry) for entry in generated_ds]
print(heights_pa[:5])
```

##     1.4 Various I/O operations

Use the converter helpers when you need a tabular view for reporting or quick inspection.

###         1.4.1 Dataset <> Pandas DataFrame

DataFrames are convenient for filtering by scalar metadata and joining with external tables.

```python
from vsbtools.materials_tools.materials_dataset.converters import ds2df, df2ds

frame = ds2df(clean_ds)
print(frame.head())

roundtrip_ds = df2ds(frame)
print(len(roundtrip_ds))
```

##     1.5 Building reports

Summary helpers collect common scalar fields into a compact table.

###         1.5.1 Summary table (text and csv)

Write a human-readable table and a CSV from the same summary frame.

```python
from vsbtools.materials_tools.materials_dataset.analysis.summary import collect_summary_df, print_pretty_df

summary_df = collect_summary_df(clean_ds)
print_pretty_df(summary_df, dump_path=WORK / 'summary.txt', pretty=True)
print_pretty_df(summary_df, dump_path=WORK / 'summary.csv', pretty=False)
summary_df.head()
```

##     1.2 Read and write Datasets on disk

The YAML/CSV/POSCAR writer creates a portable manifest plus structure files.

```python
from vsbtools.materials_tools.materials_dataset.io.yaml_csv_poscars import read, write

write(
    clean_ds,
    enforce_base_path=WORK / 'clean_dataset',
    comment='Minimal cleaned raw-generation fixture',
)
manifest = WORK / 'clean_dataset' / 'manifest.yaml'
restored_ds = read(manifest)
print(len(restored_ds), restored_ds.base_path)
```

# 2. Operations with two Datasets

Two-dataset workflows usually compare generated structures with a reference hull or reference objective values.

##     2.1 Heights of Dataset 1 above the convex hull from Dataset 2 (generated wrt reference)

Build the hull on the reference dataset, then evaluate the generated dataset entry by entry.

```python
reference_hull = PhaseDiagramTools(reference_ds)
relative_heights = {
    entry.id: reference_hull.height_above_hull_pa(entry)
    for entry in generated_ds
}
```

##     2.2 Constructing and plotting of Pareto fronts

Pareto front plotting expects `pf_*.csv` files from a postprocessing stage. This toy CSV mirrors the stage output shape.

```python
import pandas as pd
from vsbtools.materials_tools.materials_dataset.analysis.pareto_fronts import plot_pareto_fronts

pareto_stage = WORK / 'pareto_demo' / '1_add_ref_deduplicated'
pareto_stage.mkdir(parents=True, exist_ok=True)
pd.DataFrame(
    {
        'id': ['0', '1', 'ref_A'],
        'composition': ['Cu1 Si1 P1', 'Cu2 P1', 'Cu1 P3'],
        'loss': [0.08, 0.15, 0.0],
        'e_hull/at': [0.21, 0.12, 0.0],
    }
).to_csv(pareto_stage / 'pf_1.csv', index=False)

ax = plot_pareto_fronts(pareto_stage, n_fronts=1, max_loss=None, max_ehull=None, title='Cu-Si-P demo')
ax.figure.savefig(WORK / 'pareto_demo.png', dpi=200)
```

# 3. Operations with multiple datasets

Represent multiple datasets as a dictionary keyed by a short label.

```python
datasets = {
    'raw': raw_ds,
    'clean': clean_ds,
    'symmetrized': sym_ds,
}
{k: len(v) for k, v in datasets.items()}
```

##     3.1 Plot distributions of a descriptor

Descriptor distributions can be plotted from arrays, or collected from stage datasets in a processed generation repository.

###         3.1.1 If guided generation ex

For raw fixtures, classify generations by first-level archive name. To compare generations, unzip each archive temporarily and call the same directory loader on each extracted generation.

```python
generation_archives = sorted(RAW_GENERATION_DIR.glob('*.zip'))
guided_archives = [path for path in generation_archives if '_guided_' in path.stem]
nonguided_archives = [path for path in generation_archives if '_nonguided' in path.stem]
print('guided:', [path.name for path in guided_archives])
print('non-guided:', [path.name for path in nonguided_archives])

datasets_by_generation = {}
for archive in generation_archives:
    with temporary_unzip(archive) as generation_root:
        datasets_by_generation[archive.stem] = StructureDatasetIO(
            generation_root,
            id_prefix=f'{archive.stem}_',
            source_name='MatterGen',
        ).load_from_directory(elements=ELEMENTS)

{k: len(v) for k, v in datasets_by_generation.items()}
```

###         3.1.2 Histograms of a descriptor

This example plots the number of atoms per structure for each dataset.

```python
import numpy as np
from vsbtools.materials_tools.materials_dataset.analysis.guidance_statistics import values_2_histo_data, plot_multihistogram

histograms = []
for label, dataset in datasets.items():
    values = np.array([entry.natoms for entry in dataset])
    bin_centers, counts = values_2_histo_data(values, integer_bins=True)
    histograms.append({'label': label, 'bin_centers': bin_centers, 'counts': counts})

fig, ax = plot_multihistogram(histograms, xlabel='atoms per structure')
fig.savefig(WORK / 'natoms_histogram.png', dpi=200)
```

###         3.1.3 Kernel density estimation plots

KDE plots are better for continuous descriptors. Here, volume per atom is used as a fixture-friendly descriptor.

```python
from vsbtools.materials_tools.materials_dataset.analysis.guidance_statistics import plot_multi_kde

volume_pa = {
    label: np.array([entry.structure.volume / entry.natoms for entry in dataset])
    for label, dataset in datasets.items()
}
fig, ax = plot_multi_kde(volume_pa, xlabel='volume per atom / A^3')
fig.savefig(WORK / 'volume_pa_kde.png', dpi=200)
```

# 4. Pipeline for generation postprocessing

The scenario pipeline makes the usual parse, clean, reference, estimate, hull-filter, and deduplicate steps reproducible.

##     4.1 Scenarios syntax

A scenario is a small DAG. Keep external stages commented or skipped until the required optional tools are configured.

```python
from vsbtools.materials_tools.materials_dataset.analysis.scenario_pipeline import ScenarioPipeline

scenario_file = REPO_ROOT / 'vsbtools/materials_tools/materials_dataset/analysis/unittests/scenario_postprocess.yaml'
pipeline = ScenarioPipeline.from_file(scenario_file)
pipeline.ctx.globals['elements'] = sorted(ELEMENTS)
pipeline.ctx.toolkit_options['structure_parser'].pop('batch_metadata_file', None)

with temporary_unzip(GENERATION_ZIP) as generation_root:
    pipeline.ctx.toolkit_options['structure_parser']['root'] = generation_root
    partial_outputs = dict(pipeline.run(targets=['parse_raw', 'symmetrize_raw']))

partial_outputs.keys()
```

Minimal scenario shape:

```yaml
version: 1
globals:
  toolkit_options:
    structure_parser:
      source_name: MatterGen
stages:
  parse_raw:
    op: parse_raw
    needs: []
    params: {}
  symmetrize_raw:
    op: symmetrize
    needs: [parse_raw]
    params: {}
```

# 5. Standard course of action

A compact working sequence is: postprocess the generation, plot descriptor distributions, then plot Pareto fronts.

```python
# 1. Postprocess raw generations into staged datasets.
# Full runs use the same temporary-unzip generation-root convention once DB/estimator tools are configured.
postprocessed = partial_outputs['symmetrize_raw']

# 2. Plot histograms/KDE for descriptors that matter for the generation target.
# Reuse `plot_multihistogram(...)` and `plot_multi_kde(...)` from section 3.

# 3. Plot Pareto fronts from the stage that contains pf_*.csv files.
# plot_pareto_fronts(WORK / 'postprocess_repo/Cu-Si-P/1_add_ref_deduplicated')
```
