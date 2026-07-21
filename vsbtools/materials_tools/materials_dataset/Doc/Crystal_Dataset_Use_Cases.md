Examples assume `vsbtools` is importable from the checkout being documented. The raw-generation fixtures live under `Examples/raw_generations`. Each first-level zip in `RAW_GENERATION_DIR` is one MatterGen generation: `Cu-Si-P_nonguided.zip` is single-run, while the guided archives contain multiple `run_*` directories. For examples, unzip one first-level archive into the tutorial work directory, call the regular `StructureDatasetIO(...).load_from_directory(...)` on the extracted directory, then delete the extracted files when the loaded dataset is no longer using them.


```python
from pathlib import Path
import sysconfig
import vsbtools

PACKAGE_ROOT = Path(vsbtools.__file__).parent
RAW_GENERATIONS = PACKAGE_ROOT / 'materials_tools/materials_dataset/Examples/raw_generations'
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

Here we unzip the results of multi-run guided generation, stored in `GENERATION_ZIP`, use the normal directory-based loader to transform them into a `CrystalDataset` object, and then delete the extracted archive contents.


```python
from pathlib import Path
from shutil import rmtree
from zipfile import ZipFile

from vsbtools.materials_tools.materials_dataset.io.structures_dataset_io import StructureDatasetIO


# 1. Choose where the archive will be extracted for this tutorial run.
generation_root = WORK / GENERATION_NAME

# 2. Start from a clean extraction directory so repeated notebook runs are reproducible.
if generation_root.exists():
    rmtree(generation_root)
generation_root.mkdir(parents=True)

# 3. Unzip the generated structures into the visible work directory.
zip_file = ZipFile(GENERATION_ZIP)
zip_file.extractall(generation_root)
zip_file.close()

# 4. Load the extracted generation exactly as a normal directory-based dataset.
raw_ds = StructureDatasetIO(
    generation_root,
    id_prefix=f'{GENERATION_NAME}_',
    source_name='MatterGen',
).load_from_directory(elements=ELEMENTS)

# 5. Delete the extracted files after loading; the dataset object is already in memory.
rmtree(generation_root)

len(raw_ds), GENERATION_ZIP.name, sorted(raw_ds.elements)

```




    (708, 'Cu-Si-P_guided_CuP3.zip', ['Cu', 'P', 'Si'])




```python
raw_ds[0]
```




    CrystalEntry(id='Cu-Si-P_guided_CuP3_0', structure=Structure Summary
    Lattice
        abc : 3.502399882181507 3.501149695799593 6.965021133422852
     angles : 75.45851524770619 75.42161054589617 59.891605704836515
     volume : 70.70674452845024
          A : 3.3896372318267822 0.0 0.8815692663192749
          B : 1.5861077546831546 2.9949196972279606 0.8790718913078308
          C : 0.0 0.0 6.965021133422852
        pbc : True True True
    PeriodicSite: Si (4.674, 2.987, 2.307) [0.9123, 0.9975, 0.08991]
    PeriodicSite: Cu (4.1, 2.644, 4.51) [0.7967, 0.8827, 0.4353]
    PeriodicSite: Si (2.181, 1.488, 4.931) [0.4109, 0.4967, 0.5933]
    PeriodicSite: P (1.62, 1.15, 7.081) [0.2982, 0.384, 0.9304], energy=None, formula='Si2 Cu1 P1', metadata={'source': 'MatterGen', 'file': 'generated_crystals.extxyz'})



###         1.1.2 Poll Databases

Database polling is optional and may need API keys, network access, or a local cache.

#### 1.1.2.1 Polling locally deployed databases

Requires local Alexandria json files and OQMD database


```python
from vsbtools.materials_tools.materials_dataset.io.preset_loaders import load_from_oqmd

# Requires an OQMD API configured environment.
reference_ds = load_from_oqmd(
    elements={'Cu', 'Si', 'P'},
    message='Cu-Si-P reference phases from Materials Project',
)
print(len(reference_ds))
```

    257



```python
# MAY TAKE UP TO 15-20 minutes to poll the jsons!!!
from vsbtools.materials_tools.materials_dataset.io.preset_loaders import load_from_alexandria

# Requires an Alexandria json files
reference_ds = load_from_alexandria(
    elements={'Cu', 'Si', 'P'},
    message='Cu-Si-P reference phases from Materials Project', force_refresh=True,
)
print(len(reference_ds))
```

    Parsing /home/vsbat/work/Alexandria/pbe_data/alexandria_050.json ... (51/51)
    Data saved to /home/vsbat/.cache/vsbtools/DB_caches/alexandria_Cu-P-Si_none
    750


#### 1.2.1.2 Polling web-interfaces

Requires internet access and working servers of materials project or optimade ports of alexandria and oqmd
**Optimade is preferred over local sources** due to the access to the latest versions of data (includes the access to the latest MaterialsProject data as well) and quicker sql-based server-side poll of alexandria database


```python
from vsbtools.materials_tools.materials_dataset.io.preset_loaders import load_from_materials_project

# Requires an MP API key/configured environment.
reference_ds = load_from_materials_project(
    elements={'Cu', 'Si', 'P'},
    message='Cu-Si-P reference phases from Materials Project',
)
print(len(reference_ds))

```

    87


**Polling using OPTIMADE interface**


```python
# This request will poll and deduplicate info sequentiall from Alexandria, Materials Project and OQMD.
from vsbtools.materials_tools.materials_dataset.io.preset_loaders import load_from_optimade

providers = ("mp", "oqmd")  #, "alexandria")
reference_ds = load_from_optimade(
    elements={'Cu','Si', 'P'},
    message=f'Cu-Si-P from {", ".join(providers)}',
    providers=providers,
    timeout=1500, force_refresh=True)

print(len(reference_ds))
```

##     1.2 Clean the Dataset

Cleaning usually combines geometry checks, duplicate removal, and optional symmetry normalization.

###         1.2.1 Remove structures with too short bonds

Filter entries by the package geometry validator.


```python
from vsbtools.materials_tools.materials_dataset.geom_utils.structure_checks import check_min_dist_pmg

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

symmetry = SymmetryToolkit(a_sym_prec=1e-3,  # precision of analysis
                           e_sym_prec=1e-2)  # precision of enforcing a symmetry
sym_ds = symmetry.get_symmetrized_dataset(clean_ds)

changed_symmetry_count = sum(
    symmetry.sym_group_no(original) != symmetry.sym_group_no(symmetrized)
    for original, symmetrized in zip(clean_ds, sym_ds)
)
print(f'{changed_symmetry_count} / {len(sym_ds)} structures changed symmetry after symmetrization')

```

##     1.3 Estimate the energy

Energy-dependent analysis expects each `CrystalEntry.energy` to be populated, either from a database, DFT, or an estimator.

###         1.3.1 Energy estimation using GRACE


```python
from vsbtools.materials_tools.materials_dataset.energy_estimation.nn_estimator import NNEstimator
from vsbtools.materials_tools.materials_dataset.energy_estimation import grace_bridge

NNEstimator.register_model("grace", grace_bridge)
estimator = NNEstimator(default_model="grace")
reference_ds = estimator.estimate_dataset_energies(reference_ds)
```

##     1.4 Construct the convex hull and calculate heights*

Once energies are present, build a phase diagram from the reference dataset and evaluate generated entries against it.

\*special case of **2.1 Heights of Dataset 1 above the convex hull from Dataset 2 (generated wrt reference)**.


```python
from vsbtools.materials_tools.materials_dataset.analysis.phase_diagram_tools import PhaseDiagramTools

# `reference_ds` and `generated_ds` must contain energies in eV per structure.
pd_tools = PhaseDiagramTools(reference_ds)
heights_pa = [pd_tools.height_above_hull_pa(entry) for entry in reference_ds]
print(heights_pa[:5])

```

##     1.5 Various I/O operations

Use the converter helpers when you need a tabular view for reporting or quick inspection.

###         1.5.1 Dataset <> Pandas DataFrame

DataFrames are convenient for filtering by scalar metadata and joining with external tables.


```python
from vsbtools.materials_tools.materials_dataset.converters import ds2df, df2ds

frame = ds2df(clean_ds)
print(frame.head())

roundtrip_ds = df2ds(frame)
print(len(roundtrip_ds))

```

##     1.6 Building reports

Summary helpers collect common scalar fields into a compact table.

###         1.6.1 Summary table (text and csv)

Write a human-readable table and a CSV from the same summary frame.


```python
from vsbtools.materials_tools.materials_dataset.analysis.summary import collect_summary_df, print_pretty_df

summary_df = collect_summary_df(clean_ds)
summary_txt = WORK / 'summary.txt'
summary_csv = WORK / 'summary.csv'

print_pretty_df(summary_df, dump_path=summary_txt, pretty=True)
print_pretty_df(summary_df, dump_path=summary_csv, pretty=False)

for summary_path in (summary_txt, summary_csv):
    print(f'\n{summary_path.name}')
    for line in summary_path.read_text().splitlines()[:12]:
        print(line)

```

##     1.7 Read and write Datasets on disk

The YAML/CSV/POSCARS.zip writer creates a portable manifest plus zipped structure files.
Dataset is read via its `manifest.yaml` file


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
# Let's first load and estimate the generated dataset

```


```python
reference_hull = PhaseDiagramTools(reference_ds)
# generated_ds = estimator.estimate_dataset_energies(sym_ds)
relative_heights = {
    entry.id: reference_hull.height_above_hull_pa(entry)
    for entry in generated_ds
}
import random
print(random.sample(list(relative_heights.items()), 5))
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

##     3.1 Plot average coordination-number distributions

This example keeps the local `datasets` dictionary from section 3 and computes the average coordination number in the same way as the environment-guidance plots in the reference notebook.


```python
import matplotlib.pyplot as plt
import numpy as np
from pymatgen.core import Element

from vsbtools.materials_tools.materials_dataset.analysis.guidance_statistics import (
    calculate_values,
    get_target_value_fn,
    plot_multi_kde,
    plot_multihistogram,
    values_2_histo_data,
)
```

###         3.1.1 Choose the coordination environment

Set `bond` and `target` exactly as in the environment-guidance notebooks. Optional entries in `coordination_params` are passed through to `compute_mean_coordination`, for example a custom cutoff.


```python
bond = 'Cu-P'
target = 4

# Optional overrides passed to compute_mean_coordination, for example {'r_cut': 2.6}.
coordination_params = {}

type_A, type_B = (Element(el).Z for el in bond.split('-'))
callable_params = {'type_A': type_A, 'type_B': type_B, **coordination_params}

coordination_fn = get_target_value_fn(
    'compute_mean_coordination',
    **callable_params,
)
```

###         3.1.2 Histograms of average coordination number

The important replacement is the `values` line: instead of `entry.natoms`, each value is `compute_mean_coordination(entry)` for the selected bond.


```python
histograms = []
for label, dataset in datasets.items():
    values = np.array([coordination_fn(entry) for entry in dataset])
    bin_centers, counts = values_2_histo_data(values, name=label, integer_bins=True)
    histograms.append({'label': label, 'bin_centers': bin_centers, 'counts': counts})

fig, ax = plot_multihistogram(
    histograms,
    target=target,
    xlabel=f'average {bond} coordination number',
    title=f'{bond} target = {target}',
    show_gaussian=True,
    max_bincenter=10,
)
fig.savefig(WORK / f'{bond}_{target}_coordination_histogram.png', dpi=200)
```

###         3.1.3 Same calculation through `calculate_values`

This is the compact helper form used internally by the reference notebook. It returns the raw average-coordination values for each dataset label.


```python
coordination_values = calculate_values(
    datasets,
    callable_name='compute_mean_coordination',
    callable_params=callable_params,
    filter_max_el=False,
)

{k: values[:5] for k, values in coordination_values.items()}
```

###         3.1.4 KDE from the same coordination values

Use KDE for a smoother comparison when average coordination values are not limited to integer bins.


```python
fig, ax = plot_multi_kde(
    coordination_values,
    target=target,
    xlabel=f'average {bond} coordination number',
    title=f'{bond} target = {target}',
    max_value=10,
)
fig.savefig(WORK / f'{bond}_{target}_coordination_kde.png', dpi=200)
```

# 4. Pipeline for generation postprocessing

The scenario pipeline makes the usual parse, clean, reference, estimate, hull-filter, and deduplicate steps reproducible.

##     4.1 Scenarios syntax

A scenario is a small DAG. Keep external stages commented or skipped until the required optional tools are configured.

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


```python
from vsbtools.materials_tools.materials_dataset.analysis.scenario_pipeline import ScenarioPipeline

scenario_file = PACKAGE_ROOT / 'materials_tools/materials_dataset/analysis/unittests/scenario_postprocess.yaml'
pipeline = ScenarioPipeline.from_file(scenario_file)
pipeline.ctx.globals['elements'] = sorted(ELEMENTS)
pipeline.ctx.toolkit_options['structure_parser'].pop('batch_metadata_file', None)

generation_root = WORK / GENERATION_NAME

if generation_root.exists():
    rmtree(generation_root)
generation_root.mkdir(parents=True)

zip_file = ZipFile(GENERATION_ZIP)
zip_file.extractall(generation_root)
zip_file.close()

pipeline.ctx.toolkit_options['structure_parser']['root'] = generation_root
partial_outputs = dict(pipeline.run(targets=['parse_raw', 'symmetrize_raw']))

rmtree(generation_root)

partial_outputs.keys()

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
