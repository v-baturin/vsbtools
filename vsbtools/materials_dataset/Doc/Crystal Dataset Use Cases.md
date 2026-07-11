Examples assume they are run from the repository root `/home/vsbat/my_git_projects/vsbtools`. The raw-generation fixtures live under `Examples/raw_generations`. Each first-level zip in `RAW_GENERATION_DIR` is one MatterGen generation: `Cu-Si-P_nonguided.zip` is single-run, while the guided archives contain multiple `run_*` directories. For examples, unzip one first-level archive into the tutorial work directory, call the regular `StructureDatasetIO(...).load_from_directory(...)` on the extracted directory, then delete the extracted files when the loaded dataset is no longer using them.


```python
from pathlib import Path

REPO_ROOT = Path('../../../..')
RAW_GENERATIONS = REPO_ROOT / 'vsbtools/materials_dataset/Examples/raw_generations'
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

from vsbtools.materials_dataset.io.structures_dataset_io import StructureDatasetIO


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
        abc : 3.6638589564547415 3.654074180905444 5.636622428894043
     angles : 62.445883123840034 62.31136479787251 59.47841537074956
     volume : 54.9640757737494
          A : 3.244295120239258 0.0 1.70247220993042
          B : 1.2087509941275552 3.0056584790443464 1.690324306488037
          C : 0.0 0.0 5.636622428894043
        pbc : True True True
    PeriodicSite: Si (0.89, 0.7783, 4.129) [0.1778, 0.2589, 0.6011]
    PeriodicSite: P (3.47, 2.524, 6.231) [0.7567, 0.8399, 0.6249]
    PeriodicSite: Cu (0.4084, 0.4558, 3.671) [0.06937, 0.1517, 0.5849], energy=None, formula='Si1 Cu1 P1', metadata={'source': 'MatterGen', 'file': 'generated_crystals.extxyz'})



###         1.1.2 Poll Databases

Database polling is optional and may need API keys, network access, or a local cache.

#### 1.1.2.1 Polling locally deployed databases

Requires local Alexandria json files and OQMD database


```python
from vsbtools.materials_dataset.io.preset_loaders import load_from_oqmd

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
from vsbtools.materials_dataset.io.preset_loaders import load_from_alexandria

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
from vsbtools.materials_dataset.io.preset_loaders import load_from_materials_project

# Requires an MP API key/configured environment.
reference_ds = load_from_materials_project(
    elements={'Cu', 'Si', 'P'},
    message='Cu-Si-P reference phases from Materials Project',
)
print(len(reference_ds))

```


    Retrieving ThermoDoc documents:   0%|          | 0/8 [00:00<?, ?it/s]



    Retrieving ThermoDoc documents:   0%|          | 0/14 [00:00<?, ?it/s]



    Retrieving ThermoDoc documents:   0%|          | 0/42 [00:00<?, ?it/s]



    Retrieving ThermoDoc documents:   0%|          | 0/7 [00:00<?, ?it/s]



    Retrieving ThermoDoc documents:   0%|          | 0/5 [00:00<?, ?it/s]



    Retrieving ThermoDoc documents:   0%|          | 0/5 [00:00<?, ?it/s]



    Retrieving ThermoDoc documents:   0%|          | 0/6 [00:00<?, ?it/s]



    Retrieving SummaryDoc documents:   0%|          | 0/87 [00:00<?, ?it/s]


    Data saved to /home/vsbat/.cache/vsbtools/DB_caches/materials_project_Cu-P-Si_none
    87


**Polling using OPTIMADE interface**


```python
# This request will poll and deduplicate info sequentiall from Alexandria, Materials Project and OQMD.
from vsbtools.materials_dataset.io.preset_loaders import load_from_optimade

providers = ("alexandria", "mp", "oqmd")
reference_ds = load_from_optimade(
    elements={'Si', 'P'},
    message=f'Cu-Si-P from {", ".join(providers)}',
    providers=providers,
    timeout=1500, force_refresh=True)

print(len(reference_ds))
```

    OPTIMADE provider 1/3 alexandria: finished, kept 296, total accepted 296
    Data saved to /home/vsbat/.cache/vsbtools/DB_caches/optimade_alexandria_P-Si_9f9896d263
    OPTIMADE provider 2/3 materials_project: finished, kept 61, total accepted 61
    Data saved to /home/vsbat/.cache/vsbtools/DB_caches/optimade_materials_project_P-Si_58295233ed
    OPTIMADE provider 3/3 oqmd: finished, kept 137, total accepted 137
    Data saved to /home/vsbat/.cache/vsbtools/DB_caches/optimade_oqmd_P-Si_1c54297fdd
    OPTIMADE compiling with deduplication:  3/3 oqmd: merging entry 137/137
    Data saved to /home/vsbat/.cache/vsbtools/DB_caches/optimade_P-Si_62d573e7be
    422


##     1.2 Clean the Dataset

Cleaning usually combines geometry checks, duplicate removal, and optional symmetry normalization.

###         1.2.1 Remove structures with too short bonds

Filter entries by the package geometry validator.


```python
from vsbtools.materials_dataset.geom_utils.structure_checks import check_min_dist_pmg

clean_ds = raw_ds.filter(lambda entry: check_min_dist_pmg(entry.structure)[0])
print(f'{len(clean_ds)} / {len(raw_ds)} structures passed the minimum-distance check')

```

    561 / 708 structures passed the minimum-distance check


###         1.2.2 Remove Duplicates

The default in-package duplicate workflow uses `SimilarityTools` with a DScribe-backed fingerprint distance.


```python
from vsbtools.materials_dataset.analysis.similarity_tools import SimilarityTools
from vsbtools.materials_dataset.analysis.structural_distance.dscribe_bridge import DScribeBridge

bridge = DScribeBridge(elements={'Cu', 'Si', 'P'}, tol_FP=0.04)
similarity = SimilarityTools(bridge.fp_dist, bridge.tol_FP)

dedup_ds, clusters, best_idx = similarity.deduplicate(
    clean_ds,
    enforce_compositions_separation=True,
)
print(f'{len(dedup_ds)} unique structures in {len(clusters)} clusters')

```

    i = 560, j = 560
    dist matrix saved
    Separating 536 clusters by labels list
    Separated into 541 clusters to have same labels per clusters
    clustering saved in /home/vsbat/my_git_projects/vsbtools/vsbtools/materials_dataset/Doc/x139704a0b9a541d2_clusters.pkl
    Processing 541 clusters
    541 good ones. Writing
    541 unique structures in 541 clusters


###         1.2.3 Symmetrize the Entries

Symmetrization returns a new dataset with standardized structures.


```python
from vsbtools.materials_dataset.analysis.symmetry_tools import SymmetryToolkit

symmetry = SymmetryToolkit(a_sym_prec=1e-3,  # precision of analysis
                           e_sym_prec=1e-2)  # precision of enforcing a symmetry
sym_ds = symmetry.get_symmetrized_dataset(clean_ds)

changed_symmetry_count = sum(
    symmetry.sym_group_no(original) != symmetry.sym_group_no(symmetrized)
    for original, symmetrized in zip(clean_ds, sym_ds)
)
print(f'{changed_symmetry_count} / {len(sym_ds)} structures changed symmetry after symmetrization')

```

    156 / 561 structures changed symmetry after symmetrization


##     1.3 Estimate the energy

Energy-dependent analysis expects each `CrystalEntry.energy` to be populated, either from a database, DFT, or an estimator.

###         1.3.1 Construct the convex hull and calculate heights*

Once energies are present, build a phase diagram from the reference dataset and evaluate generated entries against it.

\*special case of **2.1 Heights of Dataset 1 above the convex hull from Dataset 2 (generated wrt reference)**.


```python
from vsbtools.materials_dataset.analysis.phase_diagram_tools import PhaseDiagramTools

# `reference_ds` and `generated_ds` must contain energies in eV per structure.
pd_tools = PhaseDiagramTools(reference_ds)
heights_pa = [pd_tools.height_above_hull_pa(entry) for entry in reference_ds]
print(heights_pa[:5])

```

    [0.037360772, 0.5705513, 0.37145662, 0.10703417, 0.0828698]


##     1.4 Various I/O operations

Use the converter helpers when you need a tabular view for reporting or quick inspection.

###         1.4.1 Dataset <> Pandas DataFrame

DataFrames are convenient for filtering by scalar metadata and joining with external tables.


```python
from vsbtools.materials_dataset.converters import ds2df, df2ds

frame = ds2df(clean_ds)
print(frame.head())

roundtrip_ds = df2ds(frame)
print(len(roundtrip_ds))

```

    Dropped columns with only None values: energy
                          id                                          structure  \
    0  Cu-Si-P_guided_CuP3_1  [[3.25031517 4.52660315 6.92366093] Cu, [3.029...   
    1  Cu-Si-P_guided_CuP3_3  [[2.37990909 0.65891346 4.91428714] Si, [1.820...   
    2  Cu-Si-P_guided_CuP3_4  [[ 0.82289662  4.28663968 -0.04097639] Si, [3....   
    3  Cu-Si-P_guided_CuP3_5  [[ 1.76236093  3.55165674 -0.78063469] P, [1.6...   
    4  Cu-Si-P_guided_CuP3_6  [[2.49272326 1.87171417 6.06793435] Cu, [3.936...   
    
                                                metadata  
    0  {'source': 'MatterGen', 'file': 'generated_cry...  
    1  {'source': 'MatterGen', 'file': 'generated_cry...  
    2  {'source': 'MatterGen', 'file': 'generated_cry...  
    3  {'source': 'MatterGen', 'file': 'generated_cry...  
    4  {'source': 'MatterGen', 'file': 'generated_cry...  
    561


##     1.5 Building reports

Summary helpers collect common scalar fields into a compact table.

###         1.5.1 Summary table (text and csv)

Write a human-readable table and a CSV from the same summary frame.


```python
from vsbtools.materials_dataset.analysis.summary import collect_summary_df, print_pretty_df

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

    
    summary.txt
    +-------------------------+-------------+--------+
    | id                      | composition | energy |
    +-------------------------+-------------+--------+
    | Cu-Si-P_guided_CuP3_1   | Cu3 P8 Si1  | None   |
    | Cu-Si-P_guided_CuP3_3   | Si2 Cu1 P1  | None   |
    | Cu-Si-P_guided_CuP3_4   | Si3 P4 Cu1  | None   |
    | Cu-Si-P_guided_CuP3_5   | P3 Si2 Cu2  | None   |
    | Cu-Si-P_guided_CuP3_6   | Cu2 P5 Si3  | None   |
    | Cu-Si-P_guided_CuP3_7   | Si2 P5 Cu1  | None   |
    | Cu-Si-P_guided_CuP3_8   | Si4 P3 Cu1  | None   |
    | Cu-Si-P_guided_CuP3_9   | P7 Cu2 Si1  | None   |
    | Cu-Si-P_guided_CuP3_10  | Cu1 Si1 P2  | None   |
    
    summary.csv
    ,id,composition,energy
    0,Cu-Si-P_guided_CuP3_1,Cu3 P8 Si1,
    1,Cu-Si-P_guided_CuP3_3,Si2 Cu1 P1,
    2,Cu-Si-P_guided_CuP3_4,Si3 P4 Cu1,
    3,Cu-Si-P_guided_CuP3_5,P3 Si2 Cu2,
    4,Cu-Si-P_guided_CuP3_6,Cu2 P5 Si3,
    5,Cu-Si-P_guided_CuP3_7,Si2 P5 Cu1,
    6,Cu-Si-P_guided_CuP3_8,Si4 P3 Cu1,
    7,Cu-Si-P_guided_CuP3_9,P7 Cu2 Si1,
    8,Cu-Si-P_guided_CuP3_10,Cu1 Si1 P2,
    9,Cu-Si-P_guided_CuP3_12,Si8 Cu1 P1,
    10,Cu-Si-P_guided_CuP3_13,Cu7 P3,


##     1.2 Read and write Datasets on disk

The YAML/CSV/POSCAR writer creates a portable manifest plus structure files.
Dataset is read via its `manifest.yaml` file


```python
from vsbtools.materials_dataset.io.yaml_csv_poscars import read, write

write(
    clean_ds,
    enforce_base_path=WORK / 'clean_dataset',
    comment='Minimal cleaned raw-generation fixture',
)
manifest = WORK / 'clean_dataset' / 'manifest.yaml'
restored_ds = read(manifest)
print(len(restored_ds), restored_ds.base_path)

```

    Data saved to tmp_crystal_dataset_examples/clean_dataset
    561 tmp_crystal_dataset_examples/clean_dataset


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
from vsbtools.materials_dataset.analysis.pareto_fronts import plot_pareto_fronts

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

For raw fixtures, classify generations by first-level archive name. To compare generations, unzip each archive into the tutorial work directory, call the same directory loader on each extracted generation, then delete the extracted files.


```python
generation_archives = sorted(RAW_GENERATION_DIR.glob('*.zip'))
guided_archives = [path for path in generation_archives if '_guided_' in path.stem]
nonguided_archives = [path for path in generation_archives if '_nonguided' in path.stem]
print('guided:', [path.name for path in guided_archives])
print('non-guided:', [path.name for path in nonguided_archives])

datasets_by_generation = {}
for archive in generation_archives:
    generation_root = WORK / archive.stem

    if generation_root.exists():
        rmtree(generation_root)
    generation_root.mkdir(parents=True)

    zip_file = ZipFile(archive)
    zip_file.extractall(generation_root)
    zip_file.close()

    datasets_by_generation[archive.stem] = StructureDatasetIO(
        generation_root,
        id_prefix=f'{archive.stem}_',
        source_name='MatterGen',
    ).load_from_directory(elements=ELEMENTS)

    rmtree(generation_root)

{k: len(v) for k, v in datasets_by_generation.items()}

```

###         3.1.2 Histograms of a descriptor

This example plots the number of atoms per structure for each dataset.


```python
import numpy as np
from vsbtools.materials_dataset.analysis.guidance_statistics import values_2_histo_data, plot_multihistogram

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
from vsbtools.materials_dataset.analysis.guidance_statistics import plot_multi_kde

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
from vsbtools.materials_dataset.analysis.scenario_pipeline import ScenarioPipeline

scenario_file = REPO_ROOT / 'vsbtools/materials_dataset/analysis/unittests/scenario_postprocess.yaml'
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
