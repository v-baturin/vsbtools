######## ===== OLD BEHAVIOR ====== ########
import logging
from typing import Dict, Any
from pathlib import Path
from vsbtools.materials_tools.materials_dataset.io.structures_dataset_io import get_batch_metadata, exploded_zip_tree
from vsbtools.materials_tools.ext_software_io.mattergen_tools.parsers import input_parameters_to_dict, fname_friendly_serialize
from vsbtools.materials_tools.materials_dataset.infrastructure import repository
from vsbtools.materials_tools.materials_dataset.analysis import (
    generation_postprocess as gp,
)
from vsbtools.materials_tools.ext_software_io.mattergen_tools.format_tools import build_pairs_from_guidance

here = Path(__file__).parent

# from Pipeline_corrections import repo_path
BATCH_METADATA_FILE = "input_parameters.txt"
INBOX = here / "../../unittests_datasets"
PROCESSED_ROOT = here / "PROCESSED"
structures_source_path = INBOX / "sample_LiCoO_zipped"
DB_CACHE = here / "DB_CACHE"

max_ehull = 0.3  # Both for generated dataset and reference databases


with exploded_zip_tree(structures_source_path) as tmp_inbox:
    batch_meta = get_batch_metadata(tmp_inbox, BATCH_METADATA_FILE)
    param_dict = input_parameters_to_dict(raw=batch_meta)
    system = '-'.join(sorted(param_dict['properties_to_condition_on']['chemical_system'].split('-')))
    repo_path = PROCESSED_ROOT / system / fname_friendly_serialize(param_dict)


    repo = repository.DatasetRepo(root=repo_path)
    elements_set = set(system.split('-'))

    processed_stages = dict()
    graph = repo.build_graph()
    for node in repo.list_nodes():
        ds, meta = repo.load_node(node)
        processed_stages[gp.PostprocessStage(int(ds.metadata["pipeline_stage"]))] = ds

    polling_db_options: Dict[str, Any] = dict(do_ehull_filtering=True, do_deduplication=True, max_ehull=max_ehull)
    db_cache_name = (f"{'_'.join(sorted(elements_set))}_{int(polling_db_options['do_deduplication'])}_"
                     f"{int(polling_db_options['do_deduplication'])}_{polling_db_options['max_ehull']}")
    polling_db_options["cache_base_path"] = DB_CACHE / db_cache_name

    pipeline = gp.PPPipeline(source_path=tmp_inbox, processed_stages=processed_stages, elements=elements_set,
                             toolkit_options={"structure_parser": {"batch_metadata_file": BATCH_METADATA_FILE},
                                              "symmetry": {"a_sym_prec": 1e-3, "e_sym_prec": 1e-3},
                                              "similarity": {"tol_FP": 0.04}},
                             stages_options={gp.PostprocessStage.poll_db: polling_db_options,
                                             gp.PostprocessStage.filter_hull: {'max_ehull': max_ehull}},
                             root_source_name='Mattergen')

    history = pipeline.run()

    for stage, ds in history:
        repo.commit_node(ds, prefix=str(stage.value))
        processed_stages[stage] = ds