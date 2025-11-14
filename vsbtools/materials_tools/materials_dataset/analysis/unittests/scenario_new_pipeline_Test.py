from pathlib import Path


from typing import Dict, Any
import logging
from pathlib import Path
from vsbtools.materials_tools.materials_dataset.io.structures_dataset_io import get_batch_metadata, exploded_zip_tree
from vsbtools.materials_tools.ext_software_io.mattergen_tools.parsers import input_parameters_to_dict, fname_friendly_serialize
from vsbtools.materials_tools.materials_dataset.infrastructure import repository
from vsbtools.materials_tools.materials_dataset.analysis.scenario_pipeline import ScenarioPipeline

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


here = Path(__file__).parent

BATCH_METADATA_FILE = "input_parameters.txt"
INBOX = here / "../../unittests_datasets"
PROCESSED_ROOT = here / "PROCESSED_NEW"
structures_source_path = INBOX / "sample_LiCoO_zipped"
SCENARIO_FILE = here / "scenario_postprocess.yaml"      # scenario definition lives here

with exploded_zip_tree(structures_source_path) as tmp_inbox:
    # Read batch metadata and derive chemical system
    batch_meta = get_batch_metadata(tmp_inbox, BATCH_METADATA_FILE)
    param_dict = input_parameters_to_dict(raw=batch_meta)
    system = "-".join(sorted(param_dict["properties_to_condition_on"]["chemical_system"].split("-")))
    repo_path = PROCESSED_ROOT / system / fname_friendly_serialize(param_dict)

    repo = repository.DatasetRepo(root=repo_path)
    elements = sorted(set(system.split("-")))

    # Build pipeline from scenario
    pipeline = ScenarioPipeline.from_file(SCENARIO_FILE)
    ctx = pipeline.ctx

    # Minimal batch-specific context
    ctx.globals["elements"] = elements
    ctx.toolkit_options["structure_parser"]["root"] = tmp_inbox

    # Reuse already processed stages from repo (legacy numeric metadata)
    for node in repo.list_nodes():
        ds, meta = repo.load_node(node)
        stage_name = pipeline.scenario.resolve_stage_name_from_metadata(ds)
        ctx.outputs[stage_name] = ds

    # Run only missing stages (pipeline.run() should skip those already in ctx.outputs)
    for stage_name, ds in pipeline.run():
        idx = pipeline.scenario.legacy_name_to_index.get(stage_name)
        prefix = str(idx) if idx is not None else stage_name
        repo.commit_node(ds, prefix=prefix)
