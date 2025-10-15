import warnings
from typing import Dict, Any, Set, Generator, Tuple, MutableMapping
from enum import Enum
from pathlib import Path
from dataclasses import dataclass, field
from collections import defaultdict
from ..crystal_dataset import CrystalDataset
from ..io.structures_dataset_io import StructureDatasetIO
from .symmetry_tools import SymmetryToolkit
from ..scripts.poll_databases import poll_databases
from .similarity_tools import SimilarityTools
from ..io.uspex_bridge import USPEXBridge
from ..energy_estimation import nn_estimator, mattersim_bridge
from ..analysis.phase_diagram_tools import PhaseDiagramTools

TOOLKIT_DICT = {"structure_parser":StructureDatasetIO, "symmetry": SymmetryToolkit,
                "similarity": SimilarityTools, "uspex": USPEXBridge, "phase_diag": PhaseDiagramTools,
                "estimator": nn_estimator.NNEstimator}

KNOWN_MODELS = {"mattersim": mattersim_bridge}

MAX_EHULL_PA = 0.1

class PostprocessStage(Enum):
    parse_raw = 0
    symmetrize_raw = 1
    poll_db = 2
    augment_raw_by_db = 3
    estimate = 4
    filter_hull = 5
    deduplicate = 6
    postprocess_dft = 7


@dataclass
class PPPipeline:
    source_path: Path | None = None
    elements: Set[str] | None = None
    root_source_name: str = 'NA',
    processed_stages: Dict[PostprocessStage, CrystalDataset] = field(default_factory=dict)
    stages_options: Dict[PostprocessStage, Dict[str, Any]] = field(default_factory=dict)
    toolkits: dict = field(default_factory=dict)
    toolkit_options: MutableMapping[str, Dict[str, Any]] = field(default_factory=lambda : defaultdict(dict))

    def __post_init__(self):
        USPEXBridge.uspex_entry_from_de.cache_clear()
        if not isinstance(self.toolkit_options, defaultdict):
            self.toolkit_options = defaultdict(dict, self.toolkit_options)
        self.toolkit_options["structure_parser"].update({"root": self.source_path,
                                                         "source_name": self.root_source_name})
        self.toolkit_options["uspex"].update({"elements": self.elements})
        self.stages_options = defaultdict(dict, {PostprocessStage(k): v for k, v in self.stages_options.items()})


    def get_tool(self, toolkit_name):
        if toolkit_name not in self.toolkits:
            if toolkit_name == "similarity":
                ub = self.get_tool("uspex")
                self.toolkit_options["similarity"]["dist_fn"] = ub.fp_dist
            elif toolkit_name == "estimator":
                for label, module in KNOWN_MODELS.items():
                    TOOLKIT_DICT[toolkit_name].register_model(label, module)
            self.toolkits[toolkit_name] = TOOLKIT_DICT[toolkit_name](**self.toolkit_options.get(toolkit_name, dict()))
        return self.toolkits[toolkit_name]


    def run(self, to_process: list | None = None,
            skip_list: list | None = None
            ) -> Generator[Tuple[PostprocessStage,CrystalDataset], None, None]:
        
        skip_list = [] if skip_list is None else skip_list
        skip_list = [PostprocessStage(i) for i in skip_list]
        
        to_process = to_process if to_process is not None else \
            [stage for stage in PostprocessStage if stage not in self.processed_stages and stage not in skip_list]
        
        for stage in to_process:
            print(f"Now working on stage {stage.name} ...")

            if stage is PostprocessStage.parse_raw:
                assert self.source_path and self.source_path.exists()
                assert self.elements is not None
                raw_parse = self.get_tool("structure_parser").load_from_directory(elements=self.elements)
                assert len(raw_parse) > 0, f"Empty dataset in {self.source_path}"
                print(f"loaded {len(raw_parse)} structures")
                self.processed_stages[stage] = raw_parse

            if stage is PostprocessStage.symmetrize_raw:
                symtool = self.get_tool("symmetry")
                self.processed_stages[stage] = (
                    symtool.get_symmetrized_dataset(self.processed_stages[PostprocessStage.parse_raw]))

            if stage is PostprocessStage.poll_db:
                elements = self.processed_stages[PostprocessStage.parse_raw].elements
                if "max_ehull" not in self.stages_options[stage]:
                    self.stages_options[stage]["max_ehull"] = MAX_EHULL_PA
                self.processed_stages[stage] = poll_databases(elements, **self.stages_options[stage])
                if "estimate_energies" in self.stages_options[stage] and self.stages_options[stage]["estimate_energies"]:
                    estimator = self.get_tool("estimator")
                    estimated = estimator.estimate_dataset_energies(self.processed_stages[stage])
                    estimated.parent_ids = []
                    estimated.metadata["message"] = (self.processed_stages[stage].metadata["message"].strip() + '. '
                                                     + estimated.metadata["message"].strip())
                    self.processed_stages[stage] = estimated


            if stage is PostprocessStage.augment_raw_by_db:
                similarity_tk: SimilarityTools = self.get_tool("similarity")
                augmentation = similarity_tk.get_unseen_in_ref(self.processed_stages[PostprocessStage.symmetrize_raw],
                                                               ref_ds=self.processed_stages[PostprocessStage.poll_db])
                print(augmentation.metadata["message"])
                labeled_entries = []
                for e in self.processed_stages[PostprocessStage.poll_db]:
                    if str(e.id) in augmentation.metadata["reproduced"]:
                        new_metadata = {**e.metadata, **{"reproduced": True}}
                        labeled_entries.append(e.copy_with(**{"metadata": new_metadata}))
                    else:
                        labeled_entries.append(e)
                db_parameters = self.processed_stages[PostprocessStage.poll_db].__dict__.copy()
                meta = db_parameters.pop("metadata")
                db_parameters["base_path"] = db_parameters.pop("_base_path")
                to_drop = [k for k in db_parameters if k.startswith('_')]
                for k in to_drop: db_parameters.pop(k)
                labeled_db = CrystalDataset(labeled_entries, **db_parameters)
                labeled_db.metadata = meta
                self.processed_stages[stage] = labeled_db.merge(augmentation)
                self.processed_stages[stage].metadata["message"] = "Merged DB data into generated set"
                self.processed_stages[stage].parent_ids = \
                    [self.processed_stages[s].dataset_id for s in
                     [PostprocessStage.symmetrize_raw, PostprocessStage.poll_db]]
                self.processed_stages[stage].metadata["reproducibility"] = augmentation.metadata["reproducibility"]

            if stage is PostprocessStage.estimate:
                estimator = self.get_tool("estimator")
                self.processed_stages[stage] = estimator.estimate_dataset_energies(
                    self.processed_stages[PostprocessStage.augment_raw_by_db])

            if stage is PostprocessStage.filter_hull:
                if "base_stage" in self.stages_options[stage]:
                    base_ds = self.processed_stages[PostprocessStage(self.stages_options[stage]["base_stage"])]
                else:
                    base_ds = self.processed_stages[PostprocessStage.estimate]
                self.toolkit_options["phase_diag"] = {"dataset": base_ds}
                pd_tk: PhaseDiagramTools = self.get_tool("phase_diag")
                max_ehull = self.stages_options[stage].get("max_ehull", MAX_EHULL_PA)
                filtered = self.processed_stages[PostprocessStage.estimate].filter(
                    lambda entry: pd_tk.height_above_hull_pa(entry) < max_ehull)
                filtered.metadata["message"] = f"Selected entries with e_hull < {max_ehull:.3f}"
                self.processed_stages[stage] = filtered

            if stage is PostprocessStage.deduplicate:
                similarity_tk: SimilarityTools = self.get_tool("similarity")
                self.processed_stages[stage], _, _ =\
                    similarity_tk.deduplicate(self.processed_stages[PostprocessStage.filter_hull],
                                              check_clusters_file=True, check_dist_matrix_file=True)
                self.processed_stages[stage].parent_ids = [self.processed_stages[PostprocessStage.filter_hull].dataset_id]
                self.processed_stages[stage].metadata["message"] = \
                    f"Deduplicated version of {self.processed_stages[PostprocessStage.filter_hull].dataset_id}"

            if stage in self.processed_stages:
                print(f"Finished {stage.name} stage")
                self.processed_stages[stage].metadata["pipeline_stage"] = stage.value
                yield stage, self.processed_stages[stage]

            else:
                warnings.warn(f"Stage '{stage.name}' not yet implemented")
                break







