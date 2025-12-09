import unittest
from pathlib import Path
from ...analysis.scenario_pipeline import Scenario
from ...scripts.diffusion_analysis_scripts.tools_for_histograms import \
    get_guidance_generation_dirs
from ...io.yaml_csv_poscars import read

from ..build_tables import build_energy_vs_property_table


class yaml_csv_poscars_Test(unittest.TestCase):

    def setUp(self):
        self.environment_repo = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/"
                                     "PROCESSED/O-Si/"
                                     "O-Si__guidance_environment_mode_huber_Si-O_6__diffusion_loss_weight_1-1-True__algo_0")

        self.volume_pa_repo = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/"
                                   "PROCESSED/B/B__guidance_volume_pa_6.8__diffusion_loss_weight_1-1-True__algo_0")

        self.environment_repo_multi = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/"
            "MG_postprocess_pipelines/PROCESSED/Ca-Cu-P-Si/"
            "Ca-Cu-P-Si__guidance_environment_mode_huber_Cu-P_4-2.2_Cu-Cu_0-2.9__diffusion_loss_weight_1-1-False__algo_None",
            )

    def test_environment_repo(self):
        ds_dict = dict()
        for manifest_yaml in self.environment_repo_multi.rglob("manifest.yaml"):
            ds = read(manifest_yaml)
            ds_dict[ds.metadata["pipeline_stage"]] = ds
        print(ds_dict.keys())
        build_energy_vs_property_table(ds_dict)