import unittest
from pathlib import Path

import numpy as np
from pymatgen.core import Structure

from ....crystal_entry import CrystalEntry
from ....io.zip_handling import exploded_zip_tree
from ..mattergen_bridge import get_loss_fn, get_target_value_fn, clear_globals

FIXTURES_ROOT = Path(__file__).resolve().parents[3] / "unittests_datasets"
ZIPPED_PROCESSED_PATH = (
    FIXTURES_ROOT / "MG_postprocess_pipelines" / "PROCESSED"
)


class MattergenBridge_Test(unittest.TestCase):
    def setUp(self) -> None:
        self._fixtures_context = exploded_zip_tree(ZIPPED_PROCESSED_PATH)
        self.processed_path = self._fixtures_context.__enter__()
        self.poscars_path = (
            self.processed_path / "B-Fe-Nd"
            / "B-Fe-Nd__guidance_environment_mode_huber_B-Fe_3__diffusion_loss_weight_0.5-0.5-True__algo_0"
            / "2_x204f1f3d5ffe2c98" / "POSCARS"
        )

    def tearDown(self) -> None:
        self._fixtures_context.__exit__(None, None, None)

    def test_environment_loss_target_stability(self):
        entry1 = CrystalEntry(id="agm003592845", structure=Structure.from_file(self.poscars_path / "agm003592845POSCAR"))
        entry2 = CrystalEntry(id='mp-650968', structure=Structure.from_file(self.poscars_path / "mp-650968POSCAR"))
        mean_cn_fn = get_target_value_fn(
            "compute_mean_coordination", force_gpu=0, type_A=5, type_B=26
        )
        print(f"CN{entry1.id}(Co-O): {float(mean_cn_fn(entry1))}")
        # print(f"CN{entry2.id}(Co-O): {float(mean_cn_fn(entry2))}")
        clear_globals()

        loss_fn_target6 = get_loss_fn('environment', force_gpu=0, target={'B-Fe': 6, 'mode': 'l1'})
        print(f"loss6({entry1.id}) = {loss_fn_target6(entry1)}")
        # print(f"loss6({entry2.id}) = {loss_fn_target6(entry2)}")

        loss_fn_target3 = get_loss_fn('environment', force_gpu=0, target={'B-Fe': 3, 'mode': 'l1'})
        print(f"loss3({entry1.id}) = {loss_fn_target3(entry1)}")
        print(f"loss6({entry1.id}) = {loss_fn_target6(entry1)}")
        #
        # value_target3_second = float(loss_fn_target3(entry))
        #
        # if np.isclose(value_target3_first, value_target6):
        #     self.skipTest(
        #         f"Target=3 and target=4 produced identical losses = {value_target3_first}"
        #     )
        #
        # self.assertAlmostEqual(value_target3_second, value_target3_first, places=10)
