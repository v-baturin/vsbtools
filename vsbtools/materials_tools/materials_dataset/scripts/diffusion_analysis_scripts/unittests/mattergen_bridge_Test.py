import unittest
from pathlib import Path

import numpy as np
from pymatgen.core import Structure

from ....crystal_entry import CrystalEntry
from ..mattergen_bridge import get_loss_fn, get_target_value_fn, clear_globals

POSCARS_PATH = Path(
    "/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/PROCESSED/"
    "B-Fe-Nd/B-Fe-Nd__guidance_environment_mode_huber_B-Fe_6__diffusion_loss_weight_0.5-0.5-True__algo_0"
    "/add_ref_all_311a5df7ceb60f7e/POSCARS/"
)


class MattergenBridge_Test(unittest.TestCase):
    def test_environment_loss_target_stability(self):
        entry1 = CrystalEntry(id="agm003592845", structure=Structure.from_file(POSCARS_PATH / "agm003592845POSCAR"))
        entry2 = CrystalEntry(id='mp-650968', structure=Structure.from_file(POSCARS_PATH / "mp-650968POSCAR"))
        mean_cn_fn = get_target_value_fn(
            "compute_mean_coordination", type_A=5, type_B=26
        )
        print(f"CN{entry1.id}(Co-O): {float(mean_cn_fn(entry1))}")
        # print(f"CN{entry2.id}(Co-O): {float(mean_cn_fn(entry2))}")
        clear_globals()

        loss_fn_target6 = get_loss_fn('environment', target={'B-Fe': 6, 'mode': 'l1'})
        print(f"loss6({entry1.id}) = {loss_fn_target6(entry1)}")
        # print(f"loss6({entry2.id}) = {loss_fn_target6(entry2)}")

        loss_fn_target3 = get_loss_fn('environment', target={'B-Fe': 3, 'mode': 'l1'})
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
