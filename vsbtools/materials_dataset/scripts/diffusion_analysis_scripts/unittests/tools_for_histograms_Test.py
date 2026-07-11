from typing import Callable, Dict
import os
import unittest
import numpy as np
from pathlib import Path
from pymatgen.core import Structure, Element
from ....crystal_dataset import CrystalDataset
from ....crystal_entry import CrystalEntry
from ....io.zip_handling import exploded_zip_tree
from ..guidance_stats import (get_environment_gen_dirs, get_volume_pa_gen_dirs, collect_stage_dataset_dict,
                              histo_data_collection, get_target_value_fn, plot_multihistogram, calculate_values,
                              plot_multi_kde)
from matplotlib import pyplot as plt

FIXTURES_ROOT = Path(__file__).resolve().parents[3] / "unittests_datasets"
ZIPPED_PROCESSED_PATH = FIXTURES_ROOT / "MG_postprocess_pipelines" / "PROCESSED"
RUN_INTEGRATION_TESTS = os.getenv("VSBTOOLS_RUN_INTEGRATION_TESTS") == "1"
integration_test = unittest.skipUnless(
    RUN_INTEGRATION_TESTS,
    "Set VSBTOOLS_RUN_INTEGRATION_TESTS=1 to run integration/plot tests",
)


def count_entries_around_target(ds: CrystalDataset, function: Callable, target_value: float, half_width: float = 0.5, profile: str = "rect", **kwargs):
    if profile.startswith("rect"):
        filt_predicate = lambda e: target_value - half_width <= function(e) <= target_value + half_width and len(e.composition) == len(ds.elements)
    else:
        raise RuntimeError(f"profile: {profile} not implemented")
    return len([e for e in ds if filt_predicate(e)])


def print_ds_dict_guidance_resume(ds_dict: Dict, function: Callable, target: float, half_width: float, **kwargs):
    fractions_in_target_dict = dict()
    for k, ds in ds_dict.items():
        fractions_in_target_dict[k] = count_entries_around_target(ds, function, target, half_width, **kwargs) / len(ds)
        print(f"{fractions_in_target_dict[k]:1%} of {k} in ({target:.2f}{chr(int("00B1", 16))}{half_width:.2f})")
    return fractions_in_target_dict

class hist_tools_Test(unittest.TestCase):

    def setUp(self) -> None:
        self._fixtures_context = exploded_zip_tree(ZIPPED_PROCESSED_PATH)
        self.processed_path = self._fixtures_context.__enter__()
        sio2_poscars = (
            self.processed_path / "O-Si"
            / "O-Si__guidance_environment_mode_huber_Si-O_6__diffusion_loss_weight_1-1-True__algo_0_for_histograms_stable"
            / "2_x69650df66373f8bf" / "POSCARS"
        )
        self.sio2_stich = sio2_poscars / "agm002170463POSCAR"
        self.sio2_quartz = sio2_poscars / "agm002228342POSCAR"
        self.system_sio = "Si-O"
        self.system_licoo = "Li-Co-O"
        self.system_cupsi = "Cu-P-Si"
        self.system_bfend = "B-Fe-Nd"
        self.system_b = "B"

    def tearDown(self) -> None:
        self._fixtures_context.__exit__(None, None, None)

    def _assert_histogram_inputs(self, dirs, ds_dict, hdc):
        self.assertGreater(len(dirs), 0)
        self.assertGreater(len(ds_dict), 0)
        self.assertEqual({item["label"] for item in hdc}, set(ds_dict))

        non_empty = [item for item in hdc if item["counts"] is not None]
        self.assertGreater(len(non_empty), 0)
        for item in non_empty:
            self.assertIsNotNone(item["bin_centers"])
            self.assertTrue(np.all(np.isfinite(item["counts"])))
            self.assertAlmostEqual(float(np.sum(item["counts"])), 1.0, places=6)

    def test_get_average_cn_gen_dirs(self):
        dirs = get_environment_gen_dirs(self.processed_path, self.system_sio, guidance_name='environment', bond='Si-O', target=4)
        self.assertEqual(len(dirs), 2)

    def test_get_entry_fn(self):
        fn = get_target_value_fn(fn_name="compute_mean_coordination", force_gpu=0, type_A=14, type_B=8)
        entry = CrystalEntry(id='0', structure=Structure.from_file(self.sio2_quartz))
        self.assertAlmostEqual(float(fn(entry)), 3.9771389961242676, places=6)

        fn = get_target_value_fn(fn_name="_soft_neighbor_counts_per_A_single", force_gpu=0, type_A=14, type_B=8)
        cns = np.asarray(fn(entry), dtype=float)
        self.assertEqual(cns.shape, (2,))
        np.testing.assert_allclose(cns, [3.9771419, 3.9771361], rtol=1e-6, atol=1e-6)

    def test_collect_stage_dataset_dict(self):
        dirs = get_environment_gen_dirs(self.processed_path, self.system_sio, guidance_name='environment', bond='Si-O', target=6)
        self.assertGreater(len(dirs), 0)
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)
        self.assertIn("reference", ds_dict)
        self.assertIn("Non-guided", ds_dict)
        self.assertTrue(any("environment" in key for key in ds_dict))
        hdc = histo_data_collection(ds_dict, callable_name='compute_mean_coordination', callable_params={"type_A": 14,
                                                                                                         "type_B": 8},
                                    max_bincenter=10)
        self.assertEqual({item["label"] for item in hdc}, set(ds_dict))
        for item in hdc:
            self.assertIsNotNone(item["bin_centers"])
            self.assertIsNotNone(item["counts"])
            self.assertTrue(np.all(np.isfinite(item["counts"])))
            self.assertAlmostEqual(float(np.sum(item["counts"])), 1.0, places=6)
        fig, _ax = plot_multihistogram(multidata=hdc, target=6, max_bincenter=10)
        self.assertGreater(len(fig.axes), 0)
        plt.close(fig)

    @integration_test
    def test_licoo(self):
        dirs = get_environment_gen_dirs(self.processed_path, self.system_licoo, guidance_name='environment', bond='Co-O', target=4)
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='compute_mean_coordination', callable_params={"type_A": 27,
                                                                                                         "type_B": 8},
                                    max_bincenter=10)
        self._assert_histogram_inputs(dirs, ds_dict, hdc)
        fig, _ax = plot_multihistogram(multidata=hdc, target=4, max_bincenter=10, other_cmap='RdBu', simplified_legend=True)
        self.assertGreater(len(fig.axes), 0)
        plt.close(fig)


    @integration_test
    def test_cupsi(self):
        dirs = get_environment_gen_dirs(self.processed_path, self.system_cupsi, guidance_name='environment', bond='Cu-P', target=3)
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='compute_mean_coordination', callable_params={"type_A": 29,
                                                                                                         "type_B": 15},
                                    max_bincenter=10)
        self._assert_histogram_inputs(dirs, ds_dict, hdc)
        fig, _ax = plot_multihistogram(multidata=hdc, target=3, max_bincenter=10)
        self.assertGreater(len(fig.axes), 0)
        plt.close(fig)

    @integration_test
    def test_bfend(self):
        system = "B-Fe-Nd"
        bond = "B-Fe"
        target = 3
        type_A, type_B = (Element(e).Z for e in bond.split('-'))
        dirs = get_environment_gen_dirs(self.processed_path, system=system, guidance_name='environment', bond=bond,
                                        target=target)
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='compute_mean_coordination',
                                    callable_params={"type_A": type_A, "type_B": type_B})
        self._assert_histogram_inputs(dirs, ds_dict, hdc)
        fig, _ax = plot_multihistogram(multidata=hdc, target=target, max_bincenter=10)
        self.assertGreater(len(fig.axes), 0)
        plt.close(fig)

    @integration_test
    def test_cusipca(self):
        system = "Cu-Si-P-Ca"
        bond = "Cu-P"
        target = 4
        type_A, type_B = (Element(e).Z for e in bond.split('-'))
        dirs = get_environment_gen_dirs(self.processed_path, system=system, guidance_name='environment', bond=bond,
                                        target=target)
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='compute_mean_coordination',
                                    callable_params={"type_A": type_A, "type_B": type_B})
        self._assert_histogram_inputs(dirs, ds_dict, hdc)
        fig, ax = plot_multihistogram(multidata=hdc, target=target, max_bincenter=10)
        ax = plt.gca()
        ax.set_box_aspect(1)
        self.assertGreater(len(fig.axes), 0)
        plt.close(fig)

    @integration_test
    def test_cusipca_parse_raw(self):
        system = "Cu-Si-P-Ca"
        bond = "Cu-P"
        target = 4
        type_A, type_B = (Element(e).Z for e in bond.split('-'))
        dirs = get_environment_gen_dirs(self.processed_path, system=system, guidance_name='environment', bond=bond,
                                        target=target)
        ds_dict = collect_stage_dataset_dict(dirs, "parse_raw", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='compute_mean_coordination',
                                    callable_params={"type_A": type_A, "type_B": type_B})
        self._assert_histogram_inputs(dirs, ds_dict, hdc)
        fig, ax = plot_multihistogram(multidata=hdc, target=target, max_bincenter=10)
        ax = plt.gca()
        ax.set_box_aspect(1)
        self.assertGreater(len(fig.axes), 0)
        plt.close(fig)

    @integration_test
    def test_auto_bins(self):
        target = 6.8
        dirs = get_volume_pa_gen_dirs(self.processed_path, 'B', 'volume_pa', target=target)
        ds_dict = collect_stage_dataset_dict(dirs, "deduplicate_all", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='volume_pa', callable_params={}, auto_adjust_bins=True,
                                    n_bins=10, integer_bins=False)
        self._assert_histogram_inputs(dirs, ds_dict, hdc)
        fig, _ax = plot_multihistogram(multidata=hdc, target=target, max_bincenter=10, show_gaussian=True)
        self.assertGreater(len(fig.axes), 0)
        plt.close(fig)

    @integration_test
    def test_custom_bins(self):
        target = 6.8
        dirs = get_volume_pa_gen_dirs(self.processed_path, 'B', 'volume_pa', target=target)
        ds_dict = collect_stage_dataset_dict(dirs, "deduplicate_all", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='volume_pa', callable_params={}, auto_adjust_bins=False,
                                    bin_centers=np.linspace(7.5, 12, 10), n_bins=10, integer_bins = False)
        self._assert_histogram_inputs(dirs, ds_dict, hdc)
        fig, _ax = plot_multihistogram(multidata=hdc, target=target, max_bincenter=10, show_gaussian=True)
        self.assertGreater(len(fig.axes), 0)
        plt.close(fig)

    @integration_test
    def test_histogram_correctness(self):
        system = "Si-O"
        bond = "Si-O"
        target = 4
        dirs = get_environment_gen_dirs(self.processed_path, system=system, guidance_name='environment', bond=bond,
                                        target=target)
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)

        type_A, type_B = (Element(el).Z for el in bond.split('-'))
        callable_name = 'compute_mean_coordination'
        callable_params = {"type_A": type_A, "type_B": type_B}

        hdc = histo_data_collection(ds_dict, callable_name=callable_name, callable_params=callable_params)
        idx_non_guided = [i for i, hist_data in enumerate(hdc) if hist_data['label'] == 'Non-guided' ][0]
        target_bin_counts = hdc[idx_non_guided]['counts'][hdc[idx_non_guided]['bin_centers']==target]
        self.assertGreaterEqual(target_bin_counts.size, 1)
        fig, _ax = plot_multihistogram(multidata=hdc, target=target, max_bincenter=10)
        function = get_target_value_fn(callable_name, force_gpu=0, **callable_params)
        fract_dict = print_ds_dict_guidance_resume(ds_dict, function, target=target, half_width=0.5)
        self.assertIn("Non-guided", fract_dict)
        self.assertTrue(all(0 <= fraction <= 1 for fraction in fract_dict.values()))
        self.assertGreater(len(fig.axes), 0)
        plt.close(fig)

    @integration_test
    def test_KDE(self):
        system = "B-Fe-Nd"
        bond = "B-Fe"
        target = 3
        type_A, type_B = (Element(e).Z for e in bond.split('-'))
        dirs = get_environment_gen_dirs(self.processed_path, system=system, guidance_name='environment', bond=bond,
                                        target=target)
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)
        vals_dict = calculate_values(ds_dict, callable_name='compute_mean_coordination',
                                    callable_params={"type_A": type_A, "type_B": type_B})
        self.assertGreater(len(dirs), 0)
        self.assertEqual(set(vals_dict), set(ds_dict))
        self.assertTrue(any(len(values) > 0 for values in vals_dict.values()))
        for values in vals_dict.values():
            self.assertTrue(np.all(np.isfinite(values)))
        fig, _ax = plot_multi_kde(values_dict=vals_dict, target=target)
        self.assertGreater(len(fig.axes), 0)
        plt.close(fig)

    # def test_empty_histogram(self):
    #     system = "Cu-Si-P-Ca"
