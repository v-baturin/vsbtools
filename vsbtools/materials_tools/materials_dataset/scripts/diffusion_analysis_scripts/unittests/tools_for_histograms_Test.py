from typing import Callable, Dict
import unittest
import torch
import numpy as np
from pathlib import Path
from pymatgen.core import Structure, Element
from ase.io import read as ase_read
from ....crystal_dataset import CrystalDataset
from ....crystal_entry import CrystalEntry
from ..guidance_stats import (get_environment_gen_dirs, get_volume_pa_gen_dirs, collect_stage_dataset_dict,
                              histo_data_collection, get_target_value_fn, plot_multihistogram, plot_multihistogram_new)
from matplotlib import pyplot as plt
PROCESSED_PATH = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/PROCESSED")


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
        self.sio2_stich = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/PROCESSED/O-Si/O-Si__guidance_environment_mode_huber_Si-O_6__diffusion_loss_weight_1-1-True__algo_0/2_x381aa8ac05cce031/POSCARS/agm002170463POSCAR")
        self.sio2_quartz = Path("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/MG_postprocess_pipelines/PROCESSED/O-Si/O-Si__guidance_environment_mode_huber_Si-O_6__diffusion_loss_weight_1-1-True__algo_0/2_x381aa8ac05cce031/POSCARS/agm002228342POSCAR")
        self.system_sio = "Si-O"
        self.system_licoo = "Li-Co-O"
        self.system_cupsi = "Cu-P-Si"
        self.system_bfend = "B-Fe-Nd"
        self.system_b = "B"

    def test_get_average_cn_gen_dirs(self):
        dirs = get_environment_gen_dirs(PROCESSED_PATH, self.system_sio, guidance_name='environment', bond='Si-O', target=4)
        self.assertEqual(len(dirs), 2)

    def test_get_entry_fn(self):
        fn = get_target_value_fn(fn_name="compute_mean_coordination", type_A=14, type_B=8)
        entry = CrystalEntry(id='0', structure=Structure.from_file(self.sio2_quartz))
        print(fn(entry))

        atoms = ase_read(self.sio2_quartz)

        cell = torch.tensor(np.array(atoms.cell), dtype=torch.float32,
                            device='cuda' if torch.cuda.is_available() else 'cpu')
        pos = torch.tensor(atoms.get_scaled_positions(), dtype=torch.float32, requires_grad=True, device=cell.device)
        atomic_numbers = torch.tensor(atoms.get_atomic_numbers(), dtype=torch.int64, device=cell.device)

        fn = get_target_value_fn(fn_name="_soft_neighbor_counts_per_A_single", type_A=14, type_B=8)
        cns = fn(entry)
        print(cns)

    def test_collect_stage_dataset_dict(self):
        dirs = get_environment_gen_dirs(PROCESSED_PATH, self.system_sio, guidance_name='environment', bond='Si-O', target=6)
        print(f"found {len(dirs)} dirs")
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='compute_mean_coordination', callable_params={"type_A": 14,
                                                                                                         "type_B": 8},
                                    max_bincenter=10)
        plot_multihistogram(multidata=hdc, target=6, max_bincenter=10)
        print(hdc)
        plt.show()

    def test_licoo(self):
        dirs = get_environment_gen_dirs(PROCESSED_PATH, self.system_licoo, guidance_name='environment', bond='Co-O', target=4)
        print(f"found {len(dirs)} dirs")
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='compute_mean_coordination', callable_params={"type_A": 27,
                                                                                                         "type_B": 8},
                                    max_bincenter=10)
        plot_multihistogram_new(multidata=hdc, target=4, max_bincenter=10)
        print(hdc)
        plt.show()


    def test_cupsi(self):
        dirs = get_environment_gen_dirs(PROCESSED_PATH, self.system_cupsi, guidance_name='environment', bond='Cu-P', target=3)
        print(f"found {len(dirs)} dirs")
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='compute_mean_coordination', callable_params={"type_A": 29,
                                                                                                         "type_B": 15},
                                    max_bincenter=10)
        plot_multihistogram(multidata=hdc, target=3, max_bincenter=10)
        print(hdc)
        plt.show()

    def test_bfend(self):
        system = "B-Fe-Nd"
        bond = "B-Fe"
        target = 3
        type_A, type_B = (Element(e).Z for e in bond.split('-'))
        dirs = get_environment_gen_dirs(PROCESSED_PATH, system=system, guidance_name='environment', bond=bond,
                                        target=target)
        print(f"found {len(dirs)} dirs")
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='compute_mean_coordination',
                                    callable_params={"type_A": type_A, "type_B": type_B})
        plot_multihistogram(multidata=hdc, target=target, max_bincenter=10)
        print(hdc)
        plt.show()

    def test_cusipca(self):
        system = "Cu-Si-P-Ca"
        bond = "Cu-P"
        target = 4
        type_A, type_B = (Element(e).Z for e in bond.split('-'))
        dirs = get_environment_gen_dirs(PROCESSED_PATH, system=system, guidance_name='environment', bond=bond,
                                        target=target)
        print(f"found {len(dirs)} dirs")
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='compute_mean_coordination',
                                    callable_params={"type_A": type_A, "type_B": type_B})
        plot_multihistogram_new(multidata=hdc, target=target, max_bincenter=10)
        print(hdc)
        plt.show()

    def test_auto_bins(self):
        target = 6.8
        dirs = get_volume_pa_gen_dirs(PROCESSED_PATH, 'B', 'volume_pa', target=target)
        ds_dict = collect_stage_dataset_dict(dirs, "deduplicate_all", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='volume_pa', callable_params={}, auto_adjust_bins=True,
                                    n_bins=10, integer_bins=False)
        plot_multihistogram(multidata=hdc, target=target, max_bincenter=10, show_gaussian=True)
        plt.show()
        print(list(dirs))

    def test_custom_bins(self):
        target = 6.8
        dirs = get_volume_pa_gen_dirs(PROCESSED_PATH, 'B', 'volume_pa', target=target)
        ds_dict = collect_stage_dataset_dict(dirs, "deduplicate_all", "poll_db", add_guid_descr=True)
        hdc = histo_data_collection(ds_dict, callable_name='volume_pa', callable_params={}, auto_adjust_bins=False,
                                    bin_centers=np.linspace(7.5, 12, 10), n_bins=10, integer_bins = False)
        plot_multihistogram(multidata=hdc, target=target, max_bincenter=10, show_gaussian=True)
        plt.show()
        print(list(dirs))

    def test_histogram_correctness(self):
        system = "Si-O"
        bond = "Si-O"
        target = 4
        dirs = get_environment_gen_dirs(PROCESSED_PATH, system=system, guidance_name='environment', bond=bond,
                                        target=target)
        print(f"found {len(dirs)} dirs")
        ds_dict = collect_stage_dataset_dict(dirs, "symmetrize_raw", "poll_db", add_guid_descr=True)

        type_A, type_B = (Element(el).Z for el in bond.split('-'))
        callable_name = 'compute_mean_coordination'
        callable_params = {"type_A": type_A, "type_B": type_B}

        hdc = histo_data_collection(ds_dict, callable_name=callable_name, callable_params=callable_params)
        idx_non_guided = [i for i, hist_data in enumerate(hdc) if hist_data['label'] == 'Non-guided' ][0]
        print(hdc[idx_non_guided]['counts'][hdc[idx_non_guided]['bin_centers']==target])
        # plot_multihistogram(multidata=hdc, target=target, max_bincenter=10)
        function = get_target_value_fn(callable_name, **callable_params)
        half_width = 0.5
        ds_name = 'Non-guided'
        fract_dict = print_ds_dict_guidance_resume(ds_dict, function, target=target, half_width=0.5)
        # TODO: finalize testcase, add self.assertequal or smth
        # target_entries = [e for e in ds_dict[ds_name] if target-half_width <= function(e) <= target + half_width]
        # print(len(target_entries))
        # plt.show()

    # def test_empty_histogram(self):
    #     system = "Cu-Si-P-Ca"