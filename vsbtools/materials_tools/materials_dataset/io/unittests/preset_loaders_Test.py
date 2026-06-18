import unittest
from pathlib import Path
from unittest.mock import Mock, patch
from ...io.preset_loaders import (load_from_materials_project,
                                                                                load_from_oqmd,
                                                                                load_from_alexandria,
                                                                                load_from_optimade,
                                                                                load_from_uspex_goodstructures,
                                                                                load_from_uspex_calc_folders,
                                                                                load_mattersim_estimated_set)

PATH_WITH_TESTS = Path(__file__).parent
PATH_WITH_DATASETS = PATH_WITH_TESTS / "../../unittests_datasets"


class DSLoaders_Test(unittest.TestCase):

    def setUp(self) -> None:
        self.poscars_folder = PATH_WITH_DATASETS / "POSCARS"
        self.uspex_calcfolders = PATH_WITH_DATASETS / 'uspex_folders/calcFolders'
        self.uspex_goodStructures = PATH_WITH_DATASETS / 'uspex_folders/results1/goodStructures'
        self.csv = PATH_WITH_TESTS / 'mattersim_res_for_POSCARS.csv'

    def test_loadMP(self):
        ds = load_from_materials_project({'Mo', 'Si'})
        self.assertEqual(len(ds), 58)

    def test_load_from_materials_project_interface(self):
        query_df = Mock(name="query_df")
        dataset = Mock(name="dataset")
        client = Mock()
        client.query.return_value = query_df

        with patch("vsbtools.materials_tools.materials_dataset.io.preset_loaders.mp_client", client), \
                patch("vsbtools.materials_tools.materials_dataset.io.preset_loaders.df2ds",
                      return_value=dataset) as df2ds:
            result = load_from_materials_project.__wrapped__({'Mo', 'Si'})

        client.query.assert_called_once_with({'Mo', 'Si'})
        self.assertIn("Materials Project", df2ds.call_args.kwargs["message"])
        self.assertIs(result, dataset)

    def test_loadOQMD(self):
        ds = load_from_oqmd({'Mo', 'Si'})
        self.assertEqual(len(ds), 83)

    def test_load_from_oqmd_interface(self):
        query_df = Mock(name="query_df")
        dataset = Mock(name="dataset")
        client = Mock()
        client.query.return_value = query_df

        with patch("vsbtools.materials_tools.materials_dataset.io.preset_loaders.oqmd_client", client), \
                patch("vsbtools.materials_tools.materials_dataset.io.preset_loaders.df2ds",
                      return_value=dataset) as df2ds:
            result = load_from_oqmd.__wrapped__({'Mo', 'Si'})

        client.query.assert_called_once_with({'Mo', 'Si'})
        self.assertIn("OQMD", df2ds.call_args.kwargs["message"])
        self.assertIs(result, dataset)

    def test_loadAlexandria(self):
        ds = load_from_alexandria({'Mo', 'Si'}, pattern='alexandria_00*.json')
        self.assertEqual(len(ds), 34)

    def test_load_from_alexandria_interface(self):
        query_df = Mock(name="query_df")
        dataset = Mock(name="dataset")
        client = Mock()
        client.query.return_value = query_df

        with patch("vsbtools.materials_tools.materials_dataset.io.preset_loaders.AlexandriaClient",
                   return_value=client) as client_cls, \
                patch("vsbtools.materials_tools.materials_dataset.io.preset_loaders.df2ds",
                      return_value=dataset) as df2ds:
            result = load_from_alexandria.__wrapped__({'Mo', 'Si'}, pattern="alexandria_00*.json")

        client_cls.assert_called_once_with(pattern="alexandria_00*.json")
        client.query.assert_called_once_with({'Mo', 'Si'})
        self.assertIn("Alexandria", df2ds.call_args.kwargs["message"])
        self.assertIs(result, dataset)

    def test_load_from_optimade(self):
        query_df = Mock(name="query_df")
        dataset = Mock(name="dataset")
        client = Mock()
        client.query.return_value = query_df

        with patch("vsbtools.materials_tools.materials_dataset.io.preset_loaders.OptimadeClient",
                   return_value=client) as client_cls, \
                patch("vsbtools.materials_tools.materials_dataset.io.preset_loaders.df2ds",
                      return_value=dataset) as df2ds:
            result = load_from_optimade.__wrapped__({'Mo', 'Si'}, providers=["oqmd"], page_limit=10)

        client_cls.assert_called_once_with(providers=["oqmd"], page_limit=10)
        client.query.assert_called_once_with({'Mo', 'Si'})
        self.assertTrue(df2ds.call_args.kwargs["message"].startswith("Full "))
        self.assertIn("OPTIMADE", df2ds.call_args.kwargs["message"])
        self.assertIs(result, dataset)

    def test_load_from_uspex_goodstructures(self):
        ds = load_from_uspex_goodstructures(self.uspex_goodStructures)
        self.assertEqual(len(ds), 442)

    def test_load_from_uspex_calc_folders(self):
        ds = load_from_uspex_calc_folders(self.uspex_calcfolders)
        self.assertEqual(len(ds), 214)

    def test_load_mattersim_estimated_set(self):
        ds = load_mattersim_estimated_set(self.csv, self.poscars_folder)
        self.assertEqual(len(ds), 7142)


    # def test_elements(self):
    #     ds = CrystalDataset.from_struct_folder(self.poscars_folder, search_pattern='*.POSCAR', skip_dump=True)
    #     self.assertEqual(ds.elements, {'N', 'Na', 'Al', 'Ni', 'Fe'})
    #
    # def test_deduplication(self):
    #     ds = CrystalDataset.from_struct_folder(self.poscars_folder, search_pattern='*.POSCAR', skip_dump=True)
    #     assert len(ds) == 348, "Dataset should contain 348 entries before deduplication"
    #     ds.tol_FP = 0.07
    #     deduped_ds, _, _ = ds.deduplicated(enforce_compositions_separation=True,
    #                                         fitness_list=None,
    #                                         reset_entries=False, skip_dump=True)
    #     self.assertEqual(len(deduped_ds),300)
    #
    # def test_filtering(self):
    #     ds = CrystalDataset.from_client(OQMDClient(), elements=['Al', 'Fe', 'Ni'], skip_dump=True)
    #     fs = ds.filter(predicate_fn=lambda e: e.e_above_hull < 0.02, reset_entries=False, skip_dump=True)
    #     print(fs.metadata["message"])
    #
    # def test_rw(self):
    #     ds = CrystalDataset.from_struct_folder(self.poscars_folder, search_pattern='*.POSCAR', skip_dump=True)
    #     ds.dump()
    #     with open(ds.pkl_path, 'rb') as f:
    #         ds2 = pkl.load(f)
    #     print(f"ds1.id: {ds.id}, ds2.id: {ds2.id}")
    #     self.assertEqual(ds.id, ds2.id)
    #     Path(ds.pkl_path).unlink()
    #     Path(ds.regfile).unlink()
    #     Path(ds.treefile).unlink()
