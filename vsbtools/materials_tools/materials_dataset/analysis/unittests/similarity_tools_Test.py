import unittest
from pathlib import Path
from ..similarity_tools import SimilarityTools
from ...io.uspex_bridge import USPEXBridge
from ...io.structures_dataset_io import StructureDatasetIO

PATH_WITH_TESTS = Path(__file__).parent
PATH_WITH_DATASETS = PATH_WITH_TESTS / "../../unittests_datasets"

class SimilarityTools_Test(unittest.TestCase):

    def setUp(self):
        self.test_structures = PATH_WITH_DATASETS / "dup_test_poscars"
        self.ds = StructureDatasetIO(self.test_structures, patterns=("*.POSCAR",)).load_from_directory()
        self.ub = USPEXBridge(elements={'Fe', 'Al', 'Ni'}, legacy=True, tol_FP=0.07)
        self.sim_tools = SimilarityTools(self.ub.fp_dist, self.ub.tol_FP)


    def test_uspex_based_distance(self):
        self.assertEqual(len(self.ds), 347)
        clusters_file = PATH_WITH_TESTS / "clusters.pkl"
        dist_matrix_file = PATH_WITH_TESTS / "dm.pkl"
        deduped_ds, _, _ = self.sim_tools.deduplicate(self.ds, enforce_compositions_separation=True,
                                                      clusters_file = clusters_file,
                                                      dist_matrix_file = dist_matrix_file,
                                                      check_clusters_file=True, check_dist_matrix_file=True)
        self.assertEqual(len(deduped_ds),299)
        clusters_file.unlink()
        dist_matrix_file.unlink()

    def test_suspiciousDistance(self):
        ds_sio2 = StructureDatasetIO(PATH_WITH_DATASETS / "SiO2_two_loose_structures", patterns=('*POSCAR',)).load_from_directory()
        ubSiO = USPEXBridge(elements={'Si', 'O'}, legacy=True, tol_FP=0.07)
        simtoolsSiO = SimilarityTools(ubSiO.fp_dist, ubSiO.tol_FP)
        print(simtoolsSiO.dist(ds_sio2[0], ds_sio2[1]))




