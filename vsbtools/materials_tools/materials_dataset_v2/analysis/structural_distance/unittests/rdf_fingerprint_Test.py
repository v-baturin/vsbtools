import unittest
from pathlib import Path
from ..rdf_fingerprint import RadialFingerprint, compute_fingerprint
from ....io.structures_dataset_io import  StructureDatasetIO
from ....converters import cell_pos_atomtypes_from_pmg_structure

PATH_WITH_TESTS = Path(__file__).parent
PATH_WITH_DATASETS = PATH_WITH_TESTS / "../../../unittests_datasets"
TWO_BORONS = PATH_WITH_TESTS / "../../../analysis/structural_distance/unittests/two_borons"

class USPEXBridge_Test(unittest.TestCase):

    def setUp(self):
        # self.struc_path = PATH_WITH_DATASETS/"two_systems"
        self.struc_path = TWO_BORONS
        self.sys1, self.sys2\
            = [e.structure for e in StructureDatasetIO(self.struc_path).load_from_directory()]


    def test_rdf_fingerprint(self):
        cell1, pos1, types1 = cell_pos_atomtypes_from_pmg_structure(self.sys1)
        cell2, pos2, types2 = cell_pos_atomtypes_from_pmg_structure(self.sys2)
        fp1 = compute_fingerprint(cell1, pos1, types1, [True] * 3)
        fp2 = compute_fingerprint(cell2, pos2, types2, [True] * 3)
        print(RadialFingerprint.cosine_distance(fp1, fp2))
