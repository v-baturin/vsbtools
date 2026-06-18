import unittest
import os
from pathlib import Path
import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

from ..similarity_tools import SimilarityTools
from .. import generation_postprocess as gp
from ..scenario_pipeline import Context
from ..structural_distance.dscribe_bridge import DScribeBridge, dscribe_uspex_mbtr_kwargs
from ...io.uspex_bridge import USPEXBridge
from ...io.structures_dataset_io import StructureDatasetIO
from pymatgen.core.periodic_table import Element

PATH_WITH_TESTS = Path(__file__).parent
PATH_WITH_DATASETS = PATH_WITH_TESTS / "../../unittests_datasets"

class SimilarityTools_Test(unittest.TestCase):

    def setUp(self):
        self.test_structures = PATH_WITH_DATASETS / "dup_test_poscars"
        self.ds = StructureDatasetIO(self.test_structures, pattern="*POSCAR*").load_from_directory()
        self.default_bridge = DScribeBridge(elements={'Fe', 'Al', 'Ni'}, tol_FP=0.07)
        self.sim_tools = SimilarityTools(self.default_bridge.fp_dist, self.default_bridge.tol_FP)


    def test_default_dscribe_based_distance(self):
        self.assertEqual(len(self.ds), 347)
        clusters_file = PATH_WITH_TESTS / "clusters.pkl"
        dist_matrix_file = PATH_WITH_TESTS / "dm.pkl"
        deduped_ds, _, _ = self.sim_tools.deduplicate(self.ds, enforce_compositions_separation=True,
                                                      clusters_file = clusters_file,
                                                      dist_matrix_file = dist_matrix_file,
                                                      check_clusters_file=True, check_dist_matrix_file=True)
        self.assertEqual(len(deduped_ds), 294)
        clusters_file.unlink()
        dist_matrix_file.unlink()

    def test_suspiciousDistance(self):
        ds_sio2 = StructureDatasetIO(PATH_WITH_DATASETS / "SiO2_two_loose_structures", pattern='*POSCAR').load_from_directory()
        ubSiO = USPEXBridge(elements={'Si', 'O'}, legacy=True, tol_FP=0.07)
        simtoolsSiO = SimilarityTools(ubSiO.fp_dist, ubSiO.tol_FP)
        print(simtoolsSiO.dist(ds_sio2[0], ds_sio2[1]))

    def test_suspiciousDistance2(self):
        fe3ndb2_1 = StructureDatasetIO(PATH_WITH_DATASETS / "fe3ndb2_close_sruc", pattern='*POSCAR').load_from_directory()
        ubSiO = USPEXBridge(elements={'Fe', 'Nd', 'B'}, legacy=True, tol_FP=0.07)
        simtoolsSiO = SimilarityTools(ubSiO.fp_dist, ubSiO.tol_FP)
        print(simtoolsSiO.dist(fe3ndb2_1[0], fe3ndb2_1[1]))

    @staticmethod
    def _species_by_atomic_number(entry):
        return sorted({site.specie.symbol for site in entry.structure}, key=lambda symbol: Element(symbol).Z)

    @staticmethod
    def _dscribe_uspex_equivalent_fingerprint(entry, species, rdu):
        """
        DScribe's ValleOganov shortcut fixes the grid to [0, Rmax]. The
        DScribeBridge preset uses MBTR directly, with Valle-Oganov normalization
        sampled at USPEX bin centers.
        """
        bridge = DScribeBridge(
            elements=species,
            preset="uspex",
            Rmax=rdu.Rmax,
            sigma=rdu.sigma,
            delta=rdu.delta,
        )
        return bridge.uspex_fingerprint_channels(entry)

    def test_uspex_fingerprint_matches_dscribe_mbtr_mapping(self):
        entries = [
            StructureDatasetIO(PATH_WITH_DATASETS / "two_borons").load_from_directory()[0],
            StructureDatasetIO(PATH_WITH_DATASETS / "two_systems_MgAlO").load_from_directory()[0],
        ]

        for entry in entries:
            with self.subTest(entry_id=entry.id, formula=entry.structure.composition.reduced_formula):
                species = self._species_by_atomic_number(entry)
                uspex_bridge = USPEXBridge(elements=set(species), legacy=True)
                uspex_fingerprint = uspex_bridge.uspex_entry_from_de(entry)[
                    "radialDistributionUtility.structureFingerprint.origin"
                ]
                dscribe_fingerprint = self._dscribe_uspex_equivalent_fingerprint(
                    entry, species, uspex_bridge.rdu
                )

                for pair, dscribe_values in dscribe_fingerprint.items():
                    with self.subTest(pair=pair):
                        np.testing.assert_allclose(
                            uspex_fingerprint[pair],
                            dscribe_values,
                            atol=1.5e-4,
                            rtol=0,
                        )

    def test_dscribe_bridge_uspex_preset_distance_matches_uspex_bridge(self):
        entries = StructureDatasetIO(PATH_WITH_DATASETS / "two_borons").load_from_directory()
        species = self._species_by_atomic_number(entries[0])

        uspex_bridge = USPEXBridge(elements=set(species), legacy=True)
        dscribe_bridge = DScribeBridge(
            elements=species,
            preset="uspex",
            Rmax=uspex_bridge.rdu.Rmax,
            sigma=uspex_bridge.rdu.sigma,
            delta=uspex_bridge.rdu.delta,
        )

        self.assertAlmostEqual(
            dscribe_bridge.fp_dist(entries[0], entries[1]),
            uspex_bridge.fp_dist(entries[0], entries[1]),
            delta=2e-4,
        )

    def test_dscribe_bridge_accepts_custom_descriptor_cls(self):
        entry = StructureDatasetIO(PATH_WITH_DATASETS / "two_borons").load_from_directory()[0]
        species = self._species_by_atomic_number(entry)
        bridge = DScribeBridge(
            descriptor_cls="MBTR",
            descriptor_kwargs=dscribe_uspex_mbtr_kwargs(species),
        )

        self.assertAlmostEqual(bridge.fp_dist(entry, entry), 0.0)

    def test_scenario_similarity_defaults_to_dscribe_bridge(self):
        ctx = Context()
        ctx.globals["elements"] = {"B"}

        similarity = ctx.get_tool("similarity")

        self.assertIsInstance(ctx.toolkits["dscribe"], DScribeBridge)
        self.assertIs(similarity.dist.__self__, ctx.toolkits["dscribe"])

    def test_scenario_similarity_can_still_use_uspex_bridge(self):
        ctx = Context()
        ctx.globals["elements"] = {"B"}
        ctx.toolkit_options["similarity"]["bridge"] = "uspex"

        similarity = ctx.get_tool("similarity")

        self.assertIsInstance(ctx.toolkits["uspex"], USPEXBridge)
        self.assertIs(similarity.dist.__self__, ctx.toolkits["uspex"])

    def test_legacy_pipeline_similarity_defaults_to_dscribe_bridge(self):
        pipeline = gp.PPPipeline(elements={"B"})

        similarity = pipeline.get_tool("similarity")

        self.assertIsInstance(pipeline.toolkits["dscribe"], DScribeBridge)
        self.assertIs(similarity.dist.__self__, pipeline.toolkits["dscribe"])
