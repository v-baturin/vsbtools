import shutil
import tempfile
import time
import unittest
import multiprocessing as mp
from pathlib import Path

try:
    from ...crystal_dataset import CrystalDataset
    from ...io.structures_dataset_io import StructureDatasetIO
    from ...io.zip_handling import exploded_zip_tree
    from ...analysis.similarity_tools import SimilarityTools
    from ...analysis.structural_distance.dscribe_bridge import DScribeBridge
    from ....geom_utils.structure_checks import check_min_dist_pmg
except ImportError as exc:
    CrystalDataset = None
    IMPORT_ERROR = exc
else:
    IMPORT_ERROR = None


REPO_ROOT = Path(__file__).resolve().parents[5]
RAW_GENERATION_ZIP = (
        REPO_ROOT
        / "vsbtools/materials_dataset/Examples/raw_generations"
        / "Cu-Si-P/Cu-Si-P_guided_CuP3.zip"
)
PAIR_DISTANCE_TIMEOUT_SECONDS = 10.0


def _load_filtered_cusip_entries(tmp_root: Path):
    shutil.copy2(RAW_GENERATION_ZIP, tmp_root / RAW_GENERATION_ZIP.name)
    with exploded_zip_tree(tmp_root) as expanded_root:
        raw = StructureDatasetIO(
            expanded_root,
            pattern="*.extxyz",
            source_name="Mattergen",
        ).load_from_directory(elements={"Cu", "Si", "P"})

    filtered = raw.filter(lambda entry: check_min_dist_pmg(entry.structure)[0])
    return raw, filtered


def _time_distance_364_381(queue):
    try:
        with tempfile.TemporaryDirectory(prefix="cusip_guided_cup3_pair_") as tmp:
            _, filtered = _load_filtered_cusip_entries(Path(tmp))
            entries_by_id = {entry.id: entry for entry in filtered}
            bridge = DScribeBridge(elements={"Cu", "Si", "P"}, tol_FP=0.008)
            similarity = SimilarityTools(bridge.fp_dist, bridge.tol_FP)
            start = time.perf_counter()
            distance = similarity.dist(entries_by_id["364"], entries_by_id["381"])
            queue.put(("ok", distance, time.perf_counter() - start))
    except Exception as exc:
        queue.put(("error", repr(exc)))


@unittest.skipIf(IMPORT_ERROR is not None, f"materials dependencies unavailable: {IMPORT_ERROR}")
class CuSiPGuidedDedupTest(unittest.TestCase):

    def test_check_min_dist_then_deduplicate_cusip_guided_cup3(self):
        if not RAW_GENERATION_ZIP.exists():
            self.skipTest(f"missing fixture: {RAW_GENERATION_ZIP}")

        context = mp.get_context("spawn")
        queue = context.Queue()
        process = context.Process(target=_time_distance_364_381, args=(queue,))
        process.start()
        process.join(PAIR_DISTANCE_TIMEOUT_SECONDS)
        if process.is_alive():
            process.terminate()
            process.join()
            self.fail(
                f"distance between entries 364 and 381 exceeded "
                f"{PAIR_DISTANCE_TIMEOUT_SECONDS:.1f}s"
            )
        if queue.empty():
            self.fail(f"distance worker exited with code {process.exitcode} and no result")
        status, *payload = queue.get()
        if status == "error":
            self.fail(payload[0])
        distance_364_381, elapsed = payload
        self.assertIsInstance(distance_364_381, float)
        self.assertLess(elapsed, PAIR_DISTANCE_TIMEOUT_SECONDS)

        with tempfile.TemporaryDirectory(prefix="cusip_guided_cup3_") as tmp:
            tmp_root = Path(tmp)
            raw, filtered = _load_filtered_cusip_entries(tmp_root)

            self.assertGreater(len(raw), 381)
            self.assertLessEqual(len(filtered), len(raw))

            entries_by_id = {entry.id: entry for entry in filtered}
            self.assertIn("364", entries_by_id)
            self.assertIn("381", entries_by_id)

            bridge = DScribeBridge(elements={"Cu", "Si", "P"}, tol_FP=0.008)
            similarity = SimilarityTools(bridge.fp_dist, bridge.tol_FP)
            dedup_entries = list(entries_by_id.values())[:48]
            for entry_id in ("364", "381"):
                if entries_by_id[entry_id] not in dedup_entries:
                    dedup_entries.append(entries_by_id[entry_id])
            dedup_subset = CrystalDataset(
                dedup_entries,
                message="Cu-Si-P guided CuP3 min-distance filtered dedup subset",
            )

            dedup_dir = tmp_root / "dedup"
            dedup_dir.mkdir()
            dedup_subset.override_base_path(dedup_dir)
            deduplicated, _, _ = similarity.deduplicate(dedup_subset, tol_FP=0.008)
            self.assertLessEqual(len(deduplicated), len(dedup_subset))
            self.assertIsInstance(deduplicated, CrystalDataset)


if __name__ == "__main__":
    unittest.main()
