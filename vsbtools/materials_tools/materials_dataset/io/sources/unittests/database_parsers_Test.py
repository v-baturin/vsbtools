import time
import unittest

try:
    from ..matproj_parser import MPClient
except ImportError:
    MPClient = None
from ..alexandria_parser import AlexandriaClient
from ..oqmd_parser import OQMDClient
from ..optimade_parser import OptimadeClient


def _strip_optimade_prefix(entry_id, provider):
    return str(entry_id).removeprefix(f"optimade_{provider}_")


def _normalize_source_id(entry_id, provider):
    entry_id = _strip_optimade_prefix(entry_id, provider)
    if provider == "oqmd":
        return entry_id.removeprefix("oqmd_").removeprefix("oqmd-")
    return entry_id


def _normalized_ids(df, provider):
    return {_normalize_source_id(entry_id, provider) for entry_id in df["id"]}


def _timed_call(fn):
    start = time.perf_counter()
    result = fn()
    return result, time.perf_counter() - start


def _short_id_diff_report(provider, direct_ids, optimade_ids, direct_time, optimade_time):
    direct_only = sorted(direct_ids - optimade_ids)
    optimade_only = sorted(optimade_ids - direct_ids)
    common = direct_ids & optimade_ids
    return (
        f"{provider} direct/OPTIMADE ID report: "
        f"direct={len(direct_ids)}, optimade={len(optimade_ids)}, common={len(common)}, "
        f"direct_only={len(direct_only)}, optimade_only={len(optimade_only)}, "
        f"direct_time={direct_time:.3f}s, optimade_time={optimade_time:.3f}s, "
        f"direct_only_sample={direct_only[:10]}, optimade_only_sample={optimade_only[:10]}"
    )


def _optimade_query(provider, elements):
    return OptimadeClient(
        providers=[provider],
        do_deduplication=False,
        page_limit=500,
        timeout=12000,
    ).query(elements)


class DatabaseParsersTest(unittest.TestCase):
    def test_si_oqmd_local_ids_are_subset_of_optimade(self):
        elements = {"Si"}

        try:
            optimade_df, optimade_time = _timed_call(lambda: _optimade_query("oqmd", elements))
        except Exception as err:
            oqmd_client = OQMDClient()
            try:
                _timed_call(lambda: oqmd_client.query(elements))
            except Exception as local_err:
                self.skipTest(f"OQMD OPTIMADE endpoint and local OQMD are unavailable: {err}; {local_err}")
            finally:
                conn = getattr(oqmd_client, "_conn", None)
                if conn is not None:
                    conn.close()
            self.skipTest(f"OQMD OPTIMADE endpoint is not available; local OQMD fallback works: {err}")

        oqmd_client = OQMDClient()
        try:
            direct_df, direct_time = _timed_call(lambda: oqmd_client.query(elements))
        except Exception as err:
            self.skipTest(f"Local OQMD is not available: {err}")
        finally:
            conn = getattr(oqmd_client, "_conn", None)
            if conn is not None:
                conn.close()

        direct_ids = _normalized_ids(direct_df, "oqmd")
        optimade_ids = _normalized_ids(optimade_df, "oqmd")
        report = _short_id_diff_report("OQMD", direct_ids, optimade_ids, direct_time, optimade_time)
        print(report)

        self.assertGreater(len(direct_ids), 0, report)
        self.assertGreater(len(optimade_ids), 0, report)
        self.assertNotAlmostEqual(direct_time, optimade_time, delta=0.05, msg=report)
        self.assertTrue(direct_ids <= optimade_ids, report)

    def test_si_materials_project_direct_and_optimade_ids_match(self):
        if MPClient is None:
            self.skipTest("mp-api is not available")

        elements = {"Si"}
        try:
            optimade_df, optimade_time = _timed_call(lambda: _optimade_query("materials_project", elements))
        except Exception as err:
            try:
                _timed_call(lambda: MPClient().query(elements))
            except Exception as direct_err:
                self.skipTest(
                    f"Materials Project OPTIMADE and direct API are unavailable: {err}; {direct_err}"
                )
            self.skipTest(f"Materials Project OPTIMADE endpoint is not available; direct API fallback works: {err}")

        try:
            direct_df, direct_time = _timed_call(lambda: MPClient().query(elements))
        except Exception as err:
            self.skipTest(f"Materials Project direct API is not available: {err}")

        direct_ids = _normalized_ids(direct_df, "materials_project")
        optimade_ids = _normalized_ids(optimade_df, "materials_project")
        report = _short_id_diff_report("Materials Project", direct_ids, optimade_ids, direct_time, optimade_time)
        print(report)

        self.assertGreater(len(direct_ids), 0, report)
        self.assertGreater(len(optimade_ids), 0, report)
        self.assertNotAlmostEqual(direct_time, optimade_time, delta=0.05, msg=report)
        self.assertEqual(direct_ids, optimade_ids, report)

    def test_si_alexandria_local_ids_are_subset_of_optimade(self):
        elements = {"Si"}
        try:
            optimade_df, optimade_time = _timed_call(lambda: _optimade_query("alexandria", elements))
        except Exception as err:
            try:
                _timed_call(lambda: AlexandriaClient(prompt=False).query(elements))
            except Exception as local_err:
                self.skipTest(
                    f"Alexandria OPTIMADE endpoint and local Alexandria snapshot are unavailable: "
                    f"{err}; {local_err}"
                )
            self.skipTest(f"Alexandria OPTIMADE endpoint is not available; local fallback works: {err}")

        try:
            direct_df, direct_time = _timed_call(lambda: AlexandriaClient(prompt=False).query(elements))
        except Exception as err:
            self.skipTest(f"Local Alexandria snapshot is not available: {err}")

        direct_ids = _normalized_ids(direct_df, "alexandria")
        optimade_ids = _normalized_ids(optimade_df, "alexandria")
        report = _short_id_diff_report("Alexandria", direct_ids, optimade_ids, direct_time, optimade_time)
        print(report)

        self.assertGreater(len(direct_ids), 0, report)
        self.assertGreater(len(optimade_ids), 0, report)
        self.assertLess(optimade_time, direct_time * 0.75, report)
        self.assertTrue(direct_ids <= optimade_ids, report)


if __name__ == "__main__":
    unittest.main()
