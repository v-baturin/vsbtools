import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch
from ...io.yaml_csv_poscars import write, read
from ..poll_databases import poll_databases
from ...analysis.summary import collect_summary_df, print_pretty_df


PATH_WITH_TESTS = Path(__file__).parent


class _FakeEntry:
    def __init__(self, entry_id, energy=-1.0):
        self.id = entry_id
        self.energy = energy

    def copy_with(self, **kw):
        return _FakeEntry(self.id, energy=kw.get("energy", self.energy))


class _FakeDataset:
    def __init__(self, entries):
        self.entries = list(entries)
        self.metadata = {}

    def __iter__(self):
        return iter(self.entries)

    def __len__(self):
        return len(self.entries)


class yaml_csv_poscars_Test(unittest.TestCase):

    def setUp(self) -> None:
        self.elements = {'La', 'Te', 'C'} # , 'P'}

    def test_scrape(self):
        ds = poll_databases(self.elements, do_ehull_filtering=True,
                            loader_kwargs={"alexandria": {'pattern': 'alexandria_00*.json'}})
        write(ds, enforce_base_path=PATH_WITH_TESTS / "ds3")
        print_pretty_df(collect_summary_df(ds, native_columns=("id", "composition", "energy", "metadata.duplicates")), PATH_WITH_TESTS / 'table.txt')
        ds_read = read(PATH_WITH_TESTS / 'ds3/manifest.yaml')
        print_pretty_df(collect_summary_df(ds_read, native_columns=("id", "composition", "energy", "metadata.duplicates")),
                        PATH_WITH_TESTS / 'table2.txt')
        write(ds_read, enforce_base_path=PATH_WITH_TESTS / "ds2")
        self.assertEqual(len(ds), 65)

    def test_poll_databases_accepts_optimade_source(self):
        optimade_loader = Mock(return_value=_FakeDataset([_FakeEntry("op1")]))

        with TemporaryDirectory() as tmpdir, \
                patch("vsbtools.materials_tools.materials_dataset.scripts.poll_databases.LOADERS",
                      {"op": optimade_loader}), \
                patch("vsbtools.materials_tools.materials_dataset.scripts.poll_databases.write") as write_mock:
            ds = poll_databases(
                {"Si"},
                database_names=["optimade"],
                pref_db="op",
                do_ehull_filtering=False,
                do_deduplication=False,
                loader_kwargs={"optimade": {"providers": ["oqmd"], "page_limit": 3}},
                cache_root_path=Path(tmpdir),
            )

        optimade_loader.assert_called_once_with({"Si"}, providers=["oqmd"], page_limit=3)
        write_mock.assert_called_once()
        self.assertEqual(len(ds), 1)
        self.assertIsNone(ds[0].energy)

    def test_poll_databases_falls_back_to_optimade_when_local_source_fails(self):
        local_loader = Mock(side_effect=FileNotFoundError("local snapshot missing"))
        optimade_loader = Mock(return_value=_FakeDataset([_FakeEntry("op1")]))

        with TemporaryDirectory() as tmpdir, \
                patch("vsbtools.materials_tools.materials_dataset.scripts.poll_databases.LOADERS",
                      {"al": local_loader, "op": optimade_loader}), \
                patch("vsbtools.materials_tools.materials_dataset.scripts.poll_databases.load_from_optimade",
                      optimade_loader), \
                patch("vsbtools.materials_tools.materials_dataset.scripts.poll_databases.write"):
            ds = poll_databases(
                {"Si"},
                database_names=["alexandria"],
                pref_db="al",
                do_ehull_filtering=False,
                do_deduplication=False,
                loader_kwargs={"op": {"page_limit": 2}},
                cache_root_path=Path(tmpdir),
            )

        local_loader.assert_called_once_with({"Si"}, prompt=False)
        optimade_loader.assert_called_once_with({"Si"}, page_limit=2, providers=["alexandria"])
        self.assertEqual(len(ds), 1)
        self.assertIsNone(ds[0].energy)
