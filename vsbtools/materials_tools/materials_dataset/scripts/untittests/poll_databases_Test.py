import unittest
from pathlib import Path
from ...io.yaml_csv_poscars import write, read
from ..poll_databases import poll_databases
from ...converters import ds2df
from ...analysis.summary import collect_summary_df, print_pretty_df


PATH_WITH_TESTS = Path(__file__).parent


class yaml_csv_poscars_Test(unittest.TestCase):

    def setUp(self) -> None:
        self.elements = {'Mo', 'Si', 'B', 'P'}

    def test_scrape(self):
        ds = poll_databases(self.elements, do_ehull_filtering=True,
                            loader_kwargs={"alexandria": {'pattern': 'alexandria_00*.json'}})
        write(ds, enforce_base_path=PATH_WITH_TESTS / "ds")
        print_pretty_df(collect_summary_df(ds, native_columns=("id", "composition", "energy", "metadata.duplicates")), PATH_WITH_TESTS / 'table.txt')
        ds_read = read(PATH_WITH_TESTS / 'ds/manifest.yaml')
        print_pretty_df(collect_summary_df(ds_read, native_columns=("id", "composition", "energy", "metadata.duplicates")),
                        PATH_WITH_TESTS / 'table2.txt')
        write(ds_read, enforce_base_path=PATH_WITH_TESTS / "ds2")
        self.assertEqual(len(ds), 65)