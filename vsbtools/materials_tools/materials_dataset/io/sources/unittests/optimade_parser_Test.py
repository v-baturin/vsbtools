import unittest
from unittest.mock import patch

from ..optimade_parser import OptimadeClient


class _FakeComposition:
    formula = "Si1"


class _FakeStructure:
    def __init__(self, lattice, species, positions, coords_are_cartesian=False):
        self.lattice = lattice
        self.species = species
        self.positions = positions
        self.coords_are_cartesian = coords_are_cartesian
        self.composition = _FakeComposition()

    def __len__(self):
        return len(self.species)


def _optimade_item(entry_id="1", **attrs):
    base_attrs = {
        "chemical_formula_reduced": "Si",
        "lattice_vectors": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "cartesian_site_positions": [[0, 0, 0]],
        "species_at_sites": ["Si"],
        "species": [{"name": "Si", "chemical_symbols": ["Si"], "concentration": [1.0]}],
        "structure_features": [],
        "nsites": 1,
    }
    base_attrs.update(attrs)
    return {"id": entry_id, "type": "structures", "attributes": base_attrs}


class OptimadeClientTest(unittest.TestCase):

    def test_unknown_provider_raises(self):
        client = OptimadeClient(providers=["not-a-provider"])
        with self.assertRaisesRegex(ValueError, "Unknown OPTIMADE provider"):
            client.query({"Si"})

    def test_elements_filter_uses_exact_chemsys(self):
        self.assertEqual(
            OptimadeClient._elements_filter(("O", "Si")),
            'elements HAS ALL "O", "Si" AND nelements=2',
        )

    def test_next_url_preserves_https(self):
        next_url = OptimadeClient._next_url(
            "https://optimade.materialsproject.org/v1/structures?page_limit=1",
            {"links": {"next": "http://optimade.materialsproject.org/v1/structures?page_offset=1"}},
        )
        self.assertTrue(next_url.startswith("https://"))

    @patch("vsbtools.materials_tools.materials_dataset.io.sources.optimade_parser.Structure", _FakeStructure)
    def test_alexandria_query_maps_formation_energy_and_pagination(self):
        client = OptimadeClient(providers=["alexandria"], page_limit=1)
        payload_1 = {
            "data": [_optimade_item("agm1", _alexandria_formation_energy_per_atom=-0.25)],
            "links": {"next": "/pbe/v1/structures?page_offset=1"},
        }
        payload_2 = {
            "data": [_optimade_item("agm2", _alexandria_formation_energy_per_atom=-0.5)],
            "links": {},
        }

        with patch.object(OptimadeClient, "_get_json", autospec=True,
                          side_effect=[payload_1, payload_2]) as get_json:
            df = client.query({"Si"})

        self.assertEqual(len(df), 2)
        self.assertEqual(list(df["id"]), ["optimade_alexandria_agm1", "optimade_alexandria_agm2"])
        self.assertEqual(list(df["energy"]), [-0.25, -0.5])
        self.assertEqual(df.iloc[0]["metadata"]["energy_field"], "_alexandria_formation_energy_per_atom")
        self.assertEqual(get_json.call_count, 2)

    @patch("vsbtools.materials_tools.materials_dataset.io.sources.optimade_parser.Structure", _FakeStructure)
    def test_materials_project_query_maps_selected_thermo_type(self):
        client = OptimadeClient(providers=["mp"], mp_thermo_type="gga_gga+u_r2scan")
        payload = {
            "data": [
                _optimade_item(
                    "mp-1",
                    _mp_stability={
                        "gga_gga+u": {"formation_energy_per_atom": -0.1},
                        "gga_gga+u_r2scan": {"formation_energy_per_atom": -0.2},
                    },
                )
            ],
            "links": {},
        }

        with patch.object(OptimadeClient, "_get_json", autospec=True, return_value=payload):
            df = client.query({"Si"})

        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["energy"], -0.2)
        self.assertEqual(df.iloc[0]["metadata"]["mp_thermo_type"], "gga_gga+u_r2scan")

    @patch("vsbtools.materials_tools.materials_dataset.io.sources.optimade_parser.Structure", _FakeStructure)
    def test_materials_project_canonical_interface_maps_formation_energy(self):
        client = OptimadeClient(providers=["materials_project"])
        payload = {
            "data": [
                _optimade_item(
                    "mp-2",
                    _mp_stability={"gga_gga+u": {"formation_energy_per_atom": -0.125}},
                )
            ],
            "links": {},
        }

        with patch.object(OptimadeClient, "_get_json", autospec=True, return_value=payload) as get_json:
            df = client.query({"Si"})

        requested_url = get_json.call_args.args[1]
        self.assertIn("https://optimade.materialsproject.org/v1/structures", requested_url)
        self.assertIn("_mp_stability", requested_url)
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["id"], "optimade_materials_project_mp-2")
        self.assertEqual(df.iloc[0]["energy"], -0.125)
        self.assertEqual(df.iloc[0]["metadata"]["source"], "MaterialsProject/OPTIMADE")
        self.assertEqual(df.iloc[0]["metadata"]["energy_field"], "_mp_stability.formation_energy_per_atom")

    @patch("vsbtools.materials_tools.materials_dataset.io.sources.optimade_parser.Structure", _FakeStructure)
    def test_oqmd_interface_maps_delta_e_per_atom(self):
        client = OptimadeClient(providers=["oqmd"])
        payload = {
            "data": [
                _optimade_item(
                    "123",
                    cartesian_site_positions=[[0, 0, 0], [0.25, 0.25, 0.25]],
                    species_at_sites=["Si", "Si"],
                    _oqmd_delta_e=-0.3,
                )
            ],
            "links": {},
        }

        with patch.object(OptimadeClient, "_get_json", autospec=True, return_value=payload) as get_json:
            df = client.query({"Si"})

        requested_url = get_json.call_args.args[1]
        self.assertIn("https://oqmd.org/optimade/v1/structures", requested_url)
        self.assertIn("_oqmd_delta_e", requested_url)
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["id"], "optimade_oqmd_123")
        self.assertEqual(df.iloc[0]["energy"], -0.6)
        self.assertEqual(df.iloc[0]["metadata"]["source"], "OQMD/OPTIMADE")
        self.assertEqual(df.iloc[0]["metadata"]["energy_field"], "_oqmd_delta_e")

    @patch("vsbtools.materials_tools.materials_dataset.io.sources.optimade_parser.Structure", _FakeStructure)
    def test_total_energy_mode_keeps_only_providers_with_total_energy(self):
        client = OptimadeClient(providers=["alexandria"], energy_mode="total")
        payload = {
            "data": [_optimade_item("agm1", _alexandria_energy=-10.0)],
            "links": {},
        }

        with patch.object(OptimadeClient, "_get_json", autospec=True, return_value=payload):
            df = client.query({"Si"})

        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["energy"], -10.0)
        self.assertEqual(df.iloc[0]["metadata"]["energy_field"], "_alexandria_energy")

    @patch("vsbtools.materials_tools.materials_dataset.io.sources.optimade_parser.Structure", _FakeStructure)
    def test_all_provider_interfaces_are_queried(self):
        client = OptimadeClient(providers="all")
        payloads = [
            {
                "data": [
                    _optimade_item(
                        "mp-1",
                        _mp_stability={"gga_gga+u": {"formation_energy_per_atom": -0.1}},
                    )
                ],
                "links": {},
            },
            {
                "data": [_optimade_item("oqmd-1", _oqmd_delta_e=-0.2)],
                "links": {},
            },
            {
                "data": [_optimade_item("agm1", _alexandria_formation_energy_per_atom=-0.3)],
                "links": {},
            },
        ]

        with patch.object(OptimadeClient, "_get_json", autospec=True,
                          side_effect=payloads) as get_json:
            df = client.query({"Si"})

        self.assertEqual(get_json.call_count, 3)
        self.assertEqual(
            set(df["id"]),
            {
                "optimade_materials_project_mp-1",
                "optimade_oqmd_oqmd-1",
                "optimade_alexandria_agm1",
            },
        )
        self.assertEqual(
            set(df["metadata"].map(lambda md: md["optimade_provider"])),
            {"materials_project", "oqmd", "alexandria"},
        )


if __name__ == "__main__":
    unittest.main()
