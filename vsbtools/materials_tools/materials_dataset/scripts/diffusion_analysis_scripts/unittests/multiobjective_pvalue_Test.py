import unittest
from pathlib import Path

import yaml

from ..pvalue_utils import get_p_value, get_two_proportion_z_test


def _build_callables_from_objectives(objectives_yml: Path):
    from ..guidance_stats import guidance_vs_target_properties
    from ..mattergen_bridge import get_target_value_fn

    with open(objectives_yml, "rt") as f:
        objectives = yaml.safe_load(f)
    callables = {}
    targets = {}

    for objective_name, objective_desc in objectives.items():
        guidance_name = objective_desc["guidance_name"]
        fn_name = guidance_vs_target_properties[guidance_name]
        fn_parameters = objective_desc.get("fn_parameters", {})
        callables[objective_name] = get_target_value_fn(fn_name, **fn_parameters)
        targets[objective_name] = float(objective_desc["target"])

    return callables, targets


class MultiobjectivePValue_Test(unittest.TestCase):
    def test_guided_pvalue_lt_0p2(self):
        try:
            import pymatgen  # noqa: F401
        except ModuleNotFoundError as e:
            raise unittest.SkipTest("pymatgen not installed") from e

        from ....io.structures_dataset_io import StructureDatasetIO

        root = Path(__file__).resolve().parents[3]  # .../materials_dataset
        multiobjective_root = root / "unittests_datasets" / "multiobjective"
        objectives_yml = multiobjective_root / "objectives.yml"
        guided_extxyz = next(p for p in multiobjective_root.glob("*.extxyz") if "nonguided" not in p.name)
        nonguided_extxyz = next(p for p in multiobjective_root.glob("*nonguided*.extxyz"))

        utils = StructureDatasetIO(multiobjective_root)
        ds_guided = utils.load_from_extxyz(guided_extxyz)
        ds_non_guided = utils.load_from_extxyz(nonguided_extxyz)

        try:
            callables, targets = _build_callables_from_objectives(objectives_yml)
        except ModuleNotFoundError as e:
            raise unittest.SkipTest(f"Mattergen not configured: {e}") from e
        margins = {k: 0.5 for k in callables.keys()}

        def is_good(entry):
            for key, fn in callables.items():
                if abs(fn(entry) - targets[key]) > margins[key]:
                    return False
            return True

        good_guided = [e for e in ds_guided if is_good(e)]
        good_non_guided = [e for e in ds_non_guided if is_good(e)]

        share_guided = len(good_guided) / len(ds_guided)
        share_non_guided = len(good_non_guided) / len(ds_non_guided)
        self.assertGreaterEqual(share_guided, share_non_guided)
        p_value_exact = get_p_value(
            ds_non_guided,
            callables=callables,
            targets=targets,
            margins=margins,
            ds_guided=ds_guided,
            alternative="greater",
            method="exact",
            continuity_correction=True,
        )
        z_stats = get_two_proportion_z_test(
            ds_non_guided,
            callables=callables,
            targets=targets,
            margins=margins,
            ds_guided=ds_guided,
            alternative="greater",
        )

        print("\n=== Multiobjective Guidance Stats ===")
        print(
            f"guided:     {len(good_guided):5d} / {len(ds_guided):5d} "
            f"({share_guided:7.2%})"
        )
        print(
            f"non-guided: {len(good_non_guided):5d} / {len(ds_non_guided):5d} "
            f"({share_non_guided:7.2%})"
        )
        print(f"exact p-value:          {p_value_exact: .6e}")
        print(f"two-proportion z-score: {z_stats['z_score']: .6f}")
        print(f"two-proportion p-value: {z_stats['p_value']: .6e}")

        self.assertLess(p_value_exact, 0.2)
        self.assertLess(z_stats["p_value"], 0.2)
