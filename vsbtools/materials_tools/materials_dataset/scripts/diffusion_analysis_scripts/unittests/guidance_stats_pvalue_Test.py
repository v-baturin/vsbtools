import unittest

from ..pvalue_utils import get_p_value, get_two_proportion_z_test


class GuidanceStatsPValue_Test(unittest.TestCase):
    @staticmethod
    def _binary_case(good_ng=1040, tot_ng=9300, good_g=446, tot_g=3600):
        ds_non_guided = [1] * good_ng + [0] * (tot_ng - good_ng)
        ds_guided = [1] * good_g + [0] * (tot_g - good_g)
        callables = {"x": lambda e: e}
        targets = {"x": 1.0}
        margins = {"x": 0.0}
        return ds_non_guided, ds_guided, callables, targets, margins

    def test_normal_matches_two_proportion_api(self):
        ds_non_guided, ds_guided, callables, targets, margins = self._binary_case()

        p_norm = get_p_value(
            ds_non_guided,
            callables=callables,
            targets=targets,
            margins=margins,
            ds_guided=ds_guided,
            alternative="greater",
            method="normal",
        )

        stats = get_two_proportion_z_test(
            ds_non_guided,
            callables=callables,
            targets=targets,
            margins=margins,
            ds_guided=ds_guided,
            alternative="greater",
        )
        k0 = sum(ds_non_guided)
        n0 = len(ds_non_guided)
        k1 = sum(ds_guided)
        n1 = len(ds_guided)
        print("\n=== Pooled Z-Test Consistency ===")
        print(f"guided:     {k1:5d} / {n1:5d} ({k1 / n1:7.2%})")
        print(f"non-guided: {k0:5d} / {n0:5d} ({k0 / n0:7.2%})")
        print(f"normal p-value:         {p_norm: .6e}")
        print(f"two-proportion z-score: {stats['z_score']: .6f}")
        print(f"two-proportion p-value: {stats['p_value']: .6e}")
        self.assertAlmostEqual(p_norm, stats["p_value"], places=12)

    def test_two_proportion_z_test_returns_z_and_p(self):
        ds_non_guided, ds_guided, callables, targets, margins = self._binary_case()

        stats = get_two_proportion_z_test(
            ds_non_guided,
            callables=callables,
            targets=targets,
            margins=margins,
            ds_guided=ds_guided,
            alternative="greater",
        )
        self.assertIn("z_score", stats)
        self.assertIn("p_value", stats)
        self.assertGreater(stats["z_score"], 0.0)

    def test_normal_tail_direction(self):
        ds_non_guided, ds_guided, callables, targets, margins = self._binary_case()

        p_greater = get_p_value(
            ds_non_guided,
            callables=callables,
            targets=targets,
            margins=margins,
            ds_guided=ds_guided,
            alternative="greater",
            method="normal",
        )
        p_less = get_p_value(
            ds_non_guided,
            callables=callables,
            targets=targets,
            margins=margins,
            ds_guided=ds_guided,
            alternative="less",
            method="normal",
        )

        self.assertLess(p_greater, 0.05)
        self.assertGreater(p_less, 0.95)

    def test_exact_no_overflow_large_n(self):
        ds_non_guided = [1] * 1800 + [0] * (6000 - 1800)
        ds_guided = [1] * 1700 + [0] * (5000 - 1700)
        callables = {"x": lambda e: e}
        targets = {"x": 1.0}
        margins = {"x": 0.0}

        p_exact = get_p_value(
            ds_non_guided,
            callables=callables,
            targets=targets,
            margins=margins,
            ds_guided=ds_guided,
            alternative="greater",
            method="exact",
        )
        self.assertGreaterEqual(p_exact, 0.0)
        self.assertLessEqual(p_exact, 1.0)

    def test_invalid_method_raises(self):
        ds_non_guided = [1, 0, 0, 1]
        ds_guided = [1, 0]
        callables = {"x": lambda e: e}
        targets = {"x": 1.0}
        margins = {"x": 0.0}

        with self.assertRaises(ValueError):
            get_p_value(
                ds_non_guided,
                callables=callables,
                targets=targets,
                margins=margins,
                ds_guided=ds_guided,
                method="nope",
            )
