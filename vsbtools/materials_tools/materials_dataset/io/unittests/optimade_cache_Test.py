import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock

from ..preset_loaders import cache_loader


class OptimadeCacheTest(unittest.TestCase):

    def test_cache_is_created_lazily_with_safe_bounded_name(self):
        read_fn = Mock(return_value="cached")

        def write_fn(dataset, *, enforce_base_path, **_):
            (enforce_base_path / "manifest.yaml").touch()

        with TemporaryDirectory() as tmp:
            cache_root = Path(tmp) / "not-created-yet"

            @cache_loader(read_fn, write_fn, cache_root=cache_root)
            def load_from_optimade(elements, message=None, **kwargs):
                return "fresh"

            self.assertFalse(cache_root.exists())
            result = load_from_optimade(
                {"Si"},
                cache_label='CON:<bad>|name?*' + "x" * 200,
                do_deduplication=False,
            )

            self.assertEqual(result, "fresh")
            folders = list(cache_root.iterdir())
            self.assertEqual(len(folders), 1)
            self.assertLessEqual(len(folders[0].name), 120)
            self.assertFalse(set('<>:"/\\|?*') & set(folders[0].name))

            self.assertEqual(
                load_from_optimade(
                    {"Si"},
                    cache_label='CON:<bad>|name?*' + "x" * 200,
                    do_deduplication=False,
                ),
                "cached",
            )
            read_fn.assert_called_once()


if __name__ == "__main__":
    unittest.main()
