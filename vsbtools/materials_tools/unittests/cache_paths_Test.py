import unittest
from pathlib import Path, PureWindowsPath
from tempfile import TemporaryDirectory

from ..cache_paths import CACHE_ENV_VAR, database_cache_dir, safe_cache_component, user_cache_dir


class CachePathsTest(unittest.TestCase):

    def test_windows_uses_local_app_data(self):
        path = user_cache_dir(
            environ={"LOCALAPPDATA": r"C:\Users\tester\AppData\Local"},
            platform="win32",
            home=PureWindowsPath(r"C:\Users\tester"),
        )
        self.assertEqual(
            PureWindowsPath(path),
            PureWindowsPath(r"C:\Users\tester\AppData\Local\vsbtools\Cache"),
        )

    def test_windows_falls_back_below_home(self):
        path = user_cache_dir(environ={}, platform="win32", home=Path("home"))
        self.assertEqual(path, Path("home") / "AppData" / "Local" / "vsbtools" / "Cache")

    def test_environment_override_has_precedence(self):
        path = database_cache_dir(
            environ={CACHE_ENV_VAR: "custom-cache"},
            platform="win32",
            home=Path("ignored"),
        )
        self.assertEqual(path, Path("custom-cache") / "DB_caches")

    def test_resolving_cache_does_not_create_it(self):
        with TemporaryDirectory() as tmp:
            path = database_cache_dir(
                environ={},
                platform="linux",
                home=Path(tmp),
            )
            self.assertFalse(path.exists())

    def test_cache_component_is_windows_safe_and_bounded(self):
        component = safe_cache_component('CON:<bad>|name?*' + "x" * 200)
        self.assertLessEqual(len(component), 80)
        self.assertFalse(set('<>:"/\\|?*') & set(component))

    def test_reserved_windows_name_is_prefixed(self):
        self.assertEqual(safe_cache_component("NUL"), "_NUL")


if __name__ == "__main__":
    unittest.main()
