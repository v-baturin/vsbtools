import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ..external_paths import import_from_path_validator, resolve_external_path


def existing_directory_validator(path: Path):
    return (path.is_dir(), f"{path} is not a directory")


class ExternalPaths_Test(unittest.TestCase):

    def test_resolves_saved_config_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            configured_path = root / "resource"
            configured_path.mkdir()
            config_path = root / "external_paths.json"
            config_path.write_text(json.dumps({"resource": configured_path.as_posix()}), encoding="utf-8")

            with patch.dict(os.environ, {"VSBTOOLS_EXTERNAL_PATHS_CONFIG": config_path.as_posix()}):
                resolved = resolve_external_path(
                    name="test resource",
                    config_key="resource",
                    validator=existing_directory_validator,
                    prompt=False,
                )

            self.assertEqual(resolved, configured_path)

    def test_prompts_until_valid_and_saves_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            config_path = root / "external_paths.json"
            valid_path = root / "valid"
            valid_path.mkdir()

            with patch.dict(os.environ, {"VSBTOOLS_EXTERNAL_PATHS_CONFIG": config_path.as_posix()}):
                with patch("builtins.input", return_value=valid_path.as_posix()):
                    resolved = resolve_external_path(
                        name="test resource",
                        config_key="resource",
                        explicit_path=root / "missing",
                        validator=existing_directory_validator,
                    )

            self.assertEqual(resolved, valid_path)
            self.assertEqual(json.loads(config_path.read_text(encoding="utf-8"))["resource"], valid_path.as_posix())

    def test_resolves_hostname_preset(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            preset_path = root / "external_paths_presets.json"
            valid_path = root / "valid"
            valid_path.mkdir()
            preset_path.write_text(
                json.dumps({"testhost": {"resource": valid_path.as_posix()}}),
                encoding="utf-8",
            )

            with patch.dict(os.environ, {
                "VSBTOOLS_EXTERNAL_PATHS_CONFIG": (root / "missing_config.json").as_posix(),
                "VSBTOOLS_EXTERNAL_PATHS_PRESETS": preset_path.as_posix(),
            }):
                with patch("socket.gethostname", return_value="testhost"):
                    resolved = resolve_external_path(
                        name="test resource",
                        config_key="resource",
                        validator=existing_directory_validator,
                        prompt=False,
                    )

            self.assertEqual(resolved, valid_path)

    def test_invalid_prompt_can_quit(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            config_path = root / "external_paths.json"

            with patch.dict(os.environ, {"VSBTOOLS_EXTERNAL_PATHS_CONFIG": config_path.as_posix()}):
                with patch("builtins.input", side_effect=[(root / "missing").as_posix(), "q"]):
                    with self.assertRaises(SystemExit):
                        resolve_external_path(
                            name="test resource",
                            config_key="resource",
                            validator=existing_directory_validator,
                        )

    def test_import_from_path_validator_uses_candidate_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            package = root / "candidate_package"
            package.mkdir()
            (package / "__init__.py").write_text("VALUE = 1\n", encoding="utf-8")

            ok, message = import_from_path_validator("candidate_package")(root)

            self.assertTrue(ok, message)
