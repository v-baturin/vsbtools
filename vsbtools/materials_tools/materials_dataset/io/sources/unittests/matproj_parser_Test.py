import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ..matproj_parser import MPClient, _load_saved_key, _save_key


class MPClientAuthenticationTest(unittest.TestCase):
    def test_saved_key_is_used_without_prompting(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            key_path = Path(tmp_dir) / "materials_project.json"
            _save_key(key_path, "saved-key")
            client = MPClient(key_path=key_path, prompt_for_key=False)

            self.assertEqual(client._resolve_api_key(), "saved-key")
            self.assertEqual(_load_saved_key(key_path), "saved-key")

    def test_first_interactive_key_is_persisted(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            key_path = Path(tmp_dir) / "materials_project.json"
            client = MPClient(key_path=key_path)

            with patch("vsbtools.materials_tools.materials_dataset.io.sources.matproj_parser.getpass.getpass", return_value="new-key"):
                self.assertEqual(client._resolve_api_key(), "new-key")

            self.assertEqual(_load_saved_key(key_path), "new-key")

    def test_authentication_failure_requests_one_replacement_key(self):
        client = MPClient(prompt_for_key=False)
        client._mpr = object()

        with patch.object(client, "_prompt_for_key", return_value="replacement-key") as prompt:
            self.assertTrue(client._retry_after_authentication_failure(RuntimeError("HTTP 401 Unauthorized")))

        prompt.assert_called_once()
        self.assertIsNone(client._mpr)
        self.assertEqual(client.api_key, "replacement-key")

    def test_non_authentication_failure_is_not_retried(self):
        client = MPClient(prompt_for_key=False)

        self.assertFalse(client._retry_after_authentication_failure(RuntimeError("Connection timed out")))
