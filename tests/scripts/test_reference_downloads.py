"""Deterministic HTTP resume checks; no network requests."""

import importlib.util
from pathlib import Path
import tempfile
import sys
import os
from contextlib import nullcontext
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

try:
    import requests
except ImportError:
    requests = None


@unittest.skipIf(requests is None, "Run with --with requests")
class ReferenceDownloadTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[2] / "scripts/refs/download-reference-sources.py"
        spec = importlib.util.spec_from_file_location("reference_download", path)
        cls.download = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.download)

    def response(self, status, headers, blocks):
        response = Mock(status_code=status, headers=headers)
        response.__enter__ = Mock(return_value=response)
        response.__exit__ = Mock(return_value=False)
        response.iter_content.return_value = iter(blocks)
        return response

    def test_resume_requires_matching_strong_etag(self):
        with tempfile.TemporaryDirectory() as temp:
            partial = Path(temp) / "source.part"
            partial.write_bytes(b"abc")
            state = partial.with_suffix(".part.http.json")
            self.download.write_json(state, dict(request=dict(url="https://example.org/data", params={}), etag='"v1"'))
            response = self.response(206, {"ETag": '"v1"', "Content-Range": "bytes 3-5/6"}, [b"def"])
            with patch.object(requests, "get", return_value=response) as get:
                self.download.http_download("https://example.org/data", {}, partial)
            self.assertEqual(partial.read_bytes(), b"abcdef")
            self.assertEqual(get.call_args.kwargs["headers"], {"Range": "bytes=3-", "If-Range": '"v1"'})

    def test_unidentified_partial_restarts_and_changed_etag_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            partial = Path(temp) / "source.part"
            partial.write_bytes(b"old")
            response = self.response(200, {}, [b"new"])
            with patch.object(requests, "get", return_value=response) as get:
                self.download.http_download("https://example.org/data", {}, partial)
            self.assertEqual(partial.read_bytes(), b"new")
            self.assertEqual(get.call_args.kwargs["headers"], {})
            state = partial.with_suffix(".part.http.json")
            self.download.write_json(state, dict(request=dict(url="https://example.org/data", params={}), etag='"v1"'))
            response = self.response(206, {"ETag": '"v2"', "Content-Range": "bytes 3-5/6"}, [b"bad"])
            with patch.object(requests, "get", return_value=response):
                with self.assertRaisesRegex(ValueError, "source changed"):
                    self.download.http_download("https://example.org/data", {}, partial)
            self.assertEqual(partial.read_bytes(), b"new")

    def test_download_preparation_runs_ingestion_before_family_refresh(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            script = root / "scripts/refs/download-reference-sources.py"
            argv = [str(script), "--source", "moisture", "--output", str(root / "sources"), "--prepare"]
            with patch.object(sys, "argv", argv), patch.object(self.download, "__file__", str(script)), \
                 patch.object(self.download, "cds_download") as download, \
                 patch.dict(os.environ, {"UV": str(root / "uv.exe")}), \
                 patch.dict(sys.modules, {"filelock": SimpleNamespace(FileLock=lambda path: nullcontext())}), \
                 patch("subprocess.run") as run:
                self.download.main()
            download.assert_called_once_with("moisture", root / "sources", 2001, 2020)
            self.assertEqual(run.call_count, 2)
            ingest, refresh = [call.args[0] for call in run.call_args_list]
            self.assertEqual(ingest[0], str(root / "uv.exe"))
            self.assertEqual(refresh[0], str(root / "uv.exe"))
            self.assertIn("scripts/refs/prepare-additional-references.py", ingest)
            self.assertIn(str(root / "sources/era5-cds-2001-2020"), ingest)
            self.assertEqual(refresh[-2:], ["--additional-family", "era5"])
            self.assertTrue(all(call.kwargs["check"] and call.kwargs["cwd"] == root for call in run.call_args_list))


if __name__ == "__main__":
    unittest.main()
