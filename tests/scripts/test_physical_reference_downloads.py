"""Deterministic download safeguards; no network or real data required."""

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import MagicMock, patch

SPEC = importlib.util.spec_from_file_location('physical_download',
    Path(__file__).resolve().parents[2] / 'scripts/refs/download-physical-reference-sources.py')
DOWNLOAD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DOWNLOAD)


class PhysicalDownloadTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.entry = dict(dataset='example', file='data.json', kind='json',
                          url='https://example.org/data', max_mib=1)

    def response(self, body):
        response = MagicMock(status_code=200, headers={})
        response.__enter__.return_value = response
        response.url = 'https://example.org/data?secret=temporary-token'
        response.iter_content.return_value = iter([body])
        client = MagicMock()
        client.get.return_value = response
        return client

    def fetch(self, verify_only=False, budget=1024**2):
        return DOWNLOAD.fetch(self.entry, self.root, budget, 0, verify_only)

    def test_successful_download_is_reusable_and_receipt_redacts_redirect_query(self):
        with patch.object(DOWNLOAD, 'session', return_value=self.response(b'{"value":1}')):
            self.fetch()
        receipt = json.loads((self.root / 'example/data.json.receipt.json').read_text())
        self.assertEqual(receipt['final_url'], 'https://example.org/data')
        with patch.object(DOWNLOAD, 'session') as session:
            self.assertEqual(self.fetch(verify_only=True), 11)
            session.assert_not_called()
        (self.root / 'example/data.json').write_bytes(b'{"value":2}')
        with self.assertRaisesRegex(ValueError, 'Cached source'):
            self.fetch(verify_only=True)

    def test_stream_over_budget_is_not_promoted(self):
        with patch.object(DOWNLOAD, 'session', return_value=self.response(b'123456789')):
            with self.assertRaisesRegex(ValueError, 'storage limit'):
                self.fetch(budget=8)
        self.assertFalse(list(self.root.rglob('*.part')))
        self.assertFalse((self.root / 'example/data.json').exists())

    def test_login_html_is_not_saved_as_data(self):
        with patch.object(DOWNLOAD, 'session', return_value=self.response(b'<html>Sign in</html>')):
            with self.assertRaisesRegex(ValueError, 'HTML page'):
                self.fetch()
        self.assertFalse(list(self.root.rglob('*.part')))
        self.assertFalse((self.root / 'example/data.json').exists())

    def test_catalog_path_cannot_escape_output(self):
        self.entry['dataset'] = '..'
        with patch.object(DOWNLOAD, 'session') as session:
            with self.assertRaisesRegex(ValueError, 'escapes output'):
                self.fetch()
            session.assert_not_called()

    def test_existing_unreceipted_file_is_preserved(self):
        path = self.root / 'example/data.json'
        path.parent.mkdir()
        path.write_bytes(b'original')
        with self.assertRaisesRegex(ValueError, 'no receipt'):
            self.fetch()
        self.assertEqual(path.read_bytes(), b'original')


if __name__ == '__main__':
    unittest.main()
