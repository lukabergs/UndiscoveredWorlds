"""Exercise repository scope against real staged files and deleted history."""
import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('repository_scope', ROOT / 'scripts/repository_scope.py')
scope = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scope)


class RepositoryScopeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='uw-scope-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.env = dict(os.environ, GIT_CONFIG_NOSYSTEM='1',
                        GIT_CONFIG_GLOBAL=os.devnull,
                        GIT_AUTHOR_NAME='Scope test', GIT_AUTHOR_EMAIL='scope@example.invalid',
                        GIT_COMMITTER_NAME='Scope test', GIT_COMMITTER_EMAIL='scope@example.invalid',
                        GIT_AUTHOR_DATE='2026-01-01T00:00:00Z',
                        GIT_COMMITTER_DATE='2026-01-01T00:00:00Z')
        self.git('init', '--quiet')

    def git(self, *args):
        return subprocess.check_output(['git', *args], cwd=self.root, env=self.env,
                                       stderr=subprocess.STDOUT)

    def stage(self, path, data=b'fixture\n'):
        target = self.root / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        self.git('add', '--', path)

    def test_code_and_small_fixture_pass(self):
        self.stage('scripts/prepare.py')
        self.stage('tests/climate/fixtures/tiny.npy')
        issues, count = scope.inspect(self.root)
        self.assertEqual(issues, [])
        self.assertEqual(count, 2)
        self.git('commit', '--quiet', '-m', 'Initial code and fixture')
        self.assertEqual(scope.inspect(self.root, ['HEAD'])[0], [])

    def test_local_directories_and_archives_are_rejected(self):
        for path in ['meta/readme.md', 'docs/design.md', 'refs/prepare.py',
                     'runs/registry/climate.json', 'assets/download.zip']:
            self.stage(path)
        issues, _ = scope.inspect(self.root)
        self.assertEqual(len(issues), 5)

    def test_deleted_data_is_still_rejected_in_history(self):
        self.stage('src/main.cpp')
        self.stage('refs/source/input.nc')
        self.git('commit', '--quiet', '-m', 'Accidental data')
        self.git('rm', '--', 'refs/source/input.nc')
        self.git('commit', '--quiet', '-m', 'Remove data from current tree')
        self.assertEqual(scope.inspect(self.root)[0], [])
        issues, _ = scope.inspect(self.root, ['HEAD'])
        self.assertTrue(any('refs/source/input.nc' in issue for issue in issues))

    def test_large_blob_cannot_hide_in_fixture_directory(self):
        self.stage('tests/climate/fixtures/large.npy', b'0' * (scope.MAX_BYTES + 1))
        issues, _ = scope.inspect(self.root)
        self.assertEqual(len(issues), 1)
        self.assertIn('exceeds', issues[0])

    def test_pre_push_rejects_data_history(self):
        self.stage('runs/result.txt')
        self.git('commit', '--quiet', '-m', 'Data')
        oid = self.git('rev-parse', 'HEAD').decode().strip()
        record = f'refs/heads/main {oid} refs/heads/main {"0" * 40}\n'
        result = subprocess.run([os.sys.executable, str(ROOT / 'scripts/repository_scope.py'),
                                 '--root', str(self.root), '--pre-push'],
                                input=record, text=True, capture_output=True, env=self.env)
        self.assertEqual(result.returncode, 1)
        self.assertIn('runs/result.txt', result.stderr)


if __name__ == '__main__':
    unittest.main()
