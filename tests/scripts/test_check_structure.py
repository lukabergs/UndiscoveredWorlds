"""Path rules of scripts/check_structure.py with the tracked configs/structure.json (no filesystem writes)."""
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('check_structure', ROOT / 'scripts/check_structure.py')
structure = importlib.util.module_from_spec(spec)
spec.loader.exec_module(structure)
CONFIG = structure.load(ROOT / 'configs/structure.json')


class StructureRuleTests(unittest.TestCase):
    def assertAllowed(self, *paths):
        for path in paths:
            with self.subTest(path=path):
                self.assertIsNone(structure.judge(path, CONFIG))

    def assertRefused(self, *paths):
        for path in paths:
            with self.subTest(path=path):
                self.assertIsNotNone(structure.judge(path, CONFIG))

    def test_named_places_and_kinds_pass(self):
        self.assertAllowed('AGENTS.md', 'pyproject.toml', '.claude/settings.json', '.claude/agents/worker.md',
                           'src/simulations/terrain/terrain.cpp', 'src/simulations/new_module/CMakeLists.txt',
                           'src/core/mesh/cell_grid.hpp', 'tests/simulations/landscape/x_tests.cpp',
                           'tests/scripts/test_x.py', 'tests/io/fixtures/tiny.npy', 'scripts/uwtools/x.py',
                           'docs/tasks/new-topic/notes.md', 'docs/architecture/tectonics/x.md',
                           'configs/parameters/terrain.json', 'runs/registry/world.json',
                           'runs/ledger/world.jsonl', 'meta/anything/at all.bin', 'out/scratch/W2-W2f/probe.txt',
                           'refs/processed/climate/maps/x.png', 'runs/world/12/run.json',
                           'runs/world/12/maps/relief.png')

    def test_unnamed_folders_are_refused(self):
        self.assertRefused('tmp/x.txt', 'notes.txt', 'docs/notes/x.md', 'scripts/helpers/x.py',
                           'src/simulations/New Folder/x.cpp', 'src/gui/x.cpp', 'configs/old/x.json',
                           'out/images/x.png', 'refs/extra/x.csv', 'runs/world/latest/run.json',
                           'runs/world/run.json')
        self.assertIn('not a place', structure.judge('tmp/x.txt', CONFIG))

    def test_wrong_kinds_are_refused(self):
        self.assertRefused('scripts/sweep.ps1', 'src/core/x.py', 'configs/x.yaml', 'docs/architecture/x.txt',
                           'tests/scripts/data.json', '.claude/agents/x.json')
        self.assertIn('kind', structure.judge('scripts/sweep.ps1', CONFIG))

    def test_html_only_where_pages_are_rendered(self):
        self.assertAllowed('scripts/evaluation/templates/run.html', 'out/evaluation/climate/x/index.html',
                           'runs/reports/163.html', 'runs/world/1/evaluation/index.html',
                           'runs/world/1/evaluation-2/index.html', 'scripts/refs/physical_map_index.html')
        self.assertRefused('out/scratch/sandbox/compare.html', 'scripts/evaluation/page.html',
                           'runs/world/1/maps/index.html', 'docs/reference-guides/page.htm',
                           'meta/page.html', 'scripts/refs/other.html')

    def test_runs_hold_no_scripts_and_archive_is_sealed(self):
        self.assertRefused('runs/reports/163-checks.py', 'runs/world/1/render.py', 'runs/work/x.ps1',
                           'runs/archive/climate/1/run.json')
        self.assertIsNotNone(structure.judge('runs/archive/x.json', CONFIG, existing=True))
        self.assertIsNone(structure.judge('out/images/x.png', CONFIG, existing=True))

    def test_worktree_paths_are_judged_from_their_checkout(self):
        self.assertAllowed('.claude/worktrees/agent-a1/src/core/x.cpp', '.claude/worktrees/agent-a1/AGENTS.md')
        self.assertRefused('.claude/worktrees/agent-a1/tmp/x.txt', '.claude/worktrees/agent-a1/runs/reports/x.py')

    def test_repository_path(self):
        self.assertEqual(structure.repository_path(str(ROOT / 'src/core/x.cpp'), ROOT), 'src/core/x.cpp')
        self.assertEqual(structure.repository_path('src/core/x.cpp', ROOT, cwd=str(ROOT)), 'src/core/x.cpp')
        self.assertIsNone(structure.repository_path(str(ROOT.parent / 'elsewhere/x.txt'), ROOT))
        self.assertIsNone(structure.repository_path(str(ROOT) + '-sibling/x.txt', ROOT))
        self.assertIsNone(structure.repository_path(str(ROOT), ROOT))

    def test_hook_decision(self):
        def event(path, tool='Write'):
            return {'tool_name': tool, 'tool_input': {'file_path': path}, 'cwd': str(ROOT)}
        self.assertIsNone(structure.hook_decision(event(str(ROOT.parent / 'elsewhere/x.txt')), CONFIG, ROOT))
        self.assertIsNone(structure.hook_decision(event(str(ROOT / 'out/scratch/topic/x.txt')), CONFIG, ROOT))
        self.assertIsNone(structure.hook_decision({'tool_input': {}}, CONFIG, ROOT))
        # An existing file in a named place, edited again.
        self.assertIsNone(structure.hook_decision(event(str(ROOT / 'AGENTS.md'), 'Edit'), CONFIG, ROOT))
        message = structure.hook_decision(event(str(ROOT / 'no-such-folder-x/new.txt')), CONFIG, ROOT)
        self.assertTrue(message.startswith('BLOCKED: no-such-folder-x/new.txt'))
        notebook = {'tool_input': {'notebook_path': str(ROOT / 'no-such-folder-x/n.ipynb')}, 'cwd': str(ROOT)}
        self.assertIsNotNone(structure.hook_decision(notebook, CONFIG, ROOT))

    def test_full_check_lists_each_violation(self):
        issues = structure.full_issues(['AGENTS.md', 'tmp/x.txt', 'scripts/x.ps1'], CONFIG)
        self.assertEqual([issue.split(':')[0] for issue in issues], ['tmp/x.txt', 'scripts/x.ps1'])


if __name__ == '__main__':
    unittest.main()
