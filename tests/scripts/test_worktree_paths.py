"""Scripts find the main tree from a git worktree, and check_layout accepts a worktree's .git file."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import check_layout
from uwtools import paths


def git(cwd, *arguments):
    subprocess.run(["git", *arguments], cwd=cwd, check=True, capture_output=True, text=True)


class WorktreePathTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        base = Path(self.temporary.name).resolve()
        self.main = base / "main"
        self.main.mkdir()
        git(self.main, "init", "-q")
        git(self.main, "-c", "user.name=test", "-c", "user.email=test@example.invalid",
            "commit", "-q", "--allow-empty", "-m", "root")
        self.worktree = self.main / ".claude" / "worktrees" / "agent"
        git(self.main, "worktree", "add", "-q", "--detach", str(self.worktree))

    def tearDown(self):
        self.temporary.cleanup()

    def test_main_root_in_main_tree_is_the_tree(self):
        self.assertEqual(paths.main_root(self.main), self.main)

    def test_main_root_in_worktree_is_the_main_tree(self):
        self.assertTrue((self.worktree / ".git").is_file())
        self.assertEqual(paths.main_root(self.worktree), self.main)

    def test_layout_accepts_a_worktree_git_file(self):
        self.assertEqual(check_layout.root_item_issues(self.worktree), [])
        self.assertEqual(check_layout.root_item_issues(self.main), [])
        (self.worktree / "stray.txt").write_text("", encoding="utf-8")
        self.assertEqual(len(check_layout.root_item_issues(self.worktree)), 1)


if __name__ == "__main__":
    unittest.main()
