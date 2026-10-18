"""Build output locations: the one place scripts take executable paths from.

Presets (CMakePresets.json) configure into out/build/<preset>; binaries land in
bin/<configuration>. Scripts run the headless preset's Release binaries.
"""
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
PRESET = 'headless'
CONFIGURATION = 'Release'


def build_directory(preset=PRESET):
    return ROOT / 'out/build' / preset


def binary_directory(preset=PRESET, configuration=CONFIGURATION):
    return build_directory(preset) / 'bin' / configuration


def executable(name, preset=PRESET, configuration=CONFIGURATION):
    return binary_directory(preset, configuration) / f'{name}.exe'


def main_root(root=ROOT):
    """The main tree's root: root itself in the main tree, the main checkout in a git worktree.

    Agent worktrees hold no refs/ or runs/; scripts read those from the main tree.
    """
    common = subprocess.run(['git', 'rev-parse', '--git-common-dir'], cwd=root,
                            capture_output=True, text=True, check=True).stdout.strip()
    common = Path(common)
    if not common.is_absolute():
        common = Path(root) / common
    return common.resolve().parent
