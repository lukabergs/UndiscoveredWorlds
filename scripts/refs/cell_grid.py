"""Moved to scripts/uwtools/cell_grid.py. This name stays for the scripts that import it from here."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import cell_grid

sys.modules[__name__] = cell_grid
