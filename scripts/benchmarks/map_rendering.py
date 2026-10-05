"""Moved to scripts/uwtools/map_rendering.py. This name stays for the scripts that import it from here."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import map_rendering

sys.modules[__name__] = map_rendering
