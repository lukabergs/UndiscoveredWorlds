"""global-balanced-v1 for the legacy climate report: re-exported from scripts/evaluation/score.py."""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1] / 'evaluation'))
from score import DEFINITION, VERSION, bins, evaluate, quantile, rank  # noqa: E402,F401
