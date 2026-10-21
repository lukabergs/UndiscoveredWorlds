"""No function of the run tools makes a module it imports a local name (master plan 2.4 review).

experiment.py's main() once assigned `inputs = [...]` at its end, which made the imported `inputs` module a local of
the whole function: `inputs.resolve` near its start then raised UnboundLocalError for every profile with inputs
(climate, hydrology). symtable finds such a shadow without running anything.
"""

import importlib
from pathlib import Path
import symtable
import sys
import unittest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))

CHECKED = ("evaluation/experiment.py", "evaluation/step_check.py", "evaluation/acceptance.py", "uw.py",
           "uwtools/naming.py")


def shadows(path):
    """(function, name) of every function-local name that is also a name the module imports."""
    table = symtable.symtable(path.read_text(encoding="utf-8"), str(path), "exec")
    imported = {symbol.get_name() for symbol in table.get_symbols() if symbol.is_imported()}
    found, stack = [], list(table.get_children())
    while stack:
        child = stack.pop()
        if child.get_type() == "function":
            found += [(child.get_name(), symbol.get_name()) for symbol in child.get_symbols()
                      if symbol.get_name() in imported and symbol.is_local() and not symbol.is_global()]
        stack += child.get_children()
    return found


class ToolScopeTests(unittest.TestCase):
    def test_experiment_main_reads_the_inputs_module(self):
        experiment = importlib.import_module("experiment")
        self.assertTrue(hasattr(experiment, "inputs"))
        table = symtable.symtable(Path(experiment.__file__).read_text(encoding="utf-8"), experiment.__file__, "exec")
        main = next(child for child in table.get_children() if child.get_name() == "main")
        symbol = main.lookup("inputs")
        self.assertFalse(symbol.is_local(), "experiment.main() assigns `inputs`, shadowing the inputs module")
        self.assertTrue(symbol.is_global())

    def test_no_function_shadows_an_imported_module(self):
        for name in CHECKED:
            with self.subTest(script=name):
                self.assertEqual(shadows(SCRIPTS / name), [])


if __name__ == "__main__":
    unittest.main()
