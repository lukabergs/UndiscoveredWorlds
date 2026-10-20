"""Run a PreToolUse guard so that every failure blocks: `python hook_runner.py GUARD.py`.

Claude Code blocks a tool call only on exit code 2; any other failure of a hook (a syntax error, an import error,
an exception, a missing file) lets the command run. This runner passes the guard's own exit 0 (allow) or 2 (block)
through and turns everything else into 2. Keep it this small: it is the one piece that cannot protect itself.
"""
import runpy
import sys

try:
    guard = sys.argv[1]
    sys.argv = sys.argv[1:]
    runpy.run_path(guard, run_name='__main__')
except SystemExit as stop:
    sys.exit(stop.code if stop.code in (0, 2) else 2)
except BaseException as error:
    sys.stderr.write(f'BLOCKED: the guard {sys.argv[0]!r} failed ({type(error).__name__}: {error}).\n')
    sys.exit(2)
sys.stderr.write(f'BLOCKED: the guard {sys.argv[0]!r} ended without a decision.\n')
sys.exit(2)
