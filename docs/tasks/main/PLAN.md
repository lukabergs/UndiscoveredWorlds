# Goal

Organize the repository and local artifacts for a modular world engine.

# Status

Reorganized directly on `main`: 7,062 files moved, including ignored artifacts.
1,336 maps grouped by category; seven workbooks preserved; seven unidentified
files moved to `meta/`. README and ignore rules updated. Original tracked-file
blobs retained in staged moves; documentation additions remain uncommitted.

# Verification

- `python out/reorganization/reorganize.py verify`: 17,856 file destinations,
  155 SHA-256 hashes (including every workbook), and six junctions passed.
- `python out/reorganization/audit.py`: map names, 119 relocated tracked files,
  README links, root layout, and `main` branch passed.
- `git diff --check` and `git diff --cached --check`: passed.
- `git diff --numstat -- src`: empty; moved source contents unchanged.
- Follow-up: all 85 source files are exact `R100` renames; the legacy folder
  and its empty parents were removed, and `Test-Path .../extra` returns false.

# Remaining

- CMake/include/runtime/script paths and automatic output routing intentionally
  deferred. Build and runtime test tiers not run.
- Detailed original-to-new paths: `out/reorganization/moves.json` (local).
