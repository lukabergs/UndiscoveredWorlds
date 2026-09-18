Is it necessary to keep ripgrep inside the project instead of having it as an external dependence?

Review (2026-09-03):
No. C:/Users/Luka/bin/rg.exe supplies ripgrep 15.1.0 independently of this repo, and no active build/script depended on the bundled copy. Removed tools/ripgrep/ and its ZIP; rg --version still succeeds. Other developers can install ripgrep on PATH.
