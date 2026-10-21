# Agent workflow: roles, models, checks

Decided with the user on 2026-10-06. Model strength is allocated by the cost of being wrong: planning and review get the strongest model, implementation gets a competent one with a good brief, and verification is done by scripts that an agent only interprets. Everything that can be deterministic (builds, bundle comparisons, tests, parameter and layout checks, git, the ledger) is outside the models.

## Roles and models

Claude models only (Sonnet and Haiku are not used for now). The adversarial reviewer is GPT-6 Astra through Codex (`~/.codex/config.toml`), run by the user outside these sessions and without context from them.

| Role | Who | Model, effort | Does | Does not |
|---|---|---|---|---|
| Planner / lead | the interactive session | Opus 5.5, `high` (`xhigh` for a design question the plan cannot settle; never routinely) | understands the request, writes or revises the plan (goal, decisions, step table, acceptance checks, risks), answers the user, reviews phase results | implement steps itself when a worker can |
| Orchestrator | the session that runs a phase | Opus 5.5, `high` | holds the step list, briefs workers, runs disjoint steps in parallel, integrates, runs the step check, writes ledger entries and commits, decides retry or stop | take design decisions the plan does not settle: it records the question in the plan and continues with independent work; if the whole phase is blocked, it stops with a handoff note |
| Worker (`.claude/agents/worker.md`) | subagent | Opus 5.5, `high` | one bounded step: code, configs, tests, the step's named checks | surveys, architecture, unasked checks, large runs |
| Mechanic (`.claude/agents/mechanic.md`) | subagent | Opus 5.5, `medium` | edits by recipe: renames, rewrites by pattern, moves, generated files, list and reference updates | judgement calls (left and reported) |
| Surveyor (`.claude/agents/surveyor.md`) | subagent, read-only | Opus 5.5, `medium` | fan-out reading of code, data and sources; evidence with locations | recommendations, decisions |
| Reviewer (`.claude/agents/reviewer.md`) | subagent, read-only, fresh context | Opus 5.5, `high` | checks a finished step against its plan step and the evidence behind its ledger entry | editing, running |
| Verifier (`.claude/agents/verifier.md`) | subagent, read-only | Opus 5.5, `high` | visual checklist of an evaluated run | paused since 2026-10-05; the user says when it resumes |
| Adversarial reviewer | Codex, by the user | GPT-6 Astra, `xhigh` | assumes the step is wrong and looks for the most consequential failure, on the commit range and the plan alone | receives no reasoning from the Claude sessions |

Fable 5.1 is set aside (the user, 2026-10-17: its token usage is too high): the planner, the reviewer and any `design` agent run on Opus 5.5 at the efforts above; the `learned` notes below that name Fable are records. Effort escalation: `xhigh` or `max` only after a failed attempt at `high` on the same step, and the failure is recorded in the ledger entry (`lessons`). An agent's model and effort are fixed in its definition; the Agent tool's `model` parameter overrides the model only. When a step is harder than its definition's effort, the orchestrator does it in-session rather than inventing a new agent.

A human-facing decision (what a world should look like, a data deletion under `refs/` or `runs/`, a public contract) is the user's; the planner writes the recommendation into the plan and the orchestrator continues with what does not depend on it (`self-sufficiency`, AGENTS.md).

## A phase, step by step

Since 2026-10-20 a session runs a **wave**, not one phase: up to four C++ worktree workers from different lanes (file ownership), each taking the lowest-numbered open step of its lane in `docs/architecture/master-plan.md`, plus any number of Python-only steps and surveyors; timing steps run alone at the wave's end (master plan §2). The steps below apply to each step of the wave.

1. The planner's session ends with the plan section updated and the prompt for the orchestrating session printed in the chat (never saved in the plan; consolidation plan §8). The prompt points to documents; it repeats nothing an agent can read.
2. The orchestrator reads the plan section, the ledger entries it names and the files of the first steps, then builds the step order: steps with disjoint file sets run in parallel in background workers; steps that share files run in sequence.
3. For each step: brief → work → the step's deterministic checks → the orchestrator reads the diff → `reviewer` when the step changes behaviour, a contract or a kernel (not for a one-line fix or a documentation edit) → fixes by the same worker (its context is intact: continue it with SendMessage) → ledger entry → one commit. Duration is measured commit to commit (`AGENTS.md`).
4. The phase closes with the one judging run the plan names (if any), the plan's status note, §Open items, the memory notes, and the prompt for the next phase.

### The brief

The brief (orchestrator to worker, mechanic or surveyor) is the agent's only context. Its fields, in this order (R required, O optional; write `none` under an empty one):

- **Goal** (R): one sentence.
- **Step** (R): the step id and the plan section (`master plan 2.8`, §3 Phase 2).
- **Lane, start commit** (R): the lane and the commit to work from.
- **Files to change** / **Files to leave alone** (R): the second includes the files of workers running in parallel.
- **Ruled out** (O): what was already tried or decided against, with the ledger ids.
- **Checks** (R): the exact commands and their expected results.
- **Do not run** (R): the checks, suites, builds and runs that stay out.
- **Commit** (R): yes or no.
- **Report** (R): "the headings of 'The report'" (a surveyor's: its answers and `## Unverified`), plus any extra the step needs.

A brief that cannot be written clearly means the step is not ready to delegate: the planner refines the plan first.

Parallel workers get disjoint file sets. Two steps that must touch the same file are one step, or the second waits. When a worker must experiment destructively on shared files, run it with `isolation: worktree`; the orchestrator merges.

Learned in erosion Phase B (2026-10-06), where four C++ workers ran at once with `isolation: worktree`: each worktree builds its own `out/build/headless` (a fresh configure takes a few minutes), runs its own minimal bundles and compares them with the main tree's bundles by absolute path, so the steps do not mix in each other's builds and runs; the worktree's branch starts at the repository's first commit (`b309ae0`), not at the head, so the brief names the commit to start from and the worker fast-forwards or branches from it; `uv` creates a `.venv` in each worktree on its first call. The worker commits on its branch; the orchestrator brings each commit into main with `git cherry-pick -n`, in an order that makes each step check attributable (resolving by hand two steps' edits of one function), adds its ledger entry and commits; it removes the worktrees and branches when the phase's steps are in. Workers never run `step_check.py` (it builds the shared `out/build/headless`) or the tracking build; a step that touches only scripts and configs runs in the main tree on a copy of the executable.

Learned in erosion Phase C (2026-10-07), with up to four worktree workers at once and reviews that sent most steps back once: a worktree has no `refs/` (local-only); `git cherry-pick -n A B` stops at the first conflict and cannot continue over staged changes, so a worker's later commits (review fixes) come in as `git diff <old> <new>` applied with `git apply --3way --index`; a worker's numbers come from its own branch, so the orchestrator reruns the step check on the merged tree before the entry cites them, and a reviewer told so checks the right tree; access measured on the minimal world misses the writes of a field that stays zero there (prograded land, deltas), so a step that adds such a field gets one tracked run at 256 columns; a 256-column run of a step can restart from a checkpoint after `terrain.roughness` up to `hydrology.deltas_and_wetlands` (about a minute) but not through `layers.finalize`, whose tables a bundle cannot load. Fable 5.1's limit was reached at the first review; every Phase C review ran on Opus 5.5.

Learned in erosion Phase D (2026-10-07), by data loss: a worktree must never hold a junction or link into the main tree. `git worktree remove --force` on a worktree whose `refs` was a junction followed it and deleted most of the main tree's `refs/processed`; the `rmdir` meant to remove the link first had failed silently on a forward-slash path. Since then: a worktree has no `refs/`, and nothing links one in; a worker passes absolute paths into the main tree (`C:\dev\undiscovered-worlds\refs\...`) as command-line arguments, and a check that needs scripts reading `refs/` through the repository root (Earth recipes, validations, evaluations) runs from the main tree with the worker's executable given by absolute path. The orchestrator removes a worktree only after `Get-ChildItem -Recurse -Attributes ReparsePoint` finds nothing in it, and never with `--force` while a link may remain. A global hook (`~/.claude/hooks/guard_recursive_delete.py`) blocks recursive deletes that name `refs` or whose target holds a junction, and the creation of links to `refs`. `refs/` and the other untracked data are copied to `D:/dev/bak` (`scripts/uwtools/backup.py`).

Learned in erosion Phase G (2026-10-08): read a stage's code before planning around its name (`terrain.fjords` carved nothing). Every process reseeds the shared random stream from the world seed and its id (`seedprocessstream`), so removing a stochastic stage changes only what it wrote. A worker kept in the main tree while another step's change is staged works, if the orchestrator commits the staged step with explicit paths. A step that adds a field write needs two builds: the tracked run measures the undeclared write, `process_access.py` regenerates the include file, and only the rebuilt executable passes `--audit-access` (the evaluation's `written_by` comes from that build). A comparison with a run older than new fields must still report the run's side (the validation dispatcher skipped the whole kind). `uw.py <module> ledger --add` refuses an entry without `verdict`, `decision` and `do_not_retry`. Step-check labels are single-use: the next step check's reference is the latest label in the directory the previous phase used (Phase G: `g3-glacial-2` in `out/scratch/later-float`; the commit after it changed only an evaluation script and the access file's comment).

Learned in erosion Phase E's first session (2026-10-08): a question put to the user at the start of a phase may go unanswered for the whole session; plan the steps so that everything not gated by it (kernels, contracts, tests) runs first, and keep the gated step and its checks in one row that can start cold (here E2d took the real-bundle check from E3). For seamless (window-independent) results, the reviews found the faults the workers' tests missed: any decision taken over the whole window (an early stop of an iteration), any tie broken by a local index, any rule that lets information cross a pinned boundary (8-connected channels let a diagonal step cut between two pinned cells) and any overlap test whose windows share an origin (a world of 32 x 16 cells clips every halo at both poles, so rows never differ); an overlap test must assert bit identity, not only the plan's failure bound. When a header is split between parallel workers, the contracts step implements the shared part (here the patch metric) so that no worker waits. Moving floating-point code into a helper keeps results bit-identical only if every product keeps its order. A guard hook on recursive deletes matches the word `refs` anywhere in a command: keep a check of `refs/` out of the command that removes worktrees. At the phase's close (same evening, after the user's approval), a correction calibrated on synthetic terrain was seven times too weak on the real bundle (block-mean error 178 m against 26 m): measure a default on the judging world in scratch before registering runs; and evaluations are never overwritten, so a map fix after a registration renders into a second evaluation folder of the run (`--output runs/<module>/<id>/evaluation-2`).

Learned in erosion Phase R (2026-10-08 night to 2026-10-09): a contracts step that lets two workers run at once must also split the test files (the shared synthetic world went into its own header, `tests/simulations/landscape/regional_test_world.hpp`), and a feature that is off in the struct defaults is not tested by the existing overlap tests: each step adds overlap variants with its feature on, on a world wider than its widened box. Steps started on an unreviewed version of a parallel step (R3 and R4 on R2's first version) were rebased or re-applied and re-measured on the merged tree, and their reviews ran there; a reviewer told which tree its numbers come from caught numbers quoted from the wrong tree twice. Helper and sweep scripts are Python only (the user, 2026-10-09; a sweep calls the executable through `subprocess`), never `.ps1`. When the user pauses, stop every background agent with TaskStop, not only the wakeup loop; a completion notification may arrive for a run cut short, so check the worktree for an uncommitted diff and resume the same worker (its context is intact) with SendMessage. `git worktree remove --force` was refused by the permission check this phase: the finished worktrees wait for the user.

Learned in erosion Phase R2 (2026-10-09, the lead session itself on Fable 5.1): a contracts step that gives every parallel step its own translation unit and test file, with the composing file only calling inert hooks, lets three worktree workers run at once (R7, R8 and R9 took 24 to 34 min each, from one commit); a hook that must be a no-op keeps the original expression textually (a re-associated sum is identical only by luck) and is not read at its defaults (an interpolated constant field is not exactly the constant), so the parity check on every judging window is the proof; a contract that reads beyond a block states its reach against the land-share margin formula (two reviews derived it: e = k + hookreach + 2 fine cells, r = 1 + ceil(max(0, e - k/2) / k) blocks per share pass), and a contract that scales by a rate times a time states the time within which its physical argument holds; a worker's reproductions of earlier runs use the final build and parameter names, and the comparison's log is saved; `regional_overlap.py` refuses bundles made with different parameter lists, so a reproduction after a step adds parameters is compared field by field; a PowerShell here-string passed after `git commit -F -` is taken as a pathspec, so a commit message with several lines goes through a Bash heredoc or a file; reviews that run on a worker's branch are told so, and every number is re-measured on the merged tree before the entry cites it; a reviewer's attribution of an effect (the sides carried the block-mean shift) is checked by turning the term off before it is written down (the solve's lowering carried it).

Learned in erosion Phase E5 (2026-10-09, the lead session on Fable 5.1): when the implementing model's session limit stops parallel worktree workers after their code and tests are complete but before their measurements, their diffs (`git -C <worktree> diff | git apply`) land cleanly on the main tree when the contracts step gave each its own files; the lead finishes the measurements in-session, commits each step by its files (`git add <paths>`; the commits then share a timestamp, so the ledger records the real figure beside the commit-to-commit one) and has each commit reviewed. A contracts review that runs while the workers build on the contracts pays for itself: it found a placement rule that read the solved elevation in the innermost halo ring (not seam-safe; only the start surface, the sea and the channel beds agree there: a placement or eligibility rule reads a pure global quantity, the upsampled base), a probe that took a climate constant for the world's lapse rate (read the field the code reads), and a zero-filled field passed as present. An overlap variant must assert the feature on the shared cells, not somewhere in the window; a number inferred from a sweep's shape (the catchments) is measured before the entry cites it; a sweep whose variable saturates at a cap reports the cap as the control. (Since master plan 2.2 `step_check.py` has no `--directory`: it writes `out/scratch/step-check/<id>-<tag>` of the main tree and resolves an earlier label in `out/scratch/phase5` or `out/scratch/sphere/step-check` itself.) Relief forms of a few hundred metres are invisible on the 10 km elevation palette: the next form-cutting step needs a hillshade or on-minus-off composite in `scripts/evaluation/maps.py` first.

Learned in track A (2026-10-21, master plan Phases 1 and 2): a step that writes or changes a deleter (anything calling `remove_tree`, `unlink`, `rmdir`) cannot be checked before its commit, because `delete_guard.py` refuses any command naming an uncommitted deleting script, its tests included; the order is worker, reviewer reading the diff, orchestrator commit, then `py_compile`, the tests and the real run (a dry run first). Two review rounds were normal for each deleter (2.3, 2.7). The Agent harness locks its worktrees (`git worktree list` shows `locked claude agent ...`) and a finished worker can leave a shell holding the folder: stop that shell, `git worktree unlock` the path, then `worktrees.py retire` resumes. The repository's `.venv` has no pytest: name `python -m unittest FILE` in briefs. A worktree branched before a tool existed failed the shared pre-commit hook; the hook now skips a check its tree does not have. Reviewer findings went back to the same worker (SendMessage), and the triage rule kept every finding inside its step except two new rows (2.9, 4.3).

### The report

A worker's or mechanic's answer to the orchestrator, under these `## ` headings in this order (R required, O optional; `none` where empty; no narrative of the process). `scripts/check_records.py --report FILE --kind worker` checks the headings.

- `Files changed` (R): one line per file, what changed (mechanic: the occurrences changed per file).
- `Checks` (R): each command, its exit code and its result (identical or not, test counts); a mechanic also gives the search for forbidden results and its output.
- `Not done` (R): what the brief asked that was not done, with the reason (mechanic: the occurrences left and why).
- `Not run` (R): the checks skipped and what the minimal world cannot show; it goes into the ledger entry's `not_run`.
- `Questions` (R): each question for the plan with a recommended option.
- `Triage` (O): issues seen outside the step, each as one proposed row with the phase that owns it (`AGENTS.md`, "Tracks and triage").
- `Commit` (R when the brief says commit): the hash and the paths.
- `Duration` (R): start and end clock times and the minutes.

### Questions to the user

Every open question goes to the user in one message, and each is also recorded in the plan's §5; no prompt for the next session is printed while any is open (`AGENTS.md`). Per question:

- `Q<n>` (R): the question in one sentence.
- Options (R): lettered, each with its consequence in one line.
- Recommendation (R): the option and a one-line reason.
- Blocks (R): which steps wait on it and what continues meanwhile.
- Evidence (O): page paths, ledger ids.

### The wave handoff

Printed in the chat when a wave closes (master plan §2, step 8), under these headings:

- `Wave` (R): the track, the wave number and the start..end commit range (for the adversarial review).
- `Steps` (R): id, status, commit and duration of each step, equal to the plan's §4 rows.
- `Checks` (R): `check_records.py` exit 0; the `step_check.py` directory name (`<id>-<tag>`) to use as the next reference.
- `Looks` (R): the paths of the pages waiting for the user's look.
- `Learned` (R): where the session wrote what it learned.
- `Cleanup` (R): worktrees retired, branches dropped, trash entries.
- `Questions` (R): `none`, or "asked; no prompt".
- `Next prompt` (R only when `Questions` is `none`): the prompt for the next session.

### Status row and ledger entry

`scripts/check_records.py` checks both (in `step_check.py` and the pre-commit hook); rows and entries written before these fields are grandfathered in the script.

Status row (master plan §4, one per closed step):

- Step: a §3 step id, optionally followed by `(word ...)` (lower-case words and digits).
- Status: `done: ...`, `part: ...` or `stopped: ...`.
- Commit: the hashes (`h, h`), or `this commit`, which the next step replaces with the hash.
- Duration: `<who> N min; N min from <hash>` (the worker's time, then the time from the previous step's commit).

Ledger entry (strict fields, beside the required ones of `scripts/uwtools/ledger.py`):

- `phase`: `master plan <id>: ...` for a master-plan step (the id is a §3 step).
- `duration`: `start <hash> (<time>); ...`, with the time span in `N min` or `N h`.
- `not_run`: a list (the report's `Not run`).

### Verification is scripts first

- Default check of a step: `scripts/evaluation/step_check.py <tag> --reference <previous id or id-tag>` (headless build, the three seed-1 bundles at the minimal configuration with `--audit-access`, bundle comparison, `parameters.py --check`, `check_layout.py`, `check_structure.py`, `check_records.py`, `process_access.py --check`); `--desktop` only when `src/app` changed. Each check gets a new `out/scratch/step-check/<id>-<tag>` (`scripts/uwtools/naming.py`); without `--reference` the newest passing, compared check with the same grid and `--social` is the reference, and without one the check stops and names the earlier labels (`t8-cube4`, `t8-ll4`).
- Outputs go to D and the previous run of a topic is archived by the tools themselves (master plan 2.4, `naming.finish`): when a run passes, every earlier run of the same (kind, module, tag) not yet archived is archived (so a failed run in between leaves nothing behind) (`archive.py`: its re-renderable map PNGs and pages dropped, with a receipt), or, for a sandbox, step-check, geology or page folder, marked `superseded_by` in `runs/registry/scratch.json`; a registry baseline, a path named in `configs/sandbox/*.json` or `configs/evaluation/*.json` and an input of the new run (reference, `--baseline`, `--parent-run`, base bundle) are never touched. Sandbox pages (`terrain_sandbox.py`, `nested_window.py`, `compare_runs.py`, `hypsometry.py`, `terrain_metrics.py`) go to `out/evaluation/<id>-<tag>/` and link the scratch images by relative path; a tool's own outputs (climate diagnostics, images) go below the run's folder (`naming.tool_environment`). `storage.py status` lists what is still on C with the command that moves it.
- New behaviour comes with its own deterministic test: a GoogleTest on a synthetic grid of a few hundred cells for a kernel (`tests/`, run as `ctest -R "^<Suite>\."`), a profile check for a run-level invariant (`configs/evaluation/<module>.json`, `scripts/evaluation/checks.py`), a Python unit test for a script. The test states the expected numbers; "looks right" is not a test.
- One build per step, not per edit. No full ctest, Python suite or desktop launch unless the step changes what they cover, and then only the affected suite (`-R`, `-k`).
- Runs: the minimal configuration (128 x 64, `--duration 50`, `--climate-resolution 32 --climate-coupling-iterations 2`) for every "did anything change" question; one 256-column run with the default duration per phase close when the phase changes mountains or terrain (the plan names it); 512 columns only when the user asks. Before a run longer than a few minutes, ask whether a smaller one answers the same question and record what the small one cannot show.
- Agents interpret outputs (`compare_bundles.py` summaries, test counts, check reports, the maps the evaluation renders); they do not eyeball arrays or write per-run scripts. Observed values go into the ledger entry as numbers.
- A claim of improvement needs the evaluation (`uw.py <module> evaluate`) and, for visual modules, the user's verdict; the verifier is paused.

### Review

`reviewer` gets: the plan step, the commit range or the dirty files, the draft ledger entry, and nothing of the author's reasoning. Its verdict comes under these `## ` headings, in this order (`scripts/check_records.py --report FILE --kind reviewer` checks them):

- `Verdict` (R): `pass`, `pass with notes` or `rework`.
- `Findings` (R): numbered by consequence, each with file:line, what is wrong, and the smallest fix or check that settles it.
- `Claims` (R): one line per claim of the draft ledger entry: supported, unsupported or contradicted.
- `Not exercised` (R): what the minimal configuration and the recorded checks could not exercise.
- `Scope` (R): what was asked and not done; what was done and not asked.

Its findings go into the ledger entry (`observed`, `not_run`, `lessons`) or back to the worker. The adversarial review by GPT-6 Astra is the user's step: when a phase closes, the orchestrator lists the commit range and the plan section in the handoff so the user can hand them over; findings the user brings back are new steps or §Open items.

### Budget

- Launch independent subagents in one message; keep working while they run; never poll.
- Read with excerpts and searches, not whole large files; delegate wide reading to `surveyor`.
- Prefer the smallest configuration, suite and build that answers the question; say what was skipped.
- A session writes what it learned where the next agent reads it (`AGENTS.md`, the plan, the ledger, `docs/reference-guides/`, memory) before it ends.

## Configuration

- Agents: `.claude/agents/*.md` (frontmatter `model`, `effort`, `tools`, `disallowedTools`); documented fields at https://code.claude.com/docs/en/sub-agents.md. No `CLAUDE.md` in the repository: `AGENTS.md` is the project file (consolidation plan §8).
- Session model and effort: the user's `~/.claude/settings.json` (`model: opus`, `modelSettings` with `effortLevel: high`); the agents' models are fixed in their definitions regardless of the session's.
- Codex (adversarial review): `~/.codex/config.toml` (`gpt-6-astra`, `xhigh`); it reads `~/.codex/AGENTS.md` and the repository's `AGENTS.md`.
