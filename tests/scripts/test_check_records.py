"""Rules of scripts/check_records.py and the strict fields of scripts/uwtools/ledger.py (fixtures, no git but the last case)."""
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / 'tests/scripts/fixtures/records'
spec = importlib.util.spec_from_file_location('check_records', ROOT / 'scripts/check_records.py')
records = importlib.util.module_from_spec(spec)
spec.loader.exec_module(records)
ledger = records.ledger

COMMITS = {'aaaaaaa', 'bbbbbbb', 'ccccccc'}
known = COMMITS.__contains__
PLAN = (FIXTURES / 'plan.md').read_text(encoding='utf-8')
STEPS = records.plan_steps(PLAN)
LINES = (FIXTURES / 'world.jsonl').read_text(encoding='utf-8').splitlines()
GRANDFATHERED, STRICT = (json.loads(line) for line in LINES)


def ledger_findings(*entries):
    text = '\n'.join(json.dumps(entry) for entry in entries) + '\n'
    return records.check_ledgers({'world.jsonl': text}, known, STEPS)[0]


DROP = object()


def strict(**changes):
    """The fixture's strict entry with keys changed, or removed when set to DROP."""
    entry = dict(STRICT, **changes)
    return {key: value for key, value in entry.items() if value is not DROP}


def plan_with(row):
    return PLAN.replace('\n\n## 5.', f'\n{row}\n\n## 5.')


def plan_findings(text, modules=('world',)):
    return records.check_plan(text, known, {'world-0087', 'world-0088'}, set(modules))


def report(name):
    return (FIXTURES / name).read_text(encoding='utf-8')


class LedgerTests(unittest.TestCase):
    def assertFails(self, *entries, text=None):
        found = ledger_findings(*entries)
        self.assertTrue(found, 'expected a finding')
        if text:
            self.assertTrue(any(text in problem for *_, problem in found), found)

    def test_fixture_passes(self):
        self.assertEqual(ledger_findings(GRANDFATHERED, STRICT), [])

    def test_grandfathered_entry_without_duration_is_skipped(self):
        self.assertNotIn('duration', GRANDFATHERED)
        self.assertFalse(ledger.strict(GRANDFATHERED))
        self.assertEqual(ledger_findings(GRANDFATHERED), [])

    def test_strict_entry_without_duration_fails(self):
        self.assertFails(GRANDFATHERED, strict(duration=DROP), text='duration')

    def test_duration_without_start_fails(self):
        self.assertFails(GRANDFATHERED, strict(duration='worker 30 min from aaaaaaa'), text="start <commit hash>")

    def test_unknown_start_hash_fails(self):
        self.assertFails(GRANDFATHERED, strict(duration='start ddddddd (18:00); 30 min'), text='not a commit')

    def test_unknown_master_plan_step_fails(self):
        self.assertFails(GRANDFATHERED, strict(phase='master plan 99.9: nowhere'), text='99.9')

    def test_duplicate_id_fails(self):
        self.assertFails(GRANDFATHERED, STRICT, dict(STRICT), text='already used')

    def test_not_run_string_fails(self):
        self.assertFails(GRANDFATHERED, strict(not_run='nothing'), text='not_run')

    def test_unparsable_line_fails(self):
        found = records.check_ledgers({'world.jsonl': LINES[0] + '\n{broken\n'}, known, STEPS)[0]
        self.assertEqual([(path, line) for path, line, *_ in found], [('world.jsonl', 2)])

    def test_validate_refuses_new_strict_entry_without_duration(self):
        entry = strict(id='world-9999', parent=None, duration=DROP)
        with self.assertRaises(ValueError) as caught:
            ledger.validate(entry, [])
        self.assertIn('duration', str(caught.exception))
        ledger.validate(strict(id='world-9999', parent=None), [])

    def test_template_has_phase_and_duration(self):
        entry = ledger.template('world')
        self.assertEqual((entry['phase'], entry['duration']), ('', ''))
        self.assertTrue(set(ledger.STRICT) <= set(entry))


class PlanTests(unittest.TestCase):
    def assertRowFails(self, row, text):
        found = plan_findings(plan_with(row))
        self.assertTrue(any(text in problem for *_, problem in found), found)

    def test_fixture_rows_pass(self):
        # 1.1 is grandfathered with an empty duration; 4.1 is a good row; 4.2 names 'this commit'.
        self.assertEqual(plan_findings(PLAN), [])
        self.assertEqual(STEPS, {'1.1', '4.1', '4.2', '4.3'})

    def test_bad_step_format_fails(self):
        self.assertRowFails('| 4.3 (Later) | done: x | aaaaaaa | 5 min from aaaaaaa |', 'step is not')

    def test_step_missing_from_section_3_fails(self):
        self.assertRowFails('| 4.9 | done: x | aaaaaaa | 5 min from aaaaaaa |', 'not a row of section 3')

    def test_duplicate_step_fails(self):
        self.assertRowFails('| 4.1 | done: again | aaaaaaa | 5 min from aaaaaaa |', 'two status rows')

    def test_bad_status_fails(self):
        self.assertRowFails('| 4.3 | finished | aaaaaaa | 5 min from aaaaaaa |', 'status')

    def test_unknown_hash_fails(self):
        self.assertRowFails('| 4.3 | done: x | aaaaaaa, ddddddd | 5 min from aaaaaaa |', 'ddddddd is not a commit')

    def test_empty_duration_fails(self):
        self.assertRowFails('| 4.3 | done: x | aaaaaaa | |', "'<n> min'")

    def test_duration_without_from_fails(self):
        self.assertRowFails('| 4.3 | done: x | aaaaaaa | worker 5 min |', "'from <commit>'")

    def test_unknown_ledger_id_fails(self):
        self.assertRowFails('| 4.3 | done: x; world-0099 | aaaaaaa | 5 min from aaaaaaa |', 'world-0099')

    def test_this_commit_in_parentheses_passes(self):
        self.assertEqual(plan_findings(plan_with('| 4.3 | done: x | (this commit) | 5 min from aaaaaaa |')), [])

    def test_code_and_emphasis_marks_are_stripped(self):
        row = '| **4.3** | **done**: x | `aaaaaaa`, `bbbbbbb` | worker 5 min; 9 min from `ccccccc` |'
        self.assertEqual(plan_findings(plan_with(row)), [])
        self.assertRowFails('| 4.3 | done: x | `ddddddd` | 5 min from aaaaaaa |', 'ddddddd is not a commit')
        self.assertRowFails('| 4.3 | done: x | aaaaaaa | 5 min from `ddddddd` |', 'ddddddd, which is not a commit')

    def test_ledger_ids_follow_the_ledger_modules(self):
        row = plan_with('| 4.3 | done: x; climate-0001 | aaaaaaa | 5 min from aaaaaaa |')
        self.assertEqual(plan_findings(row, modules=('world',)), [])
        found = plan_findings(row, modules=('world', 'climate'))
        self.assertTrue(any('climate-0001' in problem for *_, problem in found), found)

    def test_step_parentheses_take_digits(self):
        self.assertEqual(plan_findings(plan_with('| 4.3 (part 2) | part: x | aaaaaaa | 5 min from aaaaaaa |')), [])

    def test_finding_names_line_and_step(self):
        found = plan_findings(plan_with('| 4.3 | finished | aaaaaaa | 5 min from aaaaaaa |'))
        line = PLAN.splitlines().index('| 4.2 | part: closes here | this commit | 1 h 5 min from bbbbbbb |') + 2
        self.assertEqual([(path, number, step) for path, number, step, _ in found], [(records.PLAN, line, '4.3')])


class ReportTests(unittest.TestCase):
    def test_good_worker_report_passes(self):
        self.assertEqual(records.check_report(report('worker-good.md'), 'worker'), [])

    def test_missing_not_run_fails(self):
        found = records.check_report(report('worker-missing.md'), 'worker')
        self.assertTrue(any("'## Not run' is missing" in problem for _, problem in found), found)

    def test_out_of_order_fails(self):
        text = report('worker-good.md').replace('## Not done', '## Placeholder').replace('## Not run', '## Not done')
        text = text.replace('## Placeholder', '## Not run')
        found = records.check_report(text, 'worker')
        self.assertTrue(any('out of order' in problem for _, problem in found), found)

    def test_empty_section_fails(self):
        text = report('worker-good.md').replace('none\n\n## Triage', '\n## Triage')
        found = records.check_report(text, 'worker')
        self.assertTrue(any("'Questions' is empty" in problem for _, problem in found), found)

    def test_fenced_code_is_not_a_heading(self):
        text = report('worker-good.md').replace('- step_check.py (not asked).',
                                                '- step_check.py (not asked).\n\n```\n## Checks\n## Files changed\n```')
        self.assertEqual(records.check_report(text, 'worker'), [])

    def test_verdict_ok_fails(self):
        found = records.check_report(report('verdict-bad.md'), 'reviewer')
        self.assertTrue(any('verdict' in problem for _, problem in found), found)

    def test_verdict_pass_with_notes_passes(self):
        text = report('verdict-bad.md').replace('\nok\n', '\n`pass with notes`\n')
        self.assertEqual(records.check_report(text, 'reviewer'), [])
        text = report('verdict-bad.md').replace('\nok\n', '\n**rework**: two findings\n')
        self.assertEqual(records.check_report(text, 'reviewer'), [])


class RepositoryTests(unittest.TestCase):
    def test_real_tree_is_clean(self):
        """Guards the live plan and ledgers on purpose: a malformed status row or entry fails this test."""
        self.assertEqual(records.main([]), 0)


if __name__ == '__main__':
    unittest.main()
