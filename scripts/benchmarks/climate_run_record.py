"""Write a compact per-run record; visual work belongs to the shared batch gallery."""
import html
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
MARKER='standard-report-pending'

def write_run_record(run_id,elapsed_seconds=None):
    directory=ROOT/'runs/reports';directory.mkdir(exist_ok=True)
    registry=json.loads((ROOT/'runs/registry/climate.json').read_text(encoding='utf-8'))
    record=next(r for r in registry['runs'] if r.get('id',r.get('run_id'))==run_id)
    selection_path=directory/'climate-global-current-selection.json'
    selection=json.loads(selection_path.read_text()) if selection_path.exists() else {}
    baseline=selection.get('selected_run')
    log=f'../logs/climate/{run_id}.log'
    packet=dict(run_id=run_id,baseline=baseline,elapsed_seconds=elapsed_seconds,
                registry='runs/registry/climate.json',log=f'runs/logs/climate/{run_id}.log',
                status='Simulation finished; batch comparison and ranking follow numerical analysis.')
    (directory/f'{run_id}-record.json').write_text(json.dumps(packet,indent=2)+'\n')
    # No historical reference IDs, handpicked charts or copied console dump.
    body=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>{run_id}</title>
<!-- {MARKER} -->
<style>body{{font:16px system-ui;background:#101923;color:#e9f0f5;margin:2rem}}a{{color:#92cfff}}</style>
<h1>{run_id}</h1><p>{html.escape(packet['status'])}</p>
<p><a href="{run_id}.md">Run record</a> · <a href="{log}">Raw execution log</a> · <a href="climate-latest-comparison.html">Latest visual comparison</a></p></html>'''
    (directory/f'{run_id}.html').write_text(body,encoding='utf-8')
    (directory/f'{run_id}.md').write_text(f'# {run_id}\n\nBaseline: {baseline}. Output width: {record.get("horizontal_resolution","unavailable")}. Elapsed seconds: {elapsed_seconds}.\n\nThe batch development report records the hypothesis, analysis and decision. The raw log contains execution diagnostics.\n\n[Raw log]({log}) · [Machine-readable receipt]({run_id}-record.json) · [Visual comparison](climate-latest-comparison.html).\n',encoding='utf-8')
    print(f'Run {run_id}: execution record written; standard batch gallery follows analysis.')
