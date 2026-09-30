"""Save a small Markdown/JSON benchmark ledger without building HTML or animations."""
import argparse
import json
from pathlib import Path
from tectonics_experiment import write_json

ROOT=Path(__file__).resolve().parents[2]


def summarize(baseline,first,last):
    rows=[]
    for run in [baseline,*range(first,last+1)]:
        path=ROOT/f'runs/diagnostics/tectonics/{run}/batch-summary.json'
        if not path.exists():continue
        result=json.loads(path.read_text())
        plan=ROOT/f'runs/reports/tectonics/{run}-plan.json'
        if plan.exists():result['plan']=json.loads(plan.read_text())
        rows.append(result)
    target=ROOT/f'runs/reports/tectonics/{baseline}-{last}-metrics'
    write_json(target.with_suffix('.json'),rows)
    text=['# Tectonics batch measurements','',
        'All fractions use native equal-area cells; Earth uses spherical area. Runtime includes native export; concurrent compilation can affect small differences. Lower Earth curve gap is a distribution diagnostic, not a complete physical/visual score. Comparison baseline: '+str(baseline)+'. This report does not alter the registry baseline.','',
        '|Run|Seed|Seconds|Land %|Peak m|Mean depth m|Below −200 m %|Below −1000 m %|Above 2000 m %|Earth gap m|Ocean detail m|Coast p95 m|Splits|Fusions|',
        '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        f=r['fraction_below']['final'];e=r['events']
        values=[r['run'],r['seed'],f"{r['elapsed_seconds']:.1f}",f"{r['land_percent']:.2f}",f"{r['peak_m']:.0f}",
            f"{r['mean_depth_m']:.0f}",f"{100*f['-200']:.2f}",f"{100*f['-1000']:.2f}",f"{100*(1-f['2000']):.2f}",
            f"{r['earth_gap_m']:.0f}",f"{r['ocean_detail_m']:.0f}",f"{r['coast_p95_m']:.0f}",e['splits'],e['fusions']]
        text.append('|'+ '|'.join(map(str,values))+'|')
    text+=['','## Hypotheses and decisions','']
    for r in rows:
        p=r.get('plan',{})
        text.append(f"- **{r['run']}**, parent {p.get('parent',baseline)}: {p.get('hypothesis','Baseline / seed sensitivity')}")
        if p.get('decision'):text.append('  '+p['decision'])
    target.with_suffix('.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    print(f'Saved {len(rows)} run summaries to {target}.md')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--baseline',type=int,default=42);p.add_argument('--first',type=int,default=45);p.add_argument('--last',type=int,default=62)
    a=p.parse_args();summarize(a.baseline,a.first,a.last)
