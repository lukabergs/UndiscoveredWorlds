# /// script
# dependencies = ["numpy", "pillow"]
# ///
"""Compact global and regional numerical review; never a visual judgement."""
import argparse
import json
from pathlib import Path
import numpy as np
from climate_report_data import Archive, magnitude, weights
from climate_report import mask_for

def main():
    p=argparse.ArgumentParser();p.add_argument('--archive',type=Path,required=True);p.add_argument('--runs',type=int,nargs='+',required=True)
    args=p.parse_args();archive=Archive(args.archive,args.runs)
    metrics=json.loads((args.archive/'metrics.json').read_text())
    ranking=json.loads((args.archive/'ranking.json').read_text())
    for run in args.runs:
        print('\nRUN',run,'COST',round(ranking['scores'][str(run)],4))
        for field,loss in ranking['field_losses'][str(run)].items():
            baseline=ranking['field_losses'][str(ranking['baseline'])][field]
            entries=metrics[str(run)][field]
            row=lambda name: np.mean([next(value for label,value,_ in e['rows'] if label==name) for e in entries.values()])
            vector='Vector RMSE' in [x[0] for x in next(iter(entries.values()))['rows']]
            print(field, 'loss%',round(100*(loss/baseline-1),2),'rmse',round(row('Vector RMSE' if vector else 'RMSE'),3),
                'bias',round(row('Speed bias' if vector else 'Bias'),3),
                'P99/ref',round(np.mean([x['model_percentiles']['0.99'] for x in entries.values()]),3),
                round(np.mean([x['reference_percentiles']['0.99'] for x in entries.values()]),3),
                'direction',round(row('Direction error (°)'),1) if vector else '')
        fields,land,masks=archive.load(run)
        annual=fields['annual_rain_0'];refannual=archive.reference('annual_rain',0,annual.shape[-1])
        valid=np.isfinite(refannual);aw=weights(annual.shape)[valid]
        for label,lo,hi in [('moderate annual',1000,2500),('heavy annual',4500,np.inf)]:
            area=lambda x:float(np.average((x[valid]>=lo)&(x[valid]<hi),weights=aw)*100)
            print(label,'area% model/ref',round(area(annual),3),round(area(refannual),3))
        regions={'Pacific equator':(-2.5,2.5,-160,-90),'Pacific north':(3,10,-160,-90),'Pacific south':(-10,-3,-160,-90),
                 'Arabian Sea':(5,20,50,75),'India':(5,30,70,90),'Southern Ocean':(-65,-40,-180,180),
                 'Kuroshio':(20,40,125,150),'Agulhas':(-45,-20,15,40)}
        for region,(south,north,west,east) in regions.items():
            values=[]
            for field in (('rain','water') if 'Pacific' in region or region=='India' else ('surface_wind','current','sst_oisst')):
                for s in (0,2):
                    a=fields[f'{field}_{s}'];b=archive.reference(field,s,a.shape[-1]);shape=magnitude(a).shape
                    lat=90-(np.arange(shape[0])+.5)*180/shape[0];lon=-180+(np.arange(shape[1])+.5)*360/shape[1]
                    mask=mask_for(field,shape,land,masks,s)&np.isfinite(magnitude(b))
                    mask &= (lat[:,None]>=south)&(lat[:,None]<=north)&(lon[None,:]>=west)&(lon[None,:]<=east)
                    if not mask.any():continue
                    mean=lambda v:float(np.average(magnitude(v)[mask],weights=weights(shape)[mask]))
                    values.append(f'{field}/{s}: {mean(a):.2f}/{mean(b):.2f}')
            print(region,', '.join(values))

if __name__=='__main__':main()
