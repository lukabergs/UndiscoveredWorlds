"""Area-weighted comparisons and a fixed, inspectable global ranking rule."""
import numpy as np
from climate_report_data import magnitude,weights,FIELDS
from climate_palettes import DEFINITION

VERSION='global-balanced-v1'
def quantile(v,w,q):
    order=np.argsort(v);v=v[order];w=w[order]
    return float(np.interp(q*w.sum(),np.cumsum(w)-.5*w,v))
def bins(field):
    lo,hi=DEFINITION['scales'][field]
    middle=np.linspace(lo,hi,13)
    return np.concatenate(([-np.inf],middle,[np.inf]))
def evaluate(model,ref,mask,field):
    vector=model.ndim==3
    a,b=magnitude(model),magnitude(ref);w=weights(a.shape)[mask];av,bv=a[mask],b[mask]
    mean=lambda x:float(np.average(x,weights=w))
    err=magnitude(model-ref) if vector else np.abs(model-ref)
    e=err[mask];edges=bins(field)
    hist=lambda x:(np.histogram(x,bins=edges,weights=w)[0]/w.sum()*100).tolist()
    ah,bh=hist(av),hist(bv)
    aq={str(q):quantile(av,w,q) for q in (.1,.25,.5,.75,.9,.95,.99)}
    bq={str(q):quantile(bv,w,q) for q in (.1,.25,.5,.75,.9,.95,.99)}
    same=mean(np.digitize(av,edges)==np.digitize(bv,edges))*100
    rmse=np.sqrt(mean(e*e));p90=quantile(e,w,.9)
    tv=.5*np.abs(np.asarray(ah)-bh).sum()/100
    am,bm=mean(av),mean(bv);variance=mean((av-am)**2)*mean((bv-bm)**2)
    correlation=mean((av-am)*(bv-bm))/np.sqrt(variance) if variance>1e-20 else None
    rows=[['RMSE',float(rmse),0],['Error P90',p90,0],['Same-bin area (%)',same,100],
          ['Distribution distance (%)',tv*100,0],['Bias',am-bm,0]]
    if correlation is not None:rows.append(['Spatial correlation',float(correlation),1])
    if vector:
        labels={'RMSE':'Vector RMSE','Error P90':'Vector error P90','Same-bin area (%)':'Same speed-bin area (%)',
                'Distribution distance (%)':'Speed distribution distance (%)','Bias':'Speed bias','Spatial correlation':'Speed spatial correlation'}
        rows=[[labels.get(name,name),value,target] for name,value,target in rows]
        # The directional domain depends only on reference motion, never candidate speed.
        active=mask&(b>max(.05,.1*bm));aw=weights(a.shape)[active]
        dot=np.sum(model*ref,axis=0);den=np.maximum(a*b,1e-30)
        angle=np.degrees(np.arccos(np.clip(dot/den,-1,1)))
        angle=np.where(a>.01,angle,90.)
        direction=float(np.average(angle[active],weights=aw)) if aw.size else 0.
        within=float(np.average(angle[active]<=30,weights=aw)*100) if aw.size else 100.
        speed_rmse=np.sqrt(mean((av-bv)**2))
        rows.extend([['Speed RMSE',float(speed_rmse),0],['Direction error (°)',direction,0],['Direction within 30° (%)',within,100]])
        for i,name in enumerate(('East component RMSE','North component RMSE')):
            rows.append([name,float(np.sqrt(mean((model[i][mask]-ref[i][mask])**2))),0])
        scale=max(float(np.sqrt(mean(bv*bv))),1e-6)
        loss=.45*rmse/scale+.20*p90/max(bq['0.9'],scale*.1)+.20*direction/180+.15*tv
    else:
        lo,hi=DEFINITION['scales'][field]
        scale=max(quantile(bv,w,.95)-quantile(bv,w,.05),(hi-lo)*.01)
        loss=.50*rmse/scale+.20*p90/scale+.20*(1-same/100)+.10*tv
    return dict(rows=rows,model_hist=ah,reference_hist=bh,model_percentiles=aq,reference_percentiles=bq,
                loss=float(loss),cells=int(mask.sum()),area_weight=float(w.sum()))

def rank(metrics,baseline,ids):
    weights_by_field={f:info[3] for f,info in FIELDS.items() if info[3]}
    assert abs(sum(weights_by_field.values())-1)<1e-12
    scores={};losses={}
    for run in ids:
        losses[run]={}
        for field,weight in weights_by_field.items():
            seasons=['0'] if field=='annual_rain' else list(map(str,range(4)))
            entries=[metrics[str(run)][field][s]['loss'] for s in seasons]
            if not all(np.isfinite(entries)):raise ValueError(f'Invalid ranking input: {run}/{field}')
            losses[run][field]=float(np.mean(entries))
        scores[run]=100*sum(weights_by_field[f]*loss for f,loss in losses[run].items())
    candidates=[run for run in ids if run!=baseline]
    best=min(candidates,key=lambda run:(scores[run],run)) if candidates else baseline
    selected=best if scores[best]<scores[baseline]-1e-9 else baseline
    return dict(version=VERSION,baseline=baseline,best_in_batch=best,recommended_baseline=selected,
                scores=scores,field_losses=losses,field_weights=weights_by_field,
                improvement_percent=100*(1-scores[best]/scores[baseline]),
                lower_is_better=True,visual_override_allowed=True,
                formula='100 × weighted mean of dimensionless field losses; four seasons equally weighted. Scalar: .50 RMSE/reference P95–P05 spread + .20 error P90/same spread + .20 spatial-bin mismatch + .10 histogram total variation. Vector: .45 vector RMSE/reference RMS speed + .20 error P90/reference speed P90 + .20 direction error/180 + .15 speed-histogram total variation. Fixed common masks; calm candidate direction penalized 90 degrees where reference is active. Auxiliary quantities and mismatched physical definitions have zero weight. No Köppen component. Ties retain the baseline.')
