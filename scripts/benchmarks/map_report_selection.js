// Report selection is independent of simulation origin and rendering.
const MapReportSelection = {
  endpoints(states) {
    return states.map((state,index)=>({state,index})).filter(({state})=>['initial','final'].includes(state.label));
  },
  runs(runIds,a,b,all) { return all?[...runIds]:[...new Set([a.run,b.run])]; },
  snapshots(states,runIds,a,b,all) {
    const selected=a===b?[a]:[a,b];
    if(all) for(const run of runIds) {
      if(selected.some(s=>s.run===run))continue;
      const state=states.find(s=>s.run===run&&s.label===b.label);
      if(state)selected.push(state);
    }
    return selected;
  }
};
if(typeof module!=='undefined')module.exports={MapReportSelection};
