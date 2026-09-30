/* Standard visual workspace; narrative and provenance live in the linked report. */
(() => {
  'use strict';
  const data=JSON.parse(document.getElementById('comparison-data').textContent), $=id=>document.getElementById(id);
  const colours=['#adb0ed','#ffc06f','#5ed9cc'];
  const add=(parent,tag,text,attrs={})=>{const e=document.createElement(tag);if(text!==null)e.textContent=text;for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);parent.appendChild(e);return e;};
  const svg=(parent,tag,attrs={},text)=>{const e=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;parent.appendChild(e);return e;};
  const fmt=x=>Number.isFinite(x)?(Math.abs(x)>=100?x.toFixed(1):x.toFixed(3)):'Unavailable';
  for(const run of data.runs)add($('run'),'option',String(run),{value:run});
  for(const [id,f] of Object.entries(data.fields))add($('field'),'option',f.label,{value:id});
  for(const region of Object.keys(data.regions))add($('region'),'option',region,{value:region});
  const requested=new URLSearchParams(window.location.search).get('run');
  $('run').value=requested && data.runs.includes(Number(requested)) ? requested : String(data.ranking.best_in_batch);$('field').value='rain';$('season').value='2';
  $('report-link').href=data.report;
  function current(){const field=$('field').value,season=field==='annual_rain'?'0':$('season').value,run=$('run').value;
    return {field,season,run,base:String(data.baseline),f:data.fields[field]};}
  function get(run,field,season){return data.metrics[run]?.[field]?.[season];}
  function metricRows(item,mode,field){
    if(!item)return [];
    if(mode==='errors')return item.rows;
    if(mode==='percentiles')return Object.keys(item.model_percentiles).map(q=>[`P${Number(q)*100}`,item.model_percentiles[q],item.reference_percentiles[q]]);
    const edges=data.fields[field].bins;
    return item.model_hist.map((v,i)=>[`${edges[i]===null?'−∞':fmt(edges[i])} to ${edges[i+1]===null?'+∞':fmt(edges[i+1])} (%)`,v,item.reference_hist[i]]);
  }
  function makeTable(c){
    const mode=$('table-mode').value,table=$('metrics');table.replaceChildren();
    const head=add(table,'tr',null);for(const title of ['Measurement',String(c.base),String(c.run),'Reference'])add(head,'th',title);
    const br=metricRows(get(c.base,c.field,c.season),mode,c.field),sr=metricRows(get(c.run,c.field,c.season),mode,c.field);
    if(!br.length){const row=add(table,'tr',null);add(row,'td','Exact reference unavailable',{colspan:4});return;}
    for(let i=0;i<br.length;i++){
      const [label,b,target]=br[i],s=sr[i][1],row=add(table,'tr',null);add(row,'td',label);
      const distances=[Math.abs(b-target),Math.abs(s-target)],pairBest=Math.min(...distances);
      const best=Math.min(...data.runs.map(run=>{const r=metricRows(get(String(run),c.field,c.season),mode,c.field);return r[i]?Math.abs(r[i][1]-r[i][2]):Infinity;}));
      [b,s].forEach((value,j)=>{const cell=add(row,'td',null),winner=Math.abs(distances[j]-pairBest)<1e-8,allBest=Math.abs(distances[j]-best)<1e-8;
        const node=winner?add(cell,'strong',fmt(value)):add(cell,'span',fmt(value));if(allBest)node.setAttribute('class','best');});
      add(row,'td',fmt(target));
    }
  }
  function chartLegend(parent,c){[c.base,c.run,'Reference'].forEach((name,i)=>svg(parent,'text',{x:80+i*235,y:22,fill:colours[i],'font-size':14},name));}
  function histogram(c){
    const el=$('histogram');el.replaceChildren();const b=get(c.base,c.field,c.season),s=get(c.run,c.field,c.season);if(!b)return;
    chartLegend(el,c);const series=[b.model_hist,s.model_hist,b.reference_hist],n=series[0].length,max=Math.max(...series.flat(),1),left=65,top=45,w=1000,h=210;
    for(let k=0;k<=4;k++){const y=top+h-k*h/4;svg(el,'line',{x1:left,x2:left+w,y1:y,y2:y,stroke:'#314050'});svg(el,'text',{x:6,y:y+4,fill:'#bac9d7','font-size':12},(max*k/4).toFixed(1)+'%');}
    series.forEach((values,j)=>values.forEach((v,i)=>{const rect=svg(el,'rect',{x:left+i*w/n+j*w/n/3+1,y:top+h-v/max*h,width:w/n/3-2,height:v/max*h,fill:colours[j]});svg(rect,'title',{},metricRows(b,'distribution',c.field)[i][0]+': '+v.toFixed(3)+'%');}));
    svg(el,'text',{x:65,y:286,fill:'#bac9d7','font-size':13},'Value bins, low → high; end bins include values outside the display range. '+c.f.unit);
  }
  function profile(c){
    const el=$('profile');el.replaceChildren();const region=$('profile-region').value,b=get(c.base,c.field,c.season)?.profiles[region],s=get(c.run,c.field,c.season)?.profiles[region];if(!b)return;
    chartLegend(el,c);const series=[b.model,s.model,b.reference],all=series.flat(),low=Math.min(...all),high=Math.max(...all),span=Math.max(1e-9,high-low);
    for(let k=0;k<=4;k++){const y=45+k*185/4;svg(el,'line',{x1:75,x2:1060,y1:y,y2:y,stroke:'#314050'});svg(el,'text',{x:3,y:y+4,fill:'#bac9d7','font-size':11},fmt(high-span*k/4));}
    series.forEach((values,i)=>{const points=values.map((v,j)=>`${75+(b.latitude[j]+90)/180*985},${230-(v-low)/span*185}`).join(' ');svg(el,'polyline',{points,fill:'none',stroke:colours[i],'stroke-width':2});});
    for(const lat of [-90,-60,-30,0,30,60,90])svg(el,'text',{x:75+(lat+90)/180*985,y:252,fill:'#bac9d7','text-anchor':'middle','font-size':12},String(lat)+'°');
    svg(el,'text',{x:75,y:276,fill:'#bac9d7','font-size':12},region+' · '+c.f.unit+' · latitude averages can hide longitudinal displacement');
  }
  function legend(c){
    const el=$('legend');el.replaceChildren();const [lo,hi]=c.f.range;
    const ramp=add(el,'div',null,{class:'ramp'});ramp.style.background='linear-gradient(to right,'+c.f.palette.map(([v,rgb])=>`rgb(${rgb.join(',')}) ${(v-lo)/(hi-lo)*100}%`).join(',')+')';
    const stops=add(el,'div',null,{class:'stops'});for(const [v] of c.f.palette)add(stops,'span',fmt(v));add(el,'div',c.f.unit);
  }
  function paint(canvas,image,box,overlay=false){
    const [south,north,west,east]=box,span=east>west?east-west:360-west+east;
    const x=(west+180)/360*image.naturalWidth,y=(90-north)/180*image.naturalHeight,sw=span/360*image.naturalWidth,sh=(north-south)/180*image.naturalHeight;
    if(!overlay){canvas.width=2048;canvas.height=Math.max(1,Math.round(2048*sh/sw));}
    const smooth=$('smooth')?.checked!==false,ctx=canvas.getContext('2d');ctx.imageSmoothingEnabled=overlay?false:smooth;
    canvas.style.imageRendering=smooth?'auto':'pixelated';
    const first=Math.min(sw,image.naturalWidth-x),fraction=first/sw;
    ctx.drawImage(image,x,y,first,sh,0,0,canvas.width*fraction,canvas.height);
    if(first<sw)ctx.drawImage(image,0,y,sw-first,sh,canvas.width*fraction,0,canvas.width*(1-fraction),canvas.height);
  }
  function panels(c){
    const el=$('panels');el.replaceChildren();const style=$('style').value,flow=['current','surface_wind','850_wind','500_wind'].includes(c.field),field=flow&&style!=='magnitude'?c.field+'_'+style:c.field;
    const nativeCells=$('smooth')?.checked===false&&field===c.field;
    const maps=nativeCells&&data.nearest_maps?data.nearest_maps:data.maps,refs=nativeCells&&data.nearest_references?data.nearest_references:data.references;
    const box=data.regions[$('region').value],paths=[maps[`${c.base}/${field}/${c.season}`],maps[`${c.run}/${field}/${c.season}`],refs[`${field}/${c.season}`]];
    paths.forEach((path,i)=>{const figure=add(el,'figure',null);add(figure,'figcaption',[c.base,c.run,'Reference'][i]);
      if(!path){add(figure,'div','Exact reference unavailable',{class:'missing'});return;}
      const canvas=add(figure,'canvas',null,{'aria-label':`${[c.base,c.run,'Reference'][i]} ${c.f.label}`});
      const image=new Image();image.onload=()=>{paint(canvas,image,box);if($('outlines')?.checked&&data.continent_outlines){const overlay=new Image();overlay.onload=()=>paint(canvas,overlay,box,true);overlay.src=data.continent_outlines;}};image.onerror=()=>{canvas.replaceWith(document.createTextNode('Image unavailable: '+path));};image.src=path;
      const links=add(figure,'div',null,{class:'downloads'});add(links,'a','Full export',{href:path,target:'_blank'});
    });
  }
  function update(){
    const c=current();$('season').disabled=c.field==='annual_rain';$('style').disabled=!['current','surface_wind','850_wind','500_wind'].includes(c.field);
    $('score-note').textContent=`Numerical best: ${data.ranking.best_in_batch} · Recommended baseline: ${data.ranking.recommended_baseline} · Selected cost: ${fmt(data.ranking.scores[c.run])} · Baseline cost: ${fmt(data.ranking.scores[c.base])} (lower is better)`;
    $('field-note').textContent=data.reference_audit[c.field]+(c.f.kind==='vector'?' Arrows show east/north direction.':'');
    makeTable(c);histogram(c);profile(c);legend(c);panels(c);
  }
  for(const name of ['run','field','season','table-mode','region','style','profile-region','smooth','outlines'])$(name)?.addEventListener('change',update);
  window.climateComparison={data,current,metricRows,paint,update};update();
})();
