/* Headless interaction and link checks; deliberately does not inspect map appearance. */
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict');
const file=path.resolve(process.argv[2]),folder=path.dirname(file),source=fs.readFileSync(file,'utf8');
const elements=new Map(),images=[],draws=[];
class Element {
  constructor(tag='div'){this.tag=tag;this.children=[];this.attributes={};this.listeners={};this.style={};this.textContent='';this._value=null;}
  setAttribute(k,v){assert.ok(!/NaN|Infinity|undefined/.test(String(v)),`${k}: ${v}`);this.attributes[k]=String(v);if(k==='value')this._value=String(v);}
  appendChild(e){this.children.push(e);return e;}
  replaceChildren(){this.children=[];}
  replaceWith(e){throw Error('Image failed to load');}
  addEventListener(k,fn){(this.listeners[k]??=[]).push(fn);}
  dispatch(k){for(const fn of this.listeners[k]||[])fn();}
  get value(){return this._value??this.children.find(c=>c.tag==='option')?.value??this.textContent;}
  set value(v){this._value=String(v);}
  getContext(){return {drawImage(...args){assert.ok(args.slice(1).every(Number.isFinite));draws.push(args.slice(1));}};}
}
for(const m of source.matchAll(/<select id="([^"]+)">([\s\S]*?)<\/select>/g)){
  const el=new Element('select');for(const o of m[2].matchAll(/<option(?: value="([^"]*)")?>([^<]*)<\/option>/g)){
    const e=new Element('option');e.textContent=o[2];e.value=o[1]??o[2];el.appendChild(e);
  }elements.set(m[1],el);
}
const d=new Element('script');d.textContent=source.match(/<script id="comparison-data" type="application\/json">([\s\S]*?)<\/script>/)[1];elements.set('comparison-data',d);
for(const m of source.matchAll(/<input id="([^"]+)" type="checkbox"([^>]*)>/g)){const e=new Element('input');e.checked=m[2].includes('checked');elements.set(m[1],e);}
const document={getElementById(id){if(!elements.has(id))elements.set(id,new Element());return elements.get(id);},createElement:t=>new Element(t),createElementNS:(_,t)=>new Element(t)};
class FakeImage {constructor(){this.naturalWidth=2048;this.naturalHeight=1024;}set src(v){images.push(v);this.onload();}}
const requestedRun=JSON.parse(d.textContent).runs.at(-1);
const window={location:{search:`?run=${requestedRun}`}},context=vm.createContext({document,window,Image:FakeImage,URLSearchParams});
vm.runInContext(fs.readFileSync(path.resolve('scripts/benchmarks/climate_report_dashboard.js'),'utf8'),context);
const app=window.climateComparison,data=app.data;
assert.equal(Number(elements.get('run').value),requestedRun);
assert.equal(source.match(/<title>(.*?)<\/title>/)[1],`(${data.baseline}) - ${Math.min(...data.runs)}-${Math.max(...data.runs)}`);
assert.ok(elements.get('run').children.every(x=>/^\d+$/.test(x.textContent)));
let combinations=0;
for(const run of data.runs)for(const field of Object.keys(data.fields))for(const season of ['0','1','2','3'])for(const mode of ['errors','distribution','percentiles']){
  elements.get('run').value=run;elements.get('field').value=field;elements.get('season').value=season;elements.get('table-mode').value=mode;app.update();combinations++;
  const c=app.current(),br=app.metricRows(data.metrics[String(data.baseline)][field]?.[c.season],mode,field),sr=app.metricRows(data.metrics[String(run)][field]?.[c.season],mode,field);
  const rows=elements.get('metrics').children.slice(1);
  if(!br.length){assert.equal(rows.length,1);continue;}
  assert.equal(rows.length,br.length);
  rows.forEach((row,i)=>{
    const target=br[i][2],errors=[Math.abs(br[i][1]-target),Math.abs(sr[i][1]-target)];
    const best=Math.min(...data.runs.map(r=>{const x=app.metricRows(data.metrics[String(r)][field]?.[c.season],mode,field);return x[i]?Math.abs(x[i][1]-x[i][2]):Infinity;}));
    for(let k=0;k<2;k++){
      const child=row.children[k+1].children[0];
      assert.equal(child.tag==='strong',Math.abs(errors[k]-Math.min(...errors))<1e-8);
      assert.equal(child.attributes.class==='best',Math.abs(errors[k]-best)<1e-8);
    }
  });
}
const before=JSON.stringify(elements.get('metrics').children);elements.get('region').value='North Pacific';app.update();assert.equal(JSON.stringify(elements.get('metrics').children),before);
elements.get('outlines').checked=false;
let flowCombinations=0;
for(const run of data.runs)for(const field of ['current','surface_wind','850_wind','500_wind'])for(const season of ['0','1','2','3'])for(const style of ['magnitude','lic','particles']){
  elements.get('run').value=run;elements.get('field').value=field;elements.get('season').value=season;elements.get('style').value=style;images.length=0;app.update();
  assert.equal(elements.get('style').disabled,false);assert.equal(images.length,3);
  const key=style==='magnitude'?field:field+'_'+style;
  assert.equal(images[1],data.maps[`${run}/${key}/${season}`]);assert.equal(images[2],data.references[`${key}/${season}`]);flowCombinations++;
}
elements.get('field').value='rain';elements.get('style').value='magnitude';elements.get('smooth').checked=false;images.length=0;app.update();
assert.equal(images[1],data.nearest_maps[`${elements.get('run').value}/rain/${elements.get('season').value}`]);
elements.get('outlines').checked=true;images.length=0;app.update();assert.equal(images.filter(p=>p===data.continent_outlines).length,3);
draws.length=0;const canvas=new Element('canvas');app.paint(canvas,new FakeImage(),data.regions['North Pacific']);assert.equal(draws.length,2);assert.equal(draws[0][6]+draws[1][6],2048);
for(const field of Object.values(data.metrics))for(const seasons of Object.values(field))for(const stat of Object.values(seasons)){
  assert.ok(Math.abs(stat.model_hist.reduce((a,b)=>a+b,0)-100)<1e-8);
  assert.ok(Math.abs(stat.reference_hist.reduce((a,b)=>a+b,0)-100)<1e-8);
}
const links=new Set([...Object.values(data.maps),...Object.values(data.references),...Object.values(data.nearest_maps),...Object.values(data.nearest_references),data.continent_outlines]);
for(const link of links){const p=path.resolve(folder,link),b=fs.readFileSync(p);assert.equal(b.subarray(1,4).toString(),'PNG');assert.equal(b.readUInt32BE(16),2048);assert.equal(b.readUInt32BE(20),1024);}
for(const link of [data.report,...[...source.matchAll(/(?:href|src)="([^"]+)"/g)].map(m=>m[1])])if(!/^https?:/.test(link))assert.ok(fs.statSync(path.resolve(folder,link)).isFile(),link);
console.log(JSON.stringify({passed:true,combinations,flowCombinations,map_files:links.size,smoothing_toggle:true,outline_toggle:true,dateline_wrap:true,numeric_run_labels:true,highlighting:'Every metric checked independently',visual_review:'User only'}));
