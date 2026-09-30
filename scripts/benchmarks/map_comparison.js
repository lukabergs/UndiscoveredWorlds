// Origin-independent comparison of immutable little-endian float64 field assets.
class MapFieldComparison {
  constructor(load = async url => {
    const response = await fetch(url);
    if (!response.ok) throw new Error('Field request failed: ' + response.status);
    return response.arrayBuffer();
  }) { this.load = load; this.cache = new Map(); }
  async field(asset) {
    let pending = this.cache.get(asset.url);
    if (!pending) {
      pending = this.load(asset.url).then(buffer => {
        if (buffer.byteLength !== asset.cells * 8) throw new Error('Field size mismatch');
        const view = new DataView(buffer), values = new Float64Array(asset.cells);
        for (let i=0;i<values.length;i++) {
          values[i] = view.getFloat64(i*8,true);
          if (!Number.isFinite(values[i])) throw new Error('Non-finite field value');
        }
        return values;
      });
      this.cache.set(asset.url,pending);
      if (this.cache.size>4) this.cache.delete(this.cache.keys().next().value);
      pending.catch(()=>this.cache.delete(asset.url));
    }
    return pending;
  }
  async compare(first,second,categorical=false) {
    if (first.cells!==second.cells || !first.cells) throw new Error('Field grid mismatch');
    const [a,b] = await Promise.all([this.field(first),this.field(second)]);
    let changed=0,sum=0,absolute=0,squared=0,maximum=0;
    for (let i=0;i<a.length;i++) {
      if (a[i]!==b[i]) changed++;
      if (!categorical) {
        const delta=b[i]-a[i], magnitude=Math.abs(delta);
        sum+=delta; absolute+=magnitude; squared+=delta*delta; maximum=Math.max(maximum,magnitude);
      }
    }
    return categorical?{changed_cells:changed}:{changed_cells:changed,mean_delta:sum/a.length,
      mae:absolute/a.length,rmse:Math.sqrt(squared/a.length),max_abs:maximum};
  }
}
if (typeof module!=='undefined') module.exports={MapFieldComparison};
