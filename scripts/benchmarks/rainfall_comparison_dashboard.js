/* Numerical rainfall distributions for local benchmark galleries. No external services. */
(() => {
  const dataElement = document.getElementById('rain-data');
  if (!dataElement) return;
  const data = JSON.parse(dataElement.textContent);
  const node = id => document.getElementById(id);
  const colors = ['#f5ad55', '#60d6b0', '#aea5ff'];
  const svgElement = (tag, attributes, text) => {
    const element = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const [name, value] of Object.entries(attributes)) element.setAttribute(name, value);
    if (text !== undefined) element.textContent = text;
    return element;
  };
  const add = (parent, tag, attributes, text) => parent.appendChild(svgElement(tag, attributes, text));
  const format = value => Number(value).toFixed(2);
  function frame(svg, low, high, unit, labels) {
    svg.replaceChildren();
    svg.setAttribute('viewBox', '0 0 1000 350');
    svg.setAttribute('role', 'img');
    svg.setAttribute('aria-label', labels.join(', ') + '; ' + unit);
    const y = value => 275 - (value - low) / Math.max(1e-12, high - low) * 220;
    for (let i = 0; i <= 5; ++i) {
      const value = low + (high - low) * i / 5;
      add(svg, 'line', {x1: 75, x2: 970, y1: y(value), y2: y(value), stroke: '#394552'});
      add(svg, 'text', {x: 65, y: y(value) + 4, 'text-anchor': 'end', fill: '#edf2f7', 'font-size': 12}, format(value));
    }
    labels.forEach((label, i) => add(svg, 'text', {x: 75 + i * 300, y: 20, fill: colors[i], 'font-size': 14}, label));
    add(svg, 'text', {x: 75, y: 42, fill: '#c5d2df', 'font-size': 12}, unit);
    return y;
  }
  function update() {
    const selected = node('case').value;
    const control = node('rain-control').value === 'map' ? node('control').value : node('rain-control').value;
    if (!data.cases[selected] || !data.cases[control]) return;
    const mode = node('rain-period').value;
    const season = mode === 'season' ? node('season').value : mode;
    const region = node('rain-region').value;
    const a = data.cases[selected].seasons[season];
    const b = data.cases[control].seasons[season];
    const native = season === 'annual_native';
    const labels = [`Run ${data.cases[selected].run}`, native ? 'IMERG annual reference' : 'ERA5 reference', `Control ${data.cases[control].run}`];
    const histograms = [a.histograms[region].model, a.histograms[region].reference, b.histograms[region].model];
    const bins = native ? data.annual_bin_edges : data.bin_edges;
    const units = native ? 'mm/year' : 'mm/day';
    const svg = node('rain-histogram');
    const maximum = Math.max(...histograms.flatMap(h => h.histogram_area_percent)) * 1.1;
    const y = frame(svg, 0, maximum || 1, `Percent of matched area per bin; rainfall in ${units}`, labels);
    const width = 890 / (bins.length - 1);
    for (let i = 0; i < bins.length - 1; ++i) {
      const label = bins[i + 1] === null ? `${bins[i]}+` : `${bins[i]}–${bins[i + 1]}`;
      histograms.forEach((h, series) => {
        const value = h.histogram_area_percent[i];
        const bar = add(svg, 'rect', {x: 77 + i * width + series * width * .28,
          y: y(value), width: width * .25, height: 275 - y(value), fill: colors[series]});
        add(bar, 'title', {}, `${labels[series]}: ${label} ${units}; ${format(value)}% of area`);
      });
      add(svg, 'text', {x: 78 + i * width, y: 294, transform: `rotate(32 ${78 + i * width} 294)`, fill: '#edf2f7', 'font-size': 10}, label);
    }
    const thresholds = native ? ['1000to2500', '4500plus', '6000plus'] : ['1to6', '6to12', '12plus'];
    const thresholdLabels = native ? ['1,000–2,500', '≥4,500', '≥6,000'] : ['1–6', '6–12', '≥12'];
    const table = node('rain-statistics');table.replaceChildren();
    const header = document.createElement('tr');
    ['Dataset', `P50 ${units}`, `P90 ${units}`, `P99 ${units}`, ...thresholdLabels.map(x => `Area % ${x} ${units}`)].forEach(text => {
      const cell = document.createElement('th');cell.textContent = text;header.appendChild(cell);
    });table.appendChild(header);
    histograms.forEach((h, i) => {
      const row = document.createElement('tr');
      [labels[i], ...[.5, .9, .99].map(q => format(h.quantiles[String(q)])), ...thresholds.map(key => format(h.area_percent[key]))].forEach(text => {
        const cell = document.createElement('td');cell.textContent = text;row.appendChild(cell);
      });table.appendChild(row);
    });
    const sector = node('rain-sector').value;
    const field = node('rain-profile-field').value;
    const curves = [a.profiles[sector][field].model, a.profiles[sector][field].reference, b.profiles[sector][field].model];
    const profileUnits = {rain: units, sst: '°C', water: 'kg/m²', convergence: 'mm/day'}[field];
    const values = curves.flatMap(c => c.values);
    const minimum = Math.min(0, ...values), highest = Math.max(...values);
    const pad = Math.max(.01, (highest - minimum) * .08);
    const plot = node('rain-profile');
    const profileLabels = field === 'rain' ? labels : [labels[0], 'ERA5 reference', labels[2]];
    const py = frame(plot, minimum - pad, highest + pad, `${field}: ${profileUnits}; ocean-only latitude profiles`, profileLabels);
    const px = lat => 75 + (lat + 20) / 40 * 895;
    for (const lat of [-20, -10, 0, 10, 20]) {
      add(plot, 'text', {x: px(lat), y: 300, 'text-anchor': 'middle', fill: '#edf2f7', 'font-size': 12}, lat === 0 ? 'Equator' : `${Math.abs(lat)}°${lat < 0 ? 'S' : 'N'}`);
    }
    add(plot, 'line', {x1: px(0), x2: px(0), y1: 55, y2: 275, stroke: '#a1b0bd', 'stroke-dasharray': '5 5'});
    curves.forEach((curve, i) => add(plot, 'polyline', {points: curve.latitude.map((lat, j) => `${px(lat)},${py(curve.values[j])}`).join(' '),
      fill: 'none', stroke: colors[i], 'stroke-width': 2}));
    node('rain-profile-note').textContent = curves.map((curve, i) => `${profileLabels[i]}: weaker flank minus equator ${format(curve.bands.weaker_flank_minus_equator)} ${profileUnits}`).join(' · ') +
      '. For rainfall, a positive value means both flanks are wetter than the equatorial strip; compare the reference too. Sector averages can hide individual longitude features.';
    node('rain-definition').textContent = native ?
      'Annual histogram: exact numeric GeoTIFF behind the native model map versus its IMERG annual reference, on a common 1024×512 grid. The other profile fields use annual averages of the seasonal diagnostics. No color counting or histogram rescaling is applied.' :
      'Quarterly numerical rainfall versus ERA5 on the common 256×128 grid. Annual diagnostic mode combines quarterly means with 90/91/92/92-day weights; it is separate from the stored native annual map. No color counting or histogram rescaling is applied.';
  }
  for (const id of ['case', 'control', 'season', 'rain-control', 'rain-period', 'rain-region', 'rain-sector', 'rain-profile-field']) node(id).addEventListener('change', update);
  update();
})();
