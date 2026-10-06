"""Earth hydrology validation (erosion plan, Phase C, step C1): the routing of the earth-hydrology
recipe against the prepared HydroSHEDS fields of refs/processed/physical/fields/hydrology/ at the
run's width. The references validate; they never enter a simulation (world-model rule 4).

Kinds (each a function compare_<kind>(bundle, spec, baseline=None) returning
{'reference': text, 'run': {'metrics': {...}}, 'baseline': {...}}, as the other distribution kinds):

- upstream_area: drainage_area against rivers_upstream_area on the reference's river cells, best
  match within one cell; median absolute log10 ratio, share within a factor of 2, recall of the
  reference rivers above 1e5 km2;
- basins: the routing's drainage basins against basin_level3, area-weighted best-match overlap;
- internal_drainage: endorheic land of the routing against internal_drainage, overlap;
- discharge: discharge against rivers_discharge (a WaterGAP model estimate, not observations).

C0 declares the kinds; C1 implements them.
"""

KINDS = ('upstream_area', 'basins', 'internal_drainage', 'discharge')


def compare(kind, bundle, spec, baseline=None):
    raise NotImplementedError(f'validation kind {kind}: erosion plan C1')
