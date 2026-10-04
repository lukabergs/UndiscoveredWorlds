# Mineral systems

Gold is the first commodity whose deposits come from geological history rather
than present-day suitability scores. The chain is primary deposit formation,
transport with the host crust, exhumation and erosion of the ore body, and
river transport of liberated gold into placers. Stone, other metals, soils and
biological resources are not modelled yet.

## Data model

A `mineralsystems::Deposit` is a discrete ore body classified by the process
that formed it (`orogenic_gold`, `epithermal_gold`, `porphyry_copper_gold`,
`placer_gold`). Commodities are attributes: gold grade, and copper grade for
porphyries. Each deposit records:

- ore tonnage and grades: contained gold is `oreMt * goldGradeGt` tonnes;
- formation record: time, host crust class, nearest boundary segment, its type
  and age, the setting inferred for it (collision, subduction margin or
  transform), distance in km, setting intensity, inherited province;
- vertical geometry: formation depth of the ore-body top, vertical extent,
  exhumation since formation, eroded fraction;
- fate: preserved, eroded away or subducted, with time;
- gold liberated recently enough to remain in rivers;
- present surface conditions and, for placers, contributing source deposits.

Quantities are geological endowment assumptions, not reserves. Discovery,
extraction technology and value belong to society. `describedeposit()` gives
the plain-language explanation shown by the app's tectonics point inspector.

## Coupling to tectonics

`uw_mineral_systems` consumes a `TectonicFrame` after every physical update:
native ownership, topology identities, crust class/age, elevation, the
continental crust of the isostatic columns, boundaries and the displacement
supplied to the conservative column solver. `uw_mineral_tectonics` fills it
through `tectonic_api.hpp` and the engine's column model. The app drives history
from `generateworldtectonics`. `world_pipeline run --recipe tectonics` runs the same history (default;
`--mineral-history off` disables it) and writes the tables
`geosphere.mineral_deposits` (CSV: one row per deposit, enum codes and world
fractions) and `geosphere.mineral_budget` (JSON, `mineral-budget/v1`) at the end of the run.
Both use the same in-repo engine and baseline recipe, and the same seed stream.

Markers approximate point movement in the native displacement field. Ownership
changes re-home them; nearby host-crust snapping and crust-age comparisons are
heuristics. They do not reproduce conservative subcell mixing, track individual
rock parcels or establish exact subduction fate. There is no rigid-motion
fallback or administrative restart. Geological identity tracking needs further
validation before quantitative resource claims.

## Formation

For each update and native cell, two setting intensities are computed from the
nearest boundary segment, using cell-centre distances in km so that kernel
integrals do not depend on grid resolution:

| Setting | Boundary | Weight | Kernel |
|---|---|---|---|
| Orogenic | collision (both plates) / subduction (overriding plate) / transform (both) | 1.0 / 0.4 / 0.25 | Gaussian, 300 km from the boundary |
| Arc (epithermal, porphyry) | subduction (overriding plate) | 1.0 | Gaussian centred 150 km inboard, width 100 km |

Convergent segments are classified from the column crust within one cell of
the boundary on each plate. As in the column solver, continental crust does not
subduct and older ocean floor descends beneath younger: when the mean
continental crust of both edges reaches 20 km, the plates collide; otherwise the
plate with continental crust, or else the one with younger ocean floor,
overrides. The engine's segment regimes are not used. They derive from crust
classes that label much deep ocean floor continental, so Andean-type margins
were labelled collisions and only intra-oceanic margins were labelled arcs.

Activity is the segment's mean convergence (shear for transforms) divided by
50, capped at 1. The engine scores motion 0-100 over 0-5 cm/yr, so full
activity is 2.5 cm/yr. Orogenic gold also needs sustained deformation: maturity
rises smoothly from 5 to 20 Myr of segment age. In seed 1, one-half to
three-quarters of convergent boundary length is older than 20 Myr. Host crust
must stand above -1000 m: continents, arcs and accreted terranes. Deep ocean
floor hosts no deposits.

Expected deposits per update are `rate * area * dt * intensity`, with
equal native planar cell areas (`cellKm * cellKm`). Counts are
Poisson and locations are area-weighted, so refining the grid preserves the
expected endowment. Streams are seeded from the world seed.

| System | Rate (per Mkm² Myr) | Ore median (sigma ln) | Gold g/t | Top depth | Extent |
|---|---|---|---|---|---|
| Orogenic | 2.0 | 2 Mt (1.6) | 5 (0.6) | 3-15 km | 0.8-2.5 km |
| Epithermal | 3.0 | 1 Mt (1.5) | 5 (0.8) | 0.1-1 km | 0.2-0.7 km |
| Porphyry | 0.5 | 200 Mt (1.2) | 0.3 (0.7), Cu 0.45% (0.35) | 1.5-4 km | 0.8-2 km |

The grade-tonnage medians are order-of-magnitude choices after USGS
descriptive and grade-tonnage models (Cox and Singer 1986; Singer, Berger and
Moring 2008), truncated at 4 sigma. They are calibration inputs, not validated
statistics. Arc systems form about 75 t of gold per Mkm² of fully active arc
per Myr, comparable with a rough Cenozoic Earth estimate (about 40 kt over
66 Myr on about 8 Mkm² of arcs). That basis counts known deposits and is a lower
bound: the USGS estimates at least 3,100 Mt of undiscovered porphyry copper
against 1,800 Mt identified (Fact Sheet 2014-3004).

Inherited provinces represent history before the run: 0.4 provinces per Mkm²
of initial host crust, each an ellipse 150-400 km by 30-80 km (half
axes) holding on average 25 orogenic deposits dated 300-2500 Myr before the
run, with ore tops 0-3 km deep.

## Exhumation and liberation

Subaerial denudation follows `5 m/Myr * (elevation / 1 km)^3`, capped at
1 km/Myr: about 1-5 m/Myr on cratons, 40 m/Myr at 2 km and 320 m/Myr at 4 km.
Relief controls erosion, and relief grows with elevation. The law is a
reference-climate proxy until tectonics or the climate-erosion experiment
supplies per-cell erosion. Submerged crust does not erode.

Rock removed above a deposit is exhumation; once it reaches the ore body the
ore is removed uniformly through its vertical extent and the gold is
liberated. Liberated gold decays with a 5 Myr residence time, so the stock at
the end of the run is the gold still expected in the fluvial system.

## Present-day stage

`physical_layers::generate_gold_deposits` runs during finalization after
drainage basins and surface potentials. It records surface conditions at every
deposit (land or sea, elevation, 3x3 relief, annual rainfall, ice) and routes
each deposit's recently liberated gold along the final drainage network. Each
cell traps `1 - exp(-step / L)` of the passing gold. `L` is 25 km in ordinary
channels, shortened by up to 75% on depositional reaches and halved on deltas.
Lakes retain all gold; the sea receives the remainder. Cells with at least
0.5 t become placer deposits naming their largest sources; smaller amounts are
counted as dispersed. The budget satisfies
`injected = placer deposits + dispersed + lakes + exported`.

The existing 0..100 metal-ore and placer scores are unchanged, and society
still reads them. The Resources view adds gold tints for preserved lodes within
100 m of the surface and for placers.

## Verification and remaining work

`MineralSystemsTests` checks formation sampling, marker displacement and wrapping,
ownership changes, deposit/placer gold budgets, repeated-frame handling, and
ingestion of actual column displacements and deterministic topology identities
from a seed-1 native run. It does not establish exact material identity or
physical validity of inferred deposits.

The denudation model reads current engine elevation without subtracting a
private isostatic-lowering estimate. It remains an independent proxy, not the
solver's per-cell erosion budget. On seed 1, 15-18% of the engine's land stands
above 2 km after 300 Myr, and the proxy then removes two to four times the
volume that the column solver erodes. Deposit erosion is therefore high where the
engine keeps high terrain; per-cell erosion from the solver would remove that
inconsistency. Formation rates and inherited provinces remain order-of-magnitude
choices. Direct conserved material tracers, erosion histories, climate coupling
and biological resources remain future work.
