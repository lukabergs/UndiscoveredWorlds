# Simulation boundaries

The refactor preserves active entry-point signatures, stage order, seed derivation,
and serialized fields. `functions.hpp` remains a compatibility umbrella for the
desktop and regional generator; it delegates simulation declarations to their
own headers. Module-local helpers are private rather than declared globally.

| Responsibility | Owner |
| --- | --- |
| Fractal fields, coasts, ridges, elevation, shelves, volcanoes, craters | `simulations/terrain/terrain_*.cpp` |
| Imported mountain masks and climate-shaped landforms | `simulations/terrain/` |
| Energy, circulation, ocean dynamics, moisture and weather kernels | `simulations/climate/`, target `uw_climate` |
| Planet climate fields, classification and solver coupling | `climate_fields`, `climate_classification`, `climate_coupling` |
| River routing, drainage basins, lakes, lake effects, deltas, wetlands, tides | `simulations/hydrology/` |
| Tectonic force balance, columns, topology and native outputs | `simulations/geology/plate_tectonics/`, target `plate_tectonics` |
| Shared recipe and world ingestion | `tectonic_recipe`, `geology/world_tectonics` |
| Mineral deposits formed and carried through tectonic history | `simulations/geology/mineral_systems/`, targets `uw_mineral_systems`, `uw_mineral_tectonics` |
| Derived mineral and marine potential scores; present-day gold placers | `simulations/resources/` |
| Settlement sites, routes, polities, trade and history | `simulations/society/` |
| Climate/import orchestration and physical/social stage order | `pipeline/` |
| Generation stage execution (formerly in the workbench) | `pipeline/world_stages` |
| Process kinds and declared field access; recipes, scheduler, write audit | `pipeline/world_processes`, `pipeline/world_recipes` |
| World field catalogue and read tracking | `core/world/` |
| Run bundles (fields and tables as NumPy files with a manifest) | `io/run_bundle`, `wip/planet_bundle`, `wip/planet_fields` |
| Headless PNG I/O for imports and validation | `io/raster_image` |
| Climate module world: climate state plus bound terrain and planet inputs | `simulations/climate/climate_world` |
| Climate stages on a climate world (shared by both pipelines) | `pipeline/climate_stages` |
| World-independent map decoding and land/sea imports | `io/raster_imports` |
| Seeded random numbers and fractal noise | `core/procedural_noise` |
| Authored society CSV parsing | `io/social_definitions` |
| Seasonal indices, seed derivation, coordinates, row partitioning | `core/` |

Physical layers execute in their original dependency order: drainage basins,
surface erosion/deposition/fertility, mineral resources, gold placers, marine
resources. The 0..100 resource scores remain heuristics, not mineral-formation
or ecology models.

Gold is the first commodity with a formation history (see
[mineral systems](mineral-systems.md)). `uw_mineral_systems` is engine-neutral:
it consumes a `TectonicFrame` view after each update, forms discrete deposits
from the setting of each native cell, carries them with their host plate's
native column displacement, and exhumes them under a denudation proxy. Point markers approximate material transport; they do not track conservative subcell mixing.
`uw_mineral_tectonics` fills frames from the plate-tectonics C API. The
world-generation entry point drives both during `generateworldtectonics`; the planet stores
the deposit list; finalization evaluates surface conditions and routes
liberated gold into placers on the final drainage network.
Climate atmospheric water transport is distinct from the surface routing in
`hydrology/`.

The standalone solver in `simulations/terrain/climate_erosion.*` (part of `uw_terrain`)
accepts numerical terrain and climate fields, conserves rock-equivalent volume on a
cell-centred spherical grid, and has no app or renderer dependency. Nothing calls it
at present: its CLI adapter and the feedback script that sequenced it (tectonics run
63) were removed, and it waits for a `terrain.climate_erosion` evolver. The app's climate-only import mode
selects climate stages through `pipeline/` and asserts that terrain is unchanged.
This is a surface experiment with frozen plates; the tectonic state is not advanced.

The climate-only import mode runs the `earth-climate` recipe: climate processes on
imported terrain with the Geosphere frozen. The recipe omits the roughness stage,
whose seeded noise climate never reads; land and sea counts are computed when the
imported terrain is prepared.

Each source directory is one static library named after it (`uw_core`, `uw_io`,
`uw_world`, `uw_pipeline`, `uw_validation`, `uw_cli`, `uw_app`, and `uw_climate`,
`uw_geology`, `uw_terrain`, `uw_hydrology`, `uw_resources`, `uw_society` under
`simulations/`), defined in the `CMakeLists.txt` of its directory; the engine is
`uw_tectonics`. The libraries do not link each other, because world, pipeline and the
modules refer to each other; the root `CMakeLists.txt` names two sets. `uw_headless`
is everything but the desktop and is what `world_pipeline` links. `uw_climate_world`
(core, I/O, climate, pipeline, validation) leaves out `uw_world` and the planet-level
modules: `world_tests` links only it, so a climate source that needs planet fails
that link. Planet-level modules still rely on world storage and legacy utilities in
`src/world/`, but none depends on SFML or ImGui; `world_pipeline` builds with
`UW_BUILD_APP=OFF`, and the desktop adds `uw_app` (rendering, panels, workbench).
Planet persistence, regional generation and the utility umbrella remain coupled;
see [world model](world-model.md).

## Retention decisions

The [retained climate alternatives](../../src/simulations/climate/legacy/README.md)
remain compiled and tested. Nonlocal thermal response was introduced in run 119,
disabled in 135, and replaced by local hypsometric forcing in 136. Empirical
SST-to-wind feedback remains available for controlled ocean comparisons through
`OceanConfig::oneWay = false`; production sets it to true. Both alternatives
have small named operations, with the original numerical contracts preserved.

Removed code had no live caller in the repository: the abandoned whole-terrain
generator and its continent-building subtree, old latitude-temperature and
directional rain/monsoon paths, prescribed ocean current/SST maps, sequential
river-drop routing, and orphan helpers. Disabled regional debugging/spring
checks and the disabled mountain-drainage direction branch were also removed.
Active editor mountain chains, ocean features, all-land rainfall, legacy field
synchronization and save compatibility code remain.

Evidence comes from [the run registry](../../runs/registry/climate.json) and local
`docs/tasks/climate-physics/FUTURE_WORK.md` notes. There is no retained benchmark
comparison establishing that the unused procedural continent generator is worse
than tectonics; its removal is a reachability decision. Early run 53's eddy
rejection does not justify deleting the active closure: run 68 restored it.
Run 159 is a coarse baseline, and its comparison with 158 changes the coupling
iteration cap, so it is not a controlled ranking of simulation methods.
