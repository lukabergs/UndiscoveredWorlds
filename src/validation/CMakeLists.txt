# Comparisons with observed data, and the snapshot record of the tectonics
# engine (the world writes it into its bundles; the engine tests read bundles).
add_library(uw_validation STATIC
    climate/climate_circulation_diagnostics.cpp
    climate/climate_validation.cpp
    tectonics/tectonic_snapshots.cpp
)
target_link_libraries(uw_validation PUBLIC uw_includes uw_tectonics)
uw_tectonics_source_defaults(tectonics/tectonic_snapshots.cpp)

# Optional headless report-rendering accelerator; independent of the simulator.
add_library(climate_flow_texture SHARED climate/flow_texture_kernel.cpp)
