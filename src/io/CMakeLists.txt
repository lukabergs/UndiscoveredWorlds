# Run bundles, imports, TIFF and PNG, climate references, definitions.
add_library(uw_io STATIC
    climate_reference.cpp
    climate_tiff.cpp
    map_imports.cpp
    raster_image.cpp
    raster_imports.cpp
    run_bundle.cpp
    social_definitions.cpp
    tectonic_heightmap_io.cpp
    topography_tiff.cpp
)
target_link_libraries(uw_io PUBLIC uw_includes uw_tectonics PNG::PNG)
# Heightmap images and metric TIFFs of the engine tests.
uw_tectonics_source_defaults(tectonic_heightmap_io.cpp)
