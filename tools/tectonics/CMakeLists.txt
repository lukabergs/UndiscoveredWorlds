# PNG/ZLIB are already supplied by UW's SFML dependency closure.
find_package(PNG REQUIRED)

add_library(PlateTectonicsHeightmapIO STATIC heightmap_io.cpp)
target_include_directories(PlateTectonicsHeightmapIO PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(PlateTectonicsHeightmapIO PUBLIC plate_tectonics uw_topography_io PNG::PNG)
uw_tectonics_target_defaults(PlateTectonicsHeightmapIO)

add_library(PlateTectonicsPipelineSupport STATIC tectonic_pipeline_support.cpp)
target_include_directories(PlateTectonicsPipelineSupport PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(PlateTectonicsPipelineSupport PUBLIC plate_tectonics)
uw_tectonics_target_defaults(PlateTectonicsPipelineSupport)

if(UW_BUILD_TECTONICS_TOOLS)
    add_executable(tectonic_pipeline tectonic_pipeline.cpp)
    target_link_libraries(tectonic_pipeline PRIVATE PlateTectonicsPipelineSupport uw_mineral_tectonics)
    uw_tectonics_target_defaults(tectonic_pipeline)
endif()
