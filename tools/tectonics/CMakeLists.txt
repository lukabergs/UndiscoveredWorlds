# PNG/ZLIB are already supplied by UW's SFML dependency closure.
find_package(PNG REQUIRED)

add_library(PlateTectonicsHeightmapIO STATIC heightmap_io.cpp)
target_include_directories(PlateTectonicsHeightmapIO PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(PlateTectonicsHeightmapIO PUBLIC plate_tectonics PNG::PNG)
uw_tectonics_target_defaults(PlateTectonicsHeightmapIO)

add_library(PlateTectonicsPipelineSupport STATIC tectonic_pipeline_support.cpp)
target_include_directories(PlateTectonicsPipelineSupport PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(PlateTectonicsPipelineSupport PUBLIC plate_tectonics)
uw_tectonics_target_defaults(PlateTectonicsPipelineSupport)

if(UW_BUILD_TECTONICS_TOOLS)
    add_executable(tectonic_simulation simulation.cpp)
    target_compile_definitions(tectonic_simulation PRIVATE
        UW_TECTONICS_OUTPUT_ROOT="${CMAKE_BINARY_DIR}/tectonics")
    target_link_libraries(tectonic_simulation PRIVATE PlateTectonicsHeightmapIO ZLIB::ZLIB)
    uw_tectonics_target_defaults(tectonic_simulation)

    add_executable(contract_fixture contract_fixture.cpp)
    target_link_libraries(contract_fixture PRIVATE plate_tectonics)
    uw_tectonics_target_defaults(contract_fixture)

    add_executable(tectonic_pipeline tectonic_pipeline.cpp)
    target_link_libraries(tectonic_pipeline PRIVATE PlateTectonicsPipelineSupport)
    uw_tectonics_target_defaults(tectonic_pipeline)
endif()
