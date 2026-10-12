cmake_minimum_required(VERSION 3.19)

# Retain nvcc's existing default when no CUDA architecture is supplied.
cmake_policy(SET CMP0104 OLD)
project(UndiscoveredWorlds LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
# Multi-config generators append the configuration: bin/Release, bin/Debug.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

include(CTest)
option(UW_BUILD_APP "Build the SFML/ImGui desktop application" ON)
option(UW_BUILD_WORLD "Build world_pipeline" ON)
# Analysis builds only: world accessors record which fields each step reads.
option(UW_TRACK_FIELD_READS "Record world field reads per generation step" OFF)
if(UW_TRACK_FIELD_READS)
    add_compile_definitions(UW_TRACK_FIELD_READS)
endif()

find_package(PNG REQUIRED)
find_package(Threads REQUIRED)

# One static library per source directory, named uw_<directory>, each defined in
# the CMakeLists.txt of its directory. The libraries do not link each other,
# because world, pipeline and the simulation modules refer to each other; an
# executable links one of the two sets below.
# Every source includes project headers by their path from src/.
add_library(uw_includes INTERFACE)
target_include_directories(uw_includes INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_subdirectory(src)

# Climate on a climateworld, without planet: the climate processes, their stage
# runner, validation, imports, progress and environment. world_tests links this
# set alone, so a climate source that needs the world class fails its link.
add_library(uw_climate_world INTERFACE)
target_link_libraries(uw_climate_world INTERFACE
    uw_pipeline uw_validation uw_climate uw_io uw_core PNG::PNG Threads::Threads)

# The headless world: everything but the desktop application.
add_library(uw_headless INTERFACE)
target_link_libraries(uw_headless INTERFACE
    uw_cli uw_pipeline uw_world uw_validation uw_io uw_core
    uw_climate uw_geology uw_terrain uw_hydrology uw_landscape uw_regional uw_resources uw_society
    uw_mineral_tectonics uw_mineral_systems uw_tectonics
    PNG::PNG Threads::Threads)

if(UW_BUILD_WORLD OR UW_BUILD_APP)
    add_executable(world_pipeline src/cli/world_pipeline_main.cpp)
    target_link_libraries(world_pipeline PRIVATE uw_headless)
    # The physics of the terrain sandbox (scripts/sandbox/terrain_sandbox.py): the landscape kernels only.
    add_executable(terrain_sandbox src/cli/terrain_sandbox_main.cpp)
    target_link_libraries(terrain_sandbox PRIVATE uw_landscape uw_includes)
endif()

if(UW_BUILD_APP)
    # A language is enabled in the top directory.
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        find_package(CUDAToolkit REQUIRED)
    endif()
    add_subdirectory(src/app)
endif()

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
