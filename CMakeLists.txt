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
option(UW_BUILD_WORLD "Build the headless world core and world_pipeline" ON)
# Analysis builds only: world accessors record which fields each step reads.
option(UW_TRACK_FIELD_READS "Record world field reads per generation step" OFF)
if(UW_TRACK_FIELD_READS)
    add_compile_definitions(UW_TRACK_FIELD_READS)
endif()

add_subdirectory(src/simulations)
add_subdirectory(src/simulations/geology/plate_tectonics)

# The snapshot record and snapshot bundles of the tectonics engine: the world
# writes the record into its bundles; the engine tests use the bundles.
# Statistics and comparisons of snapshots are Python (scripts/evaluation/).
add_library(uw_tectonics_snapshots STATIC src/validation/tectonics/tectonic_snapshots.cpp)
target_include_directories(uw_tectonics_snapshots PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/validation/tectonics")
target_link_libraries(uw_tectonics_snapshots PUBLIC plate_tectonics)
uw_tectonics_target_defaults(uw_tectonics_snapshots)

if(BUILD_TESTING)
    # Heightmap images and metric TIFFs of the engine tests.
    find_package(PNG REQUIRED)
    add_library(uw_tectonics_heightmap_io STATIC src/io/tectonic_heightmap_io.cpp)
    target_include_directories(uw_tectonics_heightmap_io PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/io")
    target_link_libraries(uw_tectonics_heightmap_io PUBLIC plate_tectonics uw_topography_io PNG::PNG)
    uw_tectonics_target_defaults(uw_tectonics_heightmap_io)
endif()

if(UW_BUILD_WORLD OR UW_BUILD_APP)
    find_package(PNG REQUIRED)
    find_package(Threads REQUIRED)

    # Headless world core: world state, generation stages, recipes, I/O and
    # validation. No SFML or ImGui; shared by world_pipeline and the desktop.
    set(UW_WORLD_CORE_SOURCES
        src/app/command_line_generation.cpp
        src/app/tectonic_options.cpp
        src/app/diagnostics/profiler.cpp
        src/core/grid_coordinates.cpp
        src/io/map_imports.cpp
        src/io/social_definitions.cpp
        src/pipeline/climate_generation.cpp
        src/pipeline/physical_layers.cpp
        src/pipeline/social_generation.cpp
        src/pipeline/world_recipes.cpp
        src/pipeline/world_stages.cpp
        src/wip/classes.cpp
        src/wip/misc.cpp
        src/wip/planet.cpp
        src/wip/planet_bundle.cpp
        src/wip/planet_fields.cpp
        src/wip/region.cpp
        src/wip/regionalmap.cpp
    )

    set(UW_INCLUDE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/app"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/app/diagnostics"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/core"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/io"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/pipeline"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/simulations/climate"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/simulations/geology"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/simulations/hydrology"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/simulations/society"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/simulations/terrain"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/validation/climate"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/wip"
    )

    include(src/simulations/world_modules.cmake)

    # Climate on a climateworld, without planet: the climate processes and their
    # stage runner, validation, map imports, progress and environment. Linked by
    # uw_world_core; climate_world_tests links it alone to keep it planet-free.
    add_library(uw_climate_world STATIC
        src/app/app_environment.cpp
        src/app/climate_options.cpp
        src/io/raster_image.cpp
        src/io/raster_imports.cpp
        src/pipeline/climate_stages.cpp
        src/pipeline/generation_progress.cpp
        src/pipeline/world_generation_debug.cpp
        src/pipeline/world_processes.cpp
        src/simulations/climate/climate_classification.cpp
        src/simulations/climate/climate_coupling.cpp
        src/simulations/climate/climate_fields.cpp
        src/validation/climate/climate_circulation_diagnostics.cpp
        src/validation/climate/climate_validation.cpp
    )
    target_include_directories(uw_climate_world PUBLIC ${UW_INCLUDE_DIRS})
    target_link_libraries(uw_climate_world PUBLIC
        Threads::Threads
        PNG::PNG
        uw_climate
        uw_climate_io
        uw_topography_io
        uw_world_fields
        uw_run_bundle
        uw_generation_stages
        uw_procedural_noise
    )

    add_library(uw_world_core OBJECT ${UW_WORLD_CORE_SOURCES})
    target_include_directories(uw_world_core PUBLIC ${UW_INCLUDE_DIRS})
    target_link_libraries(uw_world_core PUBLIC
        uw_climate_world
        Threads::Threads
        PNG::PNG
        plate_tectonics
        uw_mineral_tectonics
        uw_tectonics_snapshots
        uw_climate
        uw_climate_io
        uw_topography_io
        uw_world_fields
        uw_run_bundle
        uw_generation_stages
        uw_procedural_noise
    )

    # Object libraries contribute their objects only to direct dependents, so
    # each executable links the core and every world module.
    function(uw_add_world_executable name entry)
        add_executable(${name} ${entry})
        target_link_libraries(${name} PRIVATE
            uw_world_core
            uw_geology
            uw_terrain
            uw_hydrology
            uw_resources
            uw_society
        )
    endfunction()

    uw_add_world_executable(world_pipeline src/app/world_pipeline_main.cpp)
endif()

if(UW_BUILD_APP)
    # Use vcpkg for finding packages
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        find_package(CUDAToolkit REQUIRED)
    endif()

    # Desktop dependencies
    find_package(SFML COMPONENTS graphics window system audio CONFIG REQUIRED)
    find_package(imgui CONFIG REQUIRED)
    find_package(ImGuiFileDialog CONFIG REQUIRED)
    find_package(imgui-sfml CONFIG REQUIRED)

    foreach(dep_target
        imgui::imgui
        ImGui-SFML::ImGui-SFML
        ImGuiFileDialog::ImGuiFileDialog
        fmt::fmt
        fmt::fmt-header-only)
        if(TARGET ${dep_target})
            set_property(TARGET ${dep_target} PROPERTY INTERFACE_COMPILE_FEATURES "")
        endif()
    endforeach()

    # Desktop presentation: rendering, windows, panels and the generation workbench.
    set(UW_DESKTOP_SOURCES
        src/app/app_windows.cpp
        src/app/create_world_screen.cpp
        src/app/rendering/appearance_settings.cpp
        src/app/rendering/cuda_renderers_stub.cpp
        src/app/rendering/map_appearance.cpp
        src/app/rendering/rendering.cpp
        src/app/ui_charts.cpp
        src/app/ui_panels.cpp
        src/app/world_property_controls.cpp
        src/wip/generation_workbench.cpp
    )

    add_library(uw_app_core OBJECT ${UW_DESKTOP_SOURCES})
    target_include_directories(uw_app_core PUBLIC ${UW_INCLUDE_DIRS} "${CMAKE_CURRENT_SOURCE_DIR}/src/app/rendering")
    target_link_libraries(uw_app_core PUBLIC
        uw_world_core
        sfml-graphics
        sfml-window
        sfml-system
        sfml-audio
        imgui::imgui
        ImGui-SFML::ImGui-SFML
        ImGuiFileDialog::ImGuiFileDialog
    )

    if(CMAKE_CUDA_COMPILER)
        target_sources(uw_app_core PRIVATE src/app/rendering/cuda_renderers.cu)
        target_link_libraries(uw_app_core PUBLIC CUDA::cudart)
        target_compile_definitions(uw_app_core PUBLIC UW_ENABLE_CUDA_RENDERERS=1)
        set_target_properties(uw_app_core PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED YES)
    endif()

    uw_add_world_executable(${PROJECT_NAME} src/app/main.cpp)
    target_link_libraries(${PROJECT_NAME} PRIVATE uw_app_core)
endif()

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

# Optional headless report-rendering accelerator; independent of the simulator.
add_library(climate_flow_texture SHARED src/validation/climate/flow_texture_kernel.cpp)
