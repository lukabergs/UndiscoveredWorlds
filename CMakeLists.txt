cmake_minimum_required(VERSION 3.15)
project(UndiscoveredWorlds LANGUAGES CXX)

include(CTest)
option(UW_BUILD_APP "Build the SFML/ImGui desktop application" ON)
option(UW_BUILD_TECTONICS_TOOLS "Build tectonics simulation, export, and inspection tools" OFF)

add_subdirectory(src/simulations)
add_subdirectory(src/simulations/geology/plate_tectonics)
if(BUILD_TESTING OR UW_BUILD_TECTONICS_TOOLS)
    add_subdirectory(tools/tectonics)
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
    find_package(Threads REQUIRED)

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

    # Desktop integration and world-level stages.
    set(SOURCES
        src/app/main.cpp
        src/app/app_environment.cpp
        src/app/create_world_screen.cpp
        src/wip/generation_workbench.cpp
        src/app/rendering/rendering.cpp
        src/app/diagnostics/world_generation_debug.cpp
        src/app/rendering/cuda_renderers_stub.cpp
        src/app/rendering/appearance_settings.cpp
        src/app/rendering/map_appearance.cpp
        src/pipeline/physical_layers.cpp
        src/pipeline/social_generation.cpp
        src/io/social_definitions.cpp
        src/core/grid_coordinates.cpp
        src/app/world_property_controls.cpp
        src/app/app_windows.cpp
        src/app/ui_panels.cpp
        src/app/ui_charts.cpp
        src/io/map_imports.cpp
        src/wip/classes.cpp
        src/validation/climate/climate_circulation_diagnostics.cpp
        src/pipeline/climate_generation.cpp
        src/wip/misc.cpp
        src/wip/planet.cpp
        src/app/diagnostics/profiler.cpp
        src/wip/region.cpp
        src/wip/regionalmap.cpp
    )

    set(UW_INCLUDE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/app"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/app/diagnostics"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/app/rendering"
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

    add_executable(${PROJECT_NAME} ${SOURCES})
    target_include_directories(${PROJECT_NAME} PRIVATE ${UW_INCLUDE_DIRS})
    set_target_properties(${PROJECT_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/bin"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin"
    )

    if(CMAKE_CUDA_COMPILER)
        target_sources(${PROJECT_NAME} PRIVATE src/app/rendering/cuda_renderers.cu)
        target_link_libraries(${PROJECT_NAME} PRIVATE CUDA::cudart)
        target_compile_definitions(${PROJECT_NAME} PRIVATE UW_ENABLE_CUDA_RENDERERS=1)
        set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED YES)
    endif()

    # Desktop and simulation modules
    target_link_libraries(${PROJECT_NAME} PRIVATE
        sfml-graphics
        sfml-window
        sfml-system
        sfml-audio
        Threads::Threads
        imgui::imgui
        ImGui-SFML::ImGui-SFML
        ImGuiFileDialog::ImGuiFileDialog
        plate_tectonics
        uw_climate
        uw_climate_io
        uw_generation_stages
        uw_geology
        uw_terrain
        uw_world_climate
        uw_hydrology
        uw_resources
        uw_society
    )

    set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

endif()

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

# Optional headless report-rendering accelerator; independent of the simulator.
add_library(climate_flow_texture SHARED src/validation/climate/flow_texture_kernel.cpp)
set_target_properties(climate_flow_texture PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
