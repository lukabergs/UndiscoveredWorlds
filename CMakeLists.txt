cmake_minimum_required(VERSION 3.15)
project(UndiscoveredWorlds LANGUAGES CXX)

# Use vcpkg for finding packages
set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
endif()

# 1. Find the dependencies via vcpkg
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

# 2. Define your source files (based on your directory listing)
set(SOURCES
    main.cpp
    app_environment.cpp
    create_world_screen.cpp
    generation_stage_registry.cpp
    generation_workbench.cpp
    rendering.cpp
    world_generation_debug.cpp
    cuda_renderers_stub.cpp
    appearance_settings.cpp
    map_appearance.cpp
    plate_tectonics_adapter.cpp
    physical_layers.cpp
    social_generation.cpp
    social_sites.cpp
    social_routes.cpp
    social_polities.cpp
    social_trade.cpp
    world_property_controls.cpp
    app_windows.cpp
    ui_panels.cpp
    ui_charts.cpp
    map_imports.cpp
    assetdata.cpp
    classes.cpp
    climate_physics.cpp
    climate_reference.cpp
    climate_ocean.cpp
    fastlem_mountains.cpp
    globalclimate.cpp
    globalterrain.cpp
    misc.cpp
    planet.cpp
    profiler.cpp
    region.cpp
    regionalmap.cpp
)

set(PLATE_TECTONICS_ROOT "C:/dev/plate-tectonics" CACHE PATH "Path to the upstream plate-tectonics repository")
set(PLATE_TECTONICS_SRC_DIR "${PLATE_TECTONICS_ROOT}/src")

if(NOT EXISTS "${PLATE_TECTONICS_SRC_DIR}/platecapi.cpp")
    message(FATAL_ERROR "PLATE_TECTONICS_ROOT must point to the upstream plate-tectonics repo. Missing: ${PLATE_TECTONICS_SRC_DIR}/platecapi.cpp")
endif()

file(GLOB PLATE_TECTONICS_SOURCES CONFIGURE_DEPENDS "${PLATE_TECTONICS_SRC_DIR}/*.cpp")

add_library(plate_tectonics_vendor STATIC ${PLATE_TECTONICS_SOURCES})
target_include_directories(plate_tectonics_vendor PUBLIC "${PLATE_TECTONICS_SRC_DIR}")
target_compile_definitions(plate_tectonics_vendor PRIVATE NOMINMAX)
set_target_properties(plate_tectonics_vendor PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED YES)

# 3. Create the executable
add_executable(${PROJECT_NAME} ${SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/bin"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin"
)

if(CMAKE_CUDA_COMPILER)
    target_sources(${PROJECT_NAME} PRIVATE cuda_renderers.cu)
    target_link_libraries(${PROJECT_NAME} PRIVATE CUDA::cudart)
    target_compile_definitions(${PROJECT_NAME} PRIVATE UW_ENABLE_CUDA_RENDERERS=1)
    set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED YES)
endif()

# 4. Link everything together
# vcpkg handles the include paths and library paths automatically here
target_link_libraries(${PROJECT_NAME} PRIVATE 
    sfml-graphics
    sfml-window
    sfml-system
    sfml-audio
    Threads::Threads
    imgui::imgui
    ImGui-SFML::ImGui-SFML
    ImGuiFileDialog::ImGuiFileDialog
    plate_tectonics_vendor
)

# 5. Set C++ Standard
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

include(CTest)

if(BUILD_TESTING)
    add_executable(climate_physics_tests
        tests/climate_physics_tests.cpp
        climate_physics.cpp
    )
    target_include_directories(climate_physics_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    set_target_properties(climate_physics_tests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
    add_test(NAME climate_physics_tests COMMAND climate_physics_tests)

    add_executable(climate_reference_tests
        tests/climate_reference_tests.cpp
        climate_reference.cpp
    )
    target_include_directories(climate_reference_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    set_target_properties(climate_reference_tests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
    add_test(NAME climate_reference_tests COMMAND climate_reference_tests)
endif()
