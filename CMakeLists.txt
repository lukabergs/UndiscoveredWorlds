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

# 2. Define your source files (based on your directory listing)
set(SOURCES
    main.cpp
    create_world_screen.cpp
    rendering.cpp
    cuda_renderers_stub.cpp
    appearance_settings.cpp
    world_property_controls.cpp
    app_windows.cpp
    ui_panels.cpp
    ui_charts.cpp
    map_imports.cpp
    assetdata.cpp
    classes.cpp
    globalclimate.cpp
    globalterrain.cpp
    misc.cpp
    planet.cpp
    profiler.cpp
    region.cpp
    regionalmap.cpp
)

# 3. Create the executable
add_executable(${PROJECT_NAME} ${SOURCES})

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
)

# 5. Set C++ Standard
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
