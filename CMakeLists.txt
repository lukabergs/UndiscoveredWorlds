cmake_minimum_required(VERSION 3.15)
project(UndiscoveredWorlds)

# Use vcpkg for finding packages
set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")

# 1. Find the dependencies via vcpkg
find_package(SFML COMPONENTS graphics window system audio CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)
find_package(ImGuiFileDialog CONFIG REQUIRED)
find_package(imgui-sfml CONFIG REQUIRED)

# 2. Define your source files (based on your directory listing)
set(SOURCES
    main.cpp
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

# 4. Link everything together
# vcpkg handles the include paths and library paths automatically here
target_link_libraries(${PROJECT_NAME} PRIVATE 
    sfml-graphics
    sfml-window
    sfml-system
    sfml-audio
    imgui::imgui
    ImGui-SFML::ImGui-SFML
    ImGuiFileDialog::ImGuiFileDialog
)

# 5. Set C++ Standard
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
