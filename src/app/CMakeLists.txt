# Desktop application: rendering, windows, panels and the generation workbench.
# The only directory that uses SFML and ImGui.
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

add_library(uw_app STATIC
    app_windows.cpp
    appearance_files.cpp
    create_world_screen.cpp
    generation_workbench.cpp
    rendering/appearance_settings.cpp
    rendering/cuda_renderers_stub.cpp
    rendering/map_appearance.cpp
    rendering/rendering.cpp
    ui_charts.cpp
    ui_panels.cpp
    world_property_controls.cpp
)
target_link_libraries(uw_app PUBLIC
    uw_includes
    sfml-graphics
    sfml-window
    sfml-system
    sfml-audio
    imgui::imgui
    ImGui-SFML::ImGui-SFML
    ImGuiFileDialog::ImGuiFileDialog
)

if(CMAKE_CUDA_COMPILER)
    target_sources(uw_app PRIVATE rendering/cuda_renderers.cu)
    target_link_libraries(uw_app PUBLIC CUDA::cudart)
    target_compile_definitions(uw_app PUBLIC UW_ENABLE_CUDA_RENDERERS=1)
    set_target_properties(uw_app PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED YES)
endif()

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE uw_app uw_headless)
