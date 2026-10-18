# UW's internal tectonics engine. Keep the source list explicit and repo-local.
option(PLATE_TECTONICS_WARNINGS_AS_ERRORS "Treat tectonics warnings as errors" ON)
option(PLATE_TECTONICS_ENABLE_ASAN "Enable AddressSanitizer for tectonics" OFF)
option(PLATE_TECTONICS_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer for tectonics" OFF)

function(uw_tectonics_target_defaults target)
    set_target_properties(${target} PROPERTIES CXX_EXTENSIONS NO)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive- /EHsc)
        if(PLATE_TECTONICS_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
        if(PLATE_TECTONICS_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()
    if(PLATE_TECTONICS_ENABLE_ASAN)
        if(MSVC)
            target_compile_options(${target} PRIVATE /fsanitize=address)
            target_link_options(${target} PUBLIC /fsanitize=address)
        else()
            target_compile_options(${target} PRIVATE -fsanitize=address -fno-omit-frame-pointer)
            target_link_options(${target} PUBLIC -fsanitize=address)
        endif()
    endif()
    if(PLATE_TECTONICS_ENABLE_UBSAN AND NOT MSVC)
        target_compile_options(${target} PRIVATE -fsanitize=undefined -fno-omit-frame-pointer)
        target_link_options(${target} PUBLIC -fsanitize=undefined)
    endif()
endfunction()

# The same warning level for single sources of another library that work on
# engine types (snapshot record, heightmap I/O).
function(uw_tectonics_source_defaults)
    if(MSVC)
        set(options /W4 /permissive-)
        if(PLATE_TECTONICS_WARNINGS_AS_ERRORS)
            list(APPEND options /WX)
        endif()
    else()
        set(options -Wall -Wextra -Wpedantic)
        if(PLATE_TECTONICS_WARNINGS_AS_ERRORS)
            list(APPEND options -Werror)
        endif()
    endif()
    set_source_files_properties(${ARGN} PROPERTIES COMPILE_OPTIONS "${options}")
endfunction()

include(tectonic_recipe.cmake)
add_library(uw_tectonics STATIC
    lithosphere.cpp
    tectonic_columns.cpp
    tectonic_initialization.cpp
    tectonic_deformation.cpp
    tectonic_ocean_relief.cpp
    tectonic_noise_pass.cpp
    tectonic_plumes.cpp
    tectonic_margins.cpp
    tectonic_transport.cpp
    tectonic_sphere_transport.cpp
    tectonic_topology.cpp
    tectonic_grid.cpp
    tectonic_api.cpp
    tectonic_recipe.cpp
    simplerandom.cpp
    simplexnoise.cpp
    tectonic_contract.cpp
    tectonic_dynamics.cpp
    tectonic_boundary_graph.cpp
    tectonic_deforming_regions.cpp
    tectonic_junction_coupling.cpp
    tectonic_junctions.cpp
    tectonic_provenance.cpp
    topography_codec.cpp
    utils.cpp
)
target_link_libraries(uw_tectonics PUBLIC uw_includes uw_core PRIVATE uw_landscape)
target_include_directories(uw_tectonics PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_compile_definitions(uw_tectonics PRIVATE NOMINMAX)
uw_tectonics_target_defaults(uw_tectonics)
