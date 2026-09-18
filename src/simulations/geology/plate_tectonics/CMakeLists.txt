# UW's internal tectonics engine. Keep the source list explicit and repo-local.
option(PLATE_TECTONICS_WARNINGS_AS_ERRORS "Treat tectonics warnings as errors" ON)
option(PLATE_TECTONICS_ENABLE_ASAN "Enable AddressSanitizer for tectonics" OFF)
option(PLATE_TECTONICS_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer for tectonics" OFF)

function(uw_tectonics_target_defaults target)
    # The engine needs C++20; its public headers remain usable by UW's C++17 app.
    set_target_properties(${target} PROPERTIES
        CXX_STANDARD 20 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    foreach(config Debug Release RelWithDebInfo MinSizeRel)
        string(TOUPPER "${config}" config_upper)
        set_target_properties(${target} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/bin")
    endforeach()
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

add_library(plate_tectonics STATIC
    bounds.cpp
    geometry.cpp
    heightmap.cpp
    lithosphere.cpp
    mass.cpp
    movement.cpp
    noise.cpp
    plate.cpp
    plate_functions.cpp
    platecapi.cpp
    rectangle.cpp
    segment_creator.cpp
    segment_data.cpp
    segments.cpp
    simplerandom.cpp
    simplexnoise.cpp
    sqrdmd.cpp
    tectonic_contract.cpp
    tectonic_boundary_graph.cpp
    tectonic_deforming_regions.cpp
    tectonic_junction_coupling.cpp
    tectonic_junctions.cpp
    tectonic_provenance.cpp
    topography_codec.cpp
    topography_tiff.cpp
    utils.cpp
    world_point.cpp
)
target_include_directories(plate_tectonics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(plate_tectonics PRIVATE NOMINMAX)
uw_tectonics_target_defaults(plate_tectonics)
