add_subdirectory(simulations/geology/plate_tectonics)

function(uw_add_test name source module)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE ${module})
    set_target_properties(${name} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

uw_add_test(climate_grid_tests core/climate_grid_tests.cpp uw_grid)
uw_add_test(climate_palette_tests validation/climate_palette_tests.cpp uw_climate_io)
uw_add_test(parallel_rows_tests core/parallel_rows_tests.cpp uw_grid)
set_tests_properties(parallel_rows_tests PROPERTIES TIMEOUT 15)
uw_add_test(climate_pipeline_tests pipeline/climate_pipeline_tests.cpp uw_generation_stages)

foreach(component atmosphere boundary_layer context energy flow hydrology transport koppen ocean_dynamics ocean_heat physics weather)
    uw_add_test(climate_${component}_tests
        simulations/climate/climate_${component}_tests.cpp uw_climate)
endforeach()

foreach(component benchmark_outputs reference tiff)
    uw_add_test(climate_${component}_tests io/climate_${component}_tests.cpp uw_climate_io)
endforeach()
