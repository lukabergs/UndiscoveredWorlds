# The command line of world_pipeline. Its main() is the source of the executable.
add_library(uw_cli STATIC
    climate_options.cpp
    command_line_generation.cpp
    regional_command.cpp
    tectonic_options.cpp
)
target_link_libraries(uw_cli PUBLIC uw_includes uw_tectonics)
