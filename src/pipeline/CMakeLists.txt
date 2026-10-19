# Process registry, recipes, stages and checkpoints.
add_library(uw_pipeline STATIC
    climate_adapter.cpp
    climate_stages.cpp
    export_routing.cpp
    generation_progress.cpp
    generation_stage_registry.cpp
    physical_layers.cpp
    recipe_parameters.cpp
    social_generation.cpp
    world_generation_debug.cpp
    world_processes.cpp
    world_recipes.cpp
    world_stages.cpp
)
target_link_libraries(uw_pipeline PUBLIC uw_includes)
