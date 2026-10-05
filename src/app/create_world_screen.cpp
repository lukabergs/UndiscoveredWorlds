#include "app/create_world_screen.hpp"

CreateWorldScreenActions drawcreateworldscreen(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags)
{
    CreateWorldScreenActions actions;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 483, main_viewport->WorkPos.y + 206), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(266, 78), ImGuiCond_FirstUseEver);
    ImGui::Begin("Create world", NULL, window_flags);

    if (ImGui::Button("New", ImVec2(74.0f, 0.0f)))
        actions.openPlateTectonicsMenu = true;

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Open the plate tectonics generation menu.");

    ImGui::SameLine();

    if (ImGui::Button("Custom", ImVec2(74.0f, 0.0f)))
        actions.openCustomWorld = true;

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Create a custom world.");

    ImGui::SameLine();

    if (ImGui::Button("Load", ImVec2(74.0f, 0.0f)))
        actions.openLoadDialog = true;

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Load a world.");

    ImGui::End();
    return actions;
}
