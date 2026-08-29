#include <windows.h>
#include <shellapi.h>

#include "create_world_screen.hpp"

CreateWorldScreenActions drawcreateworldscreen(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, float currentversion, float latestversion)
{
    CreateWorldScreenActions actions;

    if (latestversion > currentversion)
    {
        ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 10, main_viewport->WorkPos.y + 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(224, 107), ImGuiCond_FirstUseEver);
        ImGui::Begin("Update available!", NULL, window_flags);
        ImGui::Text("Click here to visit the website.");
        ImGui::Text(" ");
        ImGui::Text(" ");
        ImGui::SameLine(95.0f);

        if (ImGui::Button("Go"))
            ShellExecuteW(0, 0, L"https://undiscoveredworlds.blogspot.com/2019/01/what-is-undiscovered-worlds.html", 0, 0, SW_SHOW);

        ImGui::End();
    }

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
