#include <windows.h>
#include <shellapi.h>

#include "create_world_screen.hpp"

CreateWorldScreenActions drawcreateworldscreen(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, float currentversion, float latestversion, bool brandnew, int& seedentry, int (*randomseedfactory)())
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
    ImGui::SetNextWindowSize(ImVec2(249, 90), ImGuiCond_FirstUseEver);
    ImGui::Begin("Create world", NULL, window_flags);

    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputInt(" ", &seedentry);

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Please enter a seed number, which will be used to calculate the new world. The same number will always yield the same world.");

    if (seedentry < 0)
        seedentry = 0;

    const char* loadlabel = brandnew ? "Load" : "Cancel";
    const char* loadtooltip = brandnew ? "Load a world." : "Return to the world map.";

    if (ImGui::Button(loadlabel))
    {
        if (brandnew)
            actions.openLoadDialog = true;
        else
            actions.returnToGlobalMap = true;
    }

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip(loadtooltip);

    ImGui::SameLine();

    if (ImGui::Button("Custom"))
        actions.openCustomWorld = true;

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Create a custom world.");

    ImGui::SameLine();

    if (ImGui::Button("Random"))
        seedentry = randomseedfactory();

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Roll a random seed number.");

    ImGui::SameLine();

    if (ImGui::Button("OK"))
        actions.startWorldGeneration = true;

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Create a world from this seed number.");

    ImGui::End();
    return actions;
}
