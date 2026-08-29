#include <iomanip>
#include <sstream>

#include <windows.h>
#include <shellapi.h>

#include "app_windows.hpp"

using namespace std;

namespace
{
void settooltipifhovered(const char* text)
{
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", text);
}
}

bool drawcustomworldsizewindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, bool brandnew, WorldPropertyControls& controls)
{
    if (!show)
        return false;

    bool confirmed = false;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 427, main_viewport->WorkPos.y + 174), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(360, 220), ImGuiCond_FirstUseEver);
    ImGui::Begin("Create custom world?", NULL, window_flags);

    string introtext = "This will delete the current world.";

    if (brandnew)
        introtext = "Please select a size for the new world.";

    ImGui::Text(introtext.c_str());
    ImGui::Text(" ");

    const char* sizeitems[] = { "Small", "Medium", "Large" };
    if (ImGui::Combo("World size", &controls.size, sizeitems, IM_ARRAYSIZE(sizeitems)))
    {
        controls.mapWidth = getdefaultmapwidthforsize(controls.size);
        controls.mapHeight = getdefaultmapheightforsize(controls.size);
    }

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Earth: large; Mars: medium; Moon: small.");

    ImGui::InputInt("Map width", &controls.mapWidth);

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Actual world width in pixels for imports and generation.");

    ImGui::InputInt("Map height", &controls.mapHeight);

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Actual world height in pixels for imports and generation.");

    ImGui::Text(" ");
    ImGui::Text(" ");

    ImGui::SameLine(100.0f);

    if (ImGui::Button("OK"))
    {
        confirmed = true;
        show = false;
    }

    ImGui::SameLine(180.0f);

    if (ImGui::Button("Cancel"))
        show = false;

    ImGui::End();
    clampworldpropertycontrols(controls);
    return confirmed;
}

bool drawworldgenerationoptionswindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, WorldPropertyControls& controls)
{
    if (!show)
        return false;

    bool confirmed = false;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 398, main_viewport->WorkPos.y + 142), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(434, 208), ImGuiCond_FirstUseEver);
    ImGui::Begin("New world", NULL, window_flags);

    constexpr int presetcount = 7;
    const int presetwidths[presetcount] = { 32, 64, 128, 256, 512, 1024, 2048 };
    const int presetheights[presetcount] = { 16, 32, 64, 128, 256, 512, 1024 };
    const char* presetlabels[presetcount] =
    {
        "32 x 16",
        "64 x 32",
        "128 x 64",
        "256 x 128",
        "512 x 256",
        "1024 x 512",
        "2048 x 1024"
    };
    int selectedpreset = -1;
    for (int index = 0; index < presetcount; index++)
    {
        if (controls.mapWidth == presetwidths[index] && controls.mapHeight == presetheights[index])
        {
            selectedpreset = index;
            break;
        }
    }

    const char* previewlabel = selectedpreset >= 0 ? presetlabels[selectedpreset] : "Custom";
    if (ImGui::BeginCombo("World size", previewlabel))
    {
        for (int index = 0; index < presetcount; index++)
        {
            const bool selected = index == selectedpreset;
            if (ImGui::Selectable(presetlabels[index], selected))
            {
                controls.mapWidth = presetwidths[index];
                controls.mapHeight = presetheights[index];

                if (controls.mapWidth <= 128)
                    controls.size = 0;
                else if (controls.mapWidth <= 512)
                    controls.size = 1;
                else
                    controls.size = 2;
            }

            if (selected)
                ImGui::SetItemDefaultFocus();
        }

        ImGui::EndCombo();
    }
    settooltipifhovered("Choose a default world size, or override the dimensions below.");

    ImGui::InputInt("Map width", &controls.mapWidth);
    settooltipifhovered("World width in pixels for terrain generation and simulation.");

    ImGui::InputInt("Map height", &controls.mapHeight);
    settooltipifhovered("World height in pixels for terrain generation and simulation.");

    ImGui::Text(" ");
    ImGui::Text(" ");

    ImGui::SameLine(110.0f);

    if (ImGui::Button("OK"))
    {
        confirmed = true;
        show = false;
    }

    settooltipifhovered("Open the workbench with these world dimensions.");

    ImGui::SameLine(220.0f);

    if (ImGui::Button("Cancel"))
        show = false;

    settooltipifhovered("Return to the previous menu.");

    ImGui::End();
    clampworldpropertycontrols(controls);
    return confirmed;
}

bool drawterrainchooserwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show)
{
    if (!show)
        return false;

    bool confirmed = false;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 428, main_viewport->WorkPos.y + 167), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(405, 150), ImGuiCond_FirstUseEver);
    ImGui::Begin("Generate terrain", NULL, window_flags);

    ImGui::Text("This will overwrite any existing terrain.");
    ImGui::Text(" ");
    ImGui::TextWrapped("Initial plate tectonics now comes directly from the upstream plate-tectonics simulation seed and parameters.");
    ImGui::Text(" ");

    ImGui::SameLine(100.0f);

    if (ImGui::Button("OK"))
    {
        confirmed = true;
        show = false;
    }

    ImGui::SameLine(270.0f);

    if (ImGui::Button("Cancel"))
        show = false;

    ImGui::End();
    return confirmed;
}

bool drawareawarningwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show)
{
    if (!show)
        return false;

    bool confirmed = false;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 492, main_viewport->WorkPos.y + 156), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(218, 174), ImGuiCond_FirstUseEver);
    ImGui::Begin("Warning!##areawarning", NULL, window_flags);

    ImGui::Text("This is a very large area,");
    ImGui::Text("and may crash the program.");
    ImGui::Text(" ");
    ImGui::Text("Proceed?");
    ImGui::Text(" ");
    ImGui::Text(" ");

    ImGui::SameLine(40.0f);

    if (ImGui::Button("OK"))
    {
        confirmed = true;
        show = false;
    }

    ImGui::SameLine(130.0f);

    if (ImGui::Button("Cancel"))
        show = false;

    ImGui::End();
    return confirmed;
}

void drawaboutwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, float currentversion)
{
    if (!show)
        return;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 501, main_viewport->WorkPos.y + 111), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(343, 373), ImGuiCond_FirstUseEver);

    stringstream ss;
    ss << fixed << setprecision(2) << currentversion;

    const string title = "Undiscovered Worlds version " + ss.str();

    ImGui::Begin(title.c_str(), NULL, window_flags);
    ImGui::Text("Undiscovered Worlds is designed and written");
    ImGui::Text("in inelegant C++ by Jonathan Hill. The");
    ImGui::Text("interface uses Dear ImGUI and SFML.");
    ImGui::Text(" ");
    ImGui::Text("Special thanks to Frank Gennari for testing,");
    ImGui::Text("debugging, and advice, and also to");
    ImGui::Text("u/Iron-Phoenix2307 for the application icon.");
    ImGui::Text(" ");
    ImGui::Text("For more information and instructions for use,");
    ImGui::Text("please visit the website.");
    ImGui::Text(" ");
    ImGui::Text("The source code for this application is available");
    ImGui::Text("under the GNU General Public License.");
    ImGui::Text(" ");

    if (ImGui::Button("Website"))
        ShellExecuteW(0, 0, L"https://undiscoveredworlds.blogspot.com/2019/01/what-is-undiscovered-worlds.html", 0, 0, SW_SHOW);

    ImGui::SameLine(140.0f);

    if (ImGui::Button("Source"))
        ShellExecuteW(0, 0, L"https://github.com/JonathanCRH/Undiscovered_Worlds", 0, 0, SW_SHOW);

    ImGui::SameLine(280.0f);

    if (ImGui::Button("Close"))
        show = false;

    ImGui::End();
}

void drawworldeditpropertieswindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, WorldPropertyControls& controls, const string& degree, bool& show)
{
    if (!show)
        return;

    const int rightalign = 280;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 420, main_viewport->WorkPos.y + 50), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(415, 427), ImGuiCond_FirstUseEver);
    ImGui::Begin("World properties##2", NULL, window_flags);

    ImGui::PushItemWidth(180);

    const char* rotationitems[] = { "East to west", "West to east" };
    ImGui::Combo("Rotation", &controls.rotation, rotationitems, IM_ARRAYSIZE(rotationitems));

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Affects weather patterns. (Earth: west to east)");

    const char* perihelionitems[] = { "January", "July" };
    ImGui::Combo("Perihelion", &controls.perihelion, perihelionitems, IM_ARRAYSIZE(perihelionitems));

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("When the planet is closest to the sun. (Earth: January)");

    ImGui::InputFloat("Eccentricity", &controls.eccentricity, 0.01f, 1.0f, "%.3f");

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("How elliptical the orbit is. (Earth: 0.0167)");

    ImGui::InputFloat("Obliquity", &controls.tilt, 0.01f, 1.0f, "%.3f");

    const string tilttip = "Affects seasonal variation in temperature. (Earth: 22.5" + degree + ")";

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip(tilttip.c_str());

    ImGui::InputFloat("Surface gravity", &controls.gravity, 0.01f, 1.0f, "%.3f");

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Affects mountain and valley sizes. (Earth: 1.00g)");

    ImGui::InputFloat("Lunar pull", &controls.lunar, 0.01f, 1.0f, "%.3f");

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Affects tides and coastal regions. (Earth: 1.00)");

    ImGui::InputFloat("Moisture pickup rate", &controls.waterpickup, 0.01f, 1.0f, "%.3f");

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Affects how much moisture wind picks up from the ocean. (Earth: 1.00)");

    ImGui::InputFloat("Heat decrease per vertical km", &controls.tempdecrease, 0.01f, 1.0f, "%.3f");

    const string tempdecreasetip = "Affects how much colder it gets higher up. (Earth: 6.5" + degree + ")";

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip(tempdecreasetip.c_str());

    ImGui::InputInt("Glaciation temperature", &controls.glacialtemp);

    const string glacialtip = "Areas below this average temperature may show signs of past glaciation. (Earth: 4" + degree + ")";

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip(glacialtip.c_str());

    ImGui::InputInt("Average global temperature", &controls.averagetemp);

    const string avetip = "Earth: 14" + degree;

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip(avetip.c_str());

    ImGui::InputInt("North pole adjustment", &controls.northpolaradjust);

    const string northtip = "Adjustment to north pole temperature. (Earth: +3" + degree + ")";

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip(northtip.c_str());

    ImGui::InputInt("South pole adjustment", &controls.southpolaradjust);

    const string southtip = "Adjustment to south pole temperature. (Earth: -3" + degree + ")";

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip(southtip.c_str());

    ImGui::Text("   ");

    ImGui::Checkbox("Generate rivers", &controls.rivers);

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Untick this if you don't want rivers to be generated.");

    ImGui::SameLine((float)rightalign);

    ImGui::Checkbox("Generate lakes", &controls.lakes);

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Untick this if you don't want large lakes to be generated.");

    ImGui::Checkbox("Generate deltas", &controls.deltas);

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Untick this if you don't want river deltas to be generated.");

    ImGui::Text("   ");
    ImGui::Checkbox("Generate social world", &controls.socialEnabled);

    const char* socialmodeitems[] = { "Static ex nihilo", "Historical" };
    int socialmode = controls.socialMode == SocialGenerationOptions::Mode::historical ? 1 : 0;

    ImGui::BeginDisabled(!controls.socialEnabled);
    ImGui::Combo("Social mode", &socialmode, socialmodeitems, IM_ARRAYSIZE(socialmodeitems));
    controls.socialMode = socialmode == 1 ? SocialGenerationOptions::Mode::historical : SocialGenerationOptions::Mode::static_ex_nihilo;

    ImGui::BeginDisabled(controls.socialMode != SocialGenerationOptions::Mode::historical);
    ImGui::Checkbox("Use prehistory", &controls.usePrehistory);
    ImGui::InputInt("History years", &controls.historyYears);
    ImGui::EndDisabled();
    ImGui::EndDisabled();

    ImGui::SameLine((float)rightalign);

    if (ImGui::Button("Close", ImVec2(120.0f, 0.0f)))
        show = false;

    ImGui::End();

    clampworldpropertycontrols(controls);
}
