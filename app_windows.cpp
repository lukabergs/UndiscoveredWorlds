#include <iomanip>
#include <sstream>

#include <windows.h>
#include <shellapi.h>

#include "app_windows.hpp"

using namespace std;

bool drawcustomworldsizewindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, bool brandnew, int& currentsize)
{
    if (!show)
        return false;

    bool confirmed = false;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 427, main_viewport->WorkPos.y + 174), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(323, 152), ImGuiCond_FirstUseEver);
    ImGui::Begin("Create custom world?", NULL, window_flags);

    string introtext = "This will delete the current world.";

    if (brandnew)
        introtext = "Please select a size for the new world.";

    ImGui::Text(introtext.c_str());
    ImGui::Text(" ");

    const char* sizeitems[] = { "Small", "Medium", "Large" };
    ImGui::Combo("World size", &currentsize, sizeitems, IM_ARRAYSIZE(sizeitems));

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Earth: large; Mars: medium; Moon: small.");

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
    return confirmed;
}

bool drawtectonicchooserwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, int& landmass, int& mergefactor)
{
    if (!show)
        return false;

    bool confirmed = false;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 428, main_viewport->WorkPos.y + 167), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(405, 200), ImGuiCond_FirstUseEver);
    ImGui::Begin("Generate world terrain", NULL, window_flags);

    ImGui::Text("This will overwrite any existing terrain.");
    ImGui::Text(" ");

    ImGui::SliderInt("Continental mass", &landmass, 0, 10);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("The (very approximate) amount of land mass, compared to sea.");

    ImGui::Text(" ");

    ImGui::SliderInt("Marine flooding", &mergefactor, 1, 30);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("The degree to which sea covers the continents. The higher this is, the more inland seas and fragmented coastlines there will be.");

    ImGui::Text(" ");
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

bool drawnontectonicchooserwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, int& sealeveleditable, int& iterations)
{
    if (!show)
        return false;

    bool confirmed = false;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 428, main_viewport->WorkPos.y + 167), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(405, 200), ImGuiCond_FirstUseEver);
    ImGui::Begin("Generate world terrain##nontectonic", NULL, window_flags);

    ImGui::Text("This will overwrite any existing terrain.");
    ImGui::Text(" ");

    ImGui::SliderInt("Sea level", &sealeveleditable, 0, 10);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("The approximate sea level. The higher it is, the more sea (very roughly) there is likely to be.");

    ImGui::Text(" ");

    ImGui::SliderInt("Variation", &iterations, 1, 10);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("The amount of variation in terrain type. The higher this is, the more chaotic the terrain may become.");

    ImGui::Text(" ");
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

    ImGui::SameLine((float)rightalign);

    if (ImGui::Button("Close", ImVec2(120.0f, 0.0f)))
        show = false;

    ImGui::End();

    clampworldpropertycontrols(controls);
}
