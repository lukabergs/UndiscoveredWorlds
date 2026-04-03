#include <algorithm>
#include <array>
#include <cmath>
#include <string>

#include "imgui.h"

#include "map_appearance.hpp"

using namespace std;

namespace
{
    struct GradientPreviewStop
    {
        float position;
        const ImVec4* colour;
    };

    struct mapgradientrange
    {
        int minimum;
        int maximum;
        const char* unittext;
    };

    float clampgradientfactor(float value)
    {
        if (value < 0.0f)
            return 0.0f;
        if (value > 1.0f)
            return 1.0f;
        return value;
    }

    float gradientpositionfactor(int value, int minimum, int maximum)
    {
        if (maximum <= minimum)
            return 0.0f;

        return clampgradientfactor(static_cast<float>(value - minimum) / static_cast<float>(maximum - minimum));
    }

    ImVec4 lerpgradientcolour(const ImVec4& low, const ImVec4& high, float factor)
    {
        factor = clampgradientfactor(factor);

        return ImVec4(
            low.x + (high.x - low.x) * factor,
            low.y + (high.y - low.y) * factor,
            low.z + (high.z - low.z) * factor,
            1.0f);
    }

    void sortgradientstops(MapGradientSettings& gradient)
    {
        gradient.stopcount = std::clamp(gradient.stopcount, 2, MAPGRADIENTMAXSTOPS);

        std::sort(gradient.stops.begin(), gradient.stops.begin() + gradient.stopcount, [](const GradientStopSettings& left, const GradientStopSettings& right)
        {
            return left.position < right.position;
        });
    }

    ImVec4 samplegradientcolour(const MapGradientSettings& gradient, int value)
    {
        if (gradient.stopcount <= 0)
            return ImVec4(0.0f, 0.0f, 0.0f, 1.0f);

        if (value <= gradient.stops[0].position)
            return gradient.stops[0].colour;

        for (int i = 0; i < gradient.stopcount - 1; i++)
        {
            const GradientStopSettings& low = gradient.stops[i];
            const GradientStopSettings& high = gradient.stops[i + 1];

            if (value <= high.position)
            {
                if (gradient.discrete || high.position <= low.position)
                    return low.colour;

                const float factor = static_cast<float>(value - low.position) / static_cast<float>(high.position - low.position);
                return lerpgradientcolour(low.colour, high.colour, factor);
            }
        }

        return gradient.stops[gradient.stopcount - 1].colour;
    }

    bool drawgradienteditor(const char* id, MapGradientSettings& gradient, int minimum, int maximum, const char* unittext, int& selectedstop)
    {
        bool changed = false;

        if (maximum <= minimum)
            maximum = minimum + 1;

        sortgradientstops(gradient);
        selectedstop = std::clamp(selectedstop, 0, gradient.stopcount - 1);

        const float width = ImGui::GetContentRegionAvail().x;
        const float height = 26.0f;
        const float markerheight = 12.0f;
        const float markersize = 7.0f;
        const ImVec2 start = ImGui::GetCursorScreenPos();
        const ImVec2 size(width, height + markerheight + 6.0f);
        ImDrawList* drawlist = ImGui::GetWindowDrawList();

        ImGui::InvisibleButton(id, size);

        const bool hovered = ImGui::IsItemHovered();
        const ImVec2 mouse = ImGui::GetIO().MousePos;

        auto valuefrommouse = [&](float mousex)
        {
            const float factor = clampgradientfactor((mousex - start.x) / std::max(1.0f, width));
            return minimum + static_cast<int>(roundf(factor * static_cast<float>(maximum - minimum)));
        };

        auto stopx = [&](int index)
        {
            return start.x + width * gradientpositionfactor(gradient.stops[index].position, minimum, maximum);
        };

        int clickedstop = -1;
        float clickeddistance = 999999.0f;

        if (hovered && ImGui::IsMouseClicked(0))
        {
            for (int i = 0; i < gradient.stopcount; i++)
            {
                const float distance = fabsf(stopx(i) - mouse.x);

                if (distance < clickeddistance && distance <= 10.0f)
                {
                    clickeddistance = distance;
                    clickedstop = i;
                }
            }

            if (clickedstop >= 0)
            {
                selectedstop = clickedstop;
            }
            else if (gradient.stopcount < MAPGRADIENTMAXSTOPS)
            {
                const int newposition = valuefrommouse(mouse.x);
                const ImVec4 newcolour = samplegradientcolour(gradient, newposition);
                gradient.stops[gradient.stopcount].position = newposition;
                gradient.stops[gradient.stopcount].colour = newcolour;
                gradient.stopcount++;
                sortgradientstops(gradient);

                for (int i = 0; i < gradient.stopcount; i++)
                {
                    if (gradient.stops[i].position == newposition)
                    {
                        selectedstop = i;
                        break;
                    }
                }

                changed = true;
            }
        }

        if (ImGui::IsItemActive() && hovered && ImGui::IsMouseDragging(ImGuiMouseButton_Left) && selectedstop >= 0 && selectedstop < gradient.stopcount)
        {
            int minposition = minimum;
            int maxposition = maximum;

            if (selectedstop > 0)
                minposition = gradient.stops[selectedstop - 1].position;
            if (selectedstop < gradient.stopcount - 1)
                maxposition = gradient.stops[selectedstop + 1].position;

            const int draggedposition = std::clamp(valuefrommouse(mouse.x), minposition, maxposition);

            if (gradient.stops[selectedstop].position != draggedposition)
            {
                gradient.stops[selectedstop].position = draggedposition;
                changed = true;
            }
        }

        for (int i = 0; i < gradient.stopcount - 1; i++)
        {
            const float startx = stopx(i);
            const float endx = stopx(i + 1);
            const ImU32 startcolour = ImGui::ColorConvertFloat4ToU32(gradient.stops[i].colour);
            const ImU32 endcolour = ImGui::ColorConvertFloat4ToU32(gradient.discrete ? gradient.stops[i].colour : gradient.stops[i + 1].colour);

            drawlist->AddRectFilledMultiColor(
                ImVec2(startx, start.y),
                ImVec2(endx, start.y + height),
                startcolour,
                endcolour,
                endcolour,
                startcolour);
        }

        drawlist->AddRect(ImVec2(start.x, start.y), ImVec2(start.x + width, start.y + height), IM_COL32(180, 180, 180, 255));

        for (int i = 0; i < gradient.stopcount; i++)
        {
            const float markerx = stopx(i);
            const ImU32 markercolour = ImGui::ColorConvertFloat4ToU32(gradient.stops[i].colour);
            const ImU32 outlinecolour = i == selectedstop ? IM_COL32(255, 255, 0, 255) : IM_COL32(255, 255, 255, 235);

            drawlist->AddLine(ImVec2(markerx, start.y), ImVec2(markerx, start.y + height), outlinecolour);
            drawlist->AddTriangleFilled(
                ImVec2(markerx, start.y + height + 2.0f),
                ImVec2(markerx - markersize, start.y + height + markerheight),
                ImVec2(markerx + markersize, start.y + height + markerheight),
                markercolour);
            drawlist->AddTriangle(
                ImVec2(markerx, start.y + height + 2.0f),
                ImVec2(markerx - markersize, start.y + height + markerheight),
                ImVec2(markerx + markersize, start.y + height + markerheight),
                outlinecolour,
                2.0f);
        }

        ImGui::Text("%d%s", minimum, unittext);
        ImGui::SameLine();
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + std::max(0.0f, width - 150.0f));
        ImGui::Text("%d%s", maximum, unittext);

        if (selectedstop >= 0 && selectedstop < gradient.stopcount)
        {
            ImGui::Checkbox(("Discrete##" + std::string(id)).c_str(), &gradient.discrete);
            changed = changed || ImGui::IsItemDeactivatedAfterEdit();

            ImGui::SameLine(220.0f);
            ImGui::Text("Click gradient to add anchors.");

            int position = gradient.stops[selectedstop].position;
            if (ImGui::InputInt(("Value##" + std::string(id)).c_str(), &position))
            {
                int minposition = minimum;
                int maxposition = maximum;

                if (selectedstop > 0)
                    minposition = gradient.stops[selectedstop - 1].position;
                if (selectedstop < gradient.stopcount - 1)
                    maxposition = gradient.stops[selectedstop + 1].position;

                gradient.stops[selectedstop].position = std::clamp(position, minposition, maxposition);
                changed = true;
            }

            ImGui::SameLine(220.0f);
            changed = ImGui::ColorEdit3(("Colour##" + std::string(id)).c_str(), (float*)&gradient.stops[selectedstop].colour) || changed;

            const bool candelete = gradient.stopcount > 2;
            if (!candelete)
                ImGui::BeginDisabled();

            if (ImGui::Button(("Delete selected##" + std::string(id)).c_str()))
            {
                for (int i = selectedstop; i < gradient.stopcount - 1; i++)
                    gradient.stops[i] = gradient.stops[i + 1];

                gradient.stopcount--;
                selectedstop = std::clamp(selectedstop, 0, gradient.stopcount - 1);
                changed = true;
            }

            if (!candelete)
                ImGui::EndDisabled();
        }

        return changed;
    }

    void drawreliefindexedtab(AppearanceSettings& appearance, int colouralign, int otheralign)
    {
        ImGui::BeginChild("ReliefAppearanceTab", ImVec2(0.0f, 500.0f), false);
        ImGui::PushItemWidth(200);

        ImGui::ColorEdit3("Shallow ocean", (float*)&appearance.oceancolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Sea ice", (float*)&appearance.seaicecolour);

        ImGui::ColorEdit3("Deep ocean", (float*)&appearance.deepoceancolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Glaciers", (float*)&appearance.glaciercolour);

        ImGui::ColorEdit3("Base land", (float*)&appearance.basecolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Salt pans", (float*)&appearance.saltpancolour);

        ImGui::ColorEdit3("Grassland", (float*)&appearance.grasscolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Dunes", (float*)&appearance.ergcolour);

        ImGui::ColorEdit3("Low temperate", (float*)&appearance.basetempcolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Rivers", (float*)&appearance.rivercolour);

        ImGui::ColorEdit3("High temperate", (float*)&appearance.highbasecolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Lakes", (float*)&appearance.lakecolour);

        ImGui::ColorEdit3("Low desert", (float*)&appearance.desertcolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Wetlands", (float*)&appearance.wetlandscolour);

        ImGui::ColorEdit3("High desert", (float*)&appearance.highdesertcolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Mangroves", (float*)&appearance.mangrovecolour);

        ImGui::ColorEdit3("Cold desert", (float*)&appearance.colddesertcolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Mud", (float*)&appearance.mudcolour);

        ImGui::ColorEdit3("Mild tundra", (float*)&appearance.eqtundracolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Sand", (float*)&appearance.sandcolour);

        ImGui::ColorEdit3("Tundra", (float*)&appearance.tundracolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Shingle", (float*)&appearance.shinglecolour);

        ImGui::ColorEdit3("Arctic", (float*)&appearance.coldcolour);
        ImGui::SameLine((float)colouralign);
        ImGui::ColorEdit3("Highlights", (float*)&appearance.highlightcolour);

        ImGui::Separator();
        ImGui::TextUnformatted("Effects");
        ImGui::Text(" ");

        ImGui::Text("Shading");
        ImGui::SameLine((float)otheralign);
        ImGui::Text("Rivers");

        ImGui::SetCursorPosX(20);
        ImGui::SliderFloat("On land", &appearance.shadingland, 0.0f, 1.0f, "%.2f");

        ImGui::SameLine((float)otheralign + 20);
        ImGui::InputInt("Global map", &appearance.globalriversentry);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Only rivers greater than this size will be displayed on the global relief map.");

        ImGui::SetCursorPosX(20);
        ImGui::SliderFloat("On lakes", &appearance.shadinglake, 0.0f, 1.0f, "%.2f");

        ImGui::SameLine((float)otheralign + 20);
        ImGui::InputInt("Regional map", &appearance.regionalriversentry);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Only rivers greater than this size will be displayed on the regional relief map.");

        ImGui::SetCursorPosX(20);
        ImGui::SliderFloat("On sea", &appearance.shadingsea, 0.0f, 1.0f, "%.2f");

        ImGui::Text("Marbling");
        ImGui::SameLine((float)otheralign);
        ImGui::Text("Other effects");

        ImGui::SetCursorPosX(20);
        ImGui::SliderFloat("On land ", &appearance.marblingland, 0.0f, 1.0f, "%.2f");

        ImGui::SameLine((float)otheralign + 20);
        const char* lightdiritems[] = { "Southeast","Southwest","Northeast","Northwest" };
        ImGui::Combo("Light", &appearance.shadingdir, lightdiritems, 4);

        ImGui::SetCursorPosX(20);
        ImGui::SliderFloat("On lakes ", &appearance.marblinglake, 0.0f, 1.0f, "%.2f");

        ImGui::SameLine((float)otheralign + 20);
        const char* snowitems[] = { "Sudden","Speckled","Smooth" };
        ImGui::Combo("Snow", &appearance.snowchange, snowitems, 3);

        ImGui::SetCursorPosX(20);
        ImGui::SliderFloat("On sea ", &appearance.marblingsea, 0.0f, 1.0f, "%.2f");

        ImGui::SameLine((float)otheralign + 20);
        const char* seaiceitems[] = { "Permanent","None","All" };
        ImGui::Combo("Sea ice", &appearance.seaiceappearance, seaiceitems, 3);

        ImGui::SetCursorPosX(20);
        ImGui::Checkbox("Show mangrove forests", &appearance.mangroves);

        ImGui::SameLine((float)otheralign + 20);
        ImGui::Checkbox("Only cliffs use high base colour", &appearance.colourcliffs);

        ImGui::EndChild();
    }

    void drawclimateindexedtab(planet& world, AppearanceSettings& appearance)
    {
        ImGui::BeginChild("ClimateAppearanceTab", ImVec2(0.0f, 500.0f), false);
        ImGui::PushItemWidth(200);

        ImGui::TextUnformatted("Sea");
        ImGui::ColorEdit3("Open sea", (float*)&appearance.climateseacolours[climateopensea]);
        ImGui::SameLine(390.0f);
        ImGui::ColorEdit3("Seasonal sea ice", (float*)&appearance.climateseacolours[climateseasonalseaice]);
        ImGui::ColorEdit3("Permanent sea ice", (float*)&appearance.climateseacolours[climatepermanentseaice]);

        ImGui::Text(" ");
        ImGui::TextUnformatted("Land climates");

        if (ImGui::BeginTable("ClimateAppearanceTable", 3, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_ScrollY, ImVec2(0.0f, 360.0f)))
        {
            ImGui::TableSetupColumn("Code", ImGuiTableColumnFlags_WidthFixed, 55.0f);
            ImGui::TableSetupColumn("Climate");
            ImGui::TableSetupColumn("Colour", ImGuiTableColumnFlags_WidthFixed, 220.0f);
            ImGui::TableHeadersRow();

            for (int climateindex = 1; climateindex <= world.climatenumber() && climateindex < CLIMATEMAPCOLOURCOUNT; climateindex++)
            {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::Text("%s", getclimatecode((short)climateindex).c_str());
                ImGui::TableSetColumnIndex(1);
                ImGui::Text("%s", getclimatename((short)climateindex).c_str());
                ImGui::TableSetColumnIndex(2);
                ImGui::ColorEdit3(("##ClimateColour" + to_string(climateindex)).c_str(), (float*)&appearance.climatecolours[climateindex]);
            }

            ImGui::EndTable();
        }

        ImGui::EndChild();
    }

    void drawbiomeindexedtab(AppearanceSettings& appearance)
    {
        ImGui::BeginChild("BiomeAppearanceTab", ImVec2(0.0f, 500.0f), false);
        ImGui::PushItemWidth(220);

        ImGui::TextUnformatted("Sea uses the Climate tab sea palette.");
        ImGui::Text(" ");

        if (ImGui::BeginTable("BiomeAppearanceTable", 2, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_ScrollY, ImVec2(0.0f, 420.0f)))
        {
            ImGui::TableSetupColumn("Biome");
            ImGui::TableSetupColumn("Colour", ImGuiTableColumnFlags_WidthFixed, 220.0f);
            ImGui::TableHeadersRow();

            for (int biomeindex = 0; biomeindex < BIOMEMAPCOLOURCOUNT; biomeindex++)
            {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::Text("%s", biomemapnames[biomeindex]);
                ImGui::TableSetColumnIndex(1);
                ImGui::ColorEdit3(("##BiomeColour" + to_string(biomeindex)).c_str(), (float*)&appearance.biomecolours[biomeindex]);
            }

            ImGui::EndTable();
        }

        ImGui::EndChild();
    }

    void drawgradienttab(const char* childid, const char* editorid, const char* description, AppearanceSettings& appearance, int gradientindex, int minimum, int maximum, const char* unittext, int& selectedstop)
    {
        ImGui::BeginChild(childid, ImVec2(0.0f, 500.0f), false);
        ImGui::PushItemWidth(220);

        if (description != nullptr && description[0] != '\0')
        {
            ImGui::TextUnformatted(description);
            ImGui::Text(" ");
        }

        drawgradienteditor(editorid, appearance.mapgradients[gradientindex], minimum, maximum, unittext, selectedstop);

        ImGui::EndChild();
    }

    void drawriversgradienttab(planet& world, AppearanceSettings& appearance, int& selectedstop)
    {
        static const std::array<const char*, RIVERMAPFEATURECOUNT> featurenames = { "Lakes", "Salt pans", "Wetlands", "Mud", "Sand", "Shingle", "Volcanoes" };
        static const std::array<int, RIVERMAPFEATURECOUNT> featurecolours = { rivermaplake, rivermapsaltpan, rivermapwetlands, rivermapmud, rivermapsand, rivermapshingle, rivermapvolcano };

        ImGui::BeginChild("RiversAppearanceTab", ImVec2(0.0f, 500.0f), false);
        ImGui::PushItemWidth(220);
        ImGui::ColorEdit3("Background", (float*)&appearance.rivermapcolours[rivermapbackground]);
        ImGui::Text(" ");
        drawgradienteditor(
            "RiverGradientEditor",
            appearance.mapgradients[mapgradientriverflow],
            0,
            max(1000, world.maxriverflow()),
            "",
            selectedstop);
        ImGui::Text(" ");

        if (ImGui::BeginTable("RiverFeaturesTable", 3, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV))
        {
            ImGui::TableSetupColumn("Feature", ImGuiTableColumnFlags_WidthFixed, 120.0f);
            ImGui::TableSetupColumn("Show", ImGuiTableColumnFlags_WidthFixed, 60.0f);
            ImGui::TableSetupColumn("Colour");
            ImGui::TableHeadersRow();

            for (int featureindex = 0; featureindex < RIVERMAPFEATURECOUNT; featureindex++)
            {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::Text("%s", featurenames[featureindex]);
                ImGui::TableSetColumnIndex(1);
                ImGui::Checkbox(("##RiverFeatureToggle" + to_string(featureindex)).c_str(), &appearance.showrivermapfeatures[featureindex]);
                ImGui::TableSetColumnIndex(2);
                ImGui::ColorEdit3(("##RiverFeatureColour" + to_string(featureindex)).c_str(), (float*)&appearance.rivermapcolours[featurecolours[featureindex]]);
            }

            ImGui::EndTable();
        }

        ImGui::EndChild();
    }

    mapgradientrange getmapgradientrange(const mapviewdefinition& definition, const planet& world)
    {
        switch (definition.gradientindex)
        {
        case mapgradientelevation:
            return { 0 - world.sealevel(), max(1, world.maxelevation() - world.sealevel()), " m" };
        case mapgradienttemperature:
            return { -60, 60, " C" };
        case mapgradientprecipitation:
            return { 0, 1020, "" };
        case mapgradientriverflow:
            return { 0, max(1000, world.maxriverflow()), "" };
        default:
            return { 0, 1, "" };
        }
    }
}

void drawmapviewbuttons(const std::function<void(mapviewenum)>& onSelect)
{
    for (const mapviewdefinition& definition : allmapviewdefinitions)
    {
        if (standardbutton(definition.label))
            onSelect(definition.view);
    }
}

void drawgradientmapappearance(const mapviewdefinition& definition, planet& world, AppearanceSettings& appearance, array<int, MAPGRADIENTTYPECOUNT>& selectedgradientstops, int, int)
{
    switch (definition.gradientstyle)
    {
    case standardgradientstyle:
    {
        const mapgradientrange range = getmapgradientrange(definition, world);

        drawgradienttab(
            (string(definition.label) + "AppearanceTab").c_str(),
            (string(definition.label) + "GradientEditor").c_str(),
            definition.description,
            appearance,
            definition.gradientindex,
            range.minimum,
            range.maximum,
            range.unittext,
            selectedgradientstops[definition.gradientindex]);
        break;
    }
    case riversgradientstyle:
        drawriversgradienttab(world, appearance, selectedgradientstops[definition.gradientindex]);
        break;
    default:
        break;
    }
}

void drawindexedmapappearance(const mapviewdefinition& definition, planet& world, AppearanceSettings& appearance, array<int, MAPGRADIENTTYPECOUNT>&, int colouralign, int otheralign)
{
    switch (definition.indexedstyle)
    {
    case reliefindexedstyle:
        drawreliefindexedtab(appearance, colouralign, otheralign);
        break;
    case climateindexedstyle:
        drawclimateindexedtab(world, appearance);
        break;
    case biomeindexedstyle:
        drawbiomeindexedtab(appearance);
        break;
    default:
        break;
    }
}

void drawstaticmapappearance(const mapviewdefinition& definition, planet&, AppearanceSettings&, array<int, MAPGRADIENTTYPECOUNT>&, int, int)
{
    ImGui::BeginChild((string(definition.label) + "StaticAppearanceTab").c_str(), ImVec2(0.0f, 500.0f), false);

    if (definition.description != nullptr && definition.description[0] != '\0')
        ImGui::TextWrapped("%s", definition.description);
    else
        ImGui::TextUnformatted("This diagnostic map uses a fixed palette.");

    ImGui::Spacing();
    ImGui::TextUnformatted("Appearance settings are fixed for this view.");

    ImGui::EndChild();
}

void drawmapviewappearancetab(const mapviewdefinition& definition, planet& world, AppearanceSettings& appearance, std::array<int, MAPGRADIENTTYPECOUNT>& selectedgradientstops, int colouralign, int otheralign)
{
    definition.drawappearance(definition, world, appearance, selectedgradientstops, colouralign, otheralign);
}
