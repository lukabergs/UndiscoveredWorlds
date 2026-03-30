#include <array>
#include <cmath>
#include <string>

#include "ui_charts.hpp"
#include "functions.hpp"

using namespace std;

namespace
{
    constexpr float chartbarwidth = 40.0f;
    const array<const char*, 12> monthlabels = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

    array<float, 12> buildquarterlyseries(float jan, float apr, float jul, float oct)
    {
        array<float, 12> values = {};

        values[0] = jan;
        values[3] = apr;
        values[6] = jul;
        values[9] = oct;
        values[1] = (values[0] * 2.0f + values[3]) / 3.0f;
        values[2] = (values[0] + values[3] * 2.0f) / 3.0f;
        values[4] = (values[3] * 2.0f + values[6]) / 3.0f;
        values[5] = (values[3] + values[6] * 2.0f) / 3.0f;
        values[7] = (values[6] * 2.0f + values[9]) / 3.0f;
        values[8] = (values[6] + values[9] * 2.0f) / 3.0f;
        values[10] = (values[9] * 2.0f + values[0]) / 3.0f;
        values[11] = (values[9] + values[0] * 2.0f) / 3.0f;

        return values;
    }

    void drawchartlabels(const array<float, 12>& values)
    {
        ImGui::Text(" ");

        for (int n = 0; n < 12; n++)
        {
            const string valuetext = formatnumber((int)values[n]);

            ImGui::SameLine(chartbarwidth * (n + 1) - ImGui::CalcTextSize(valuetext.c_str()).x / 2.0f);
            ImGui::Text(valuetext.c_str());
        }

        ImGui::Text(" ");

        for (int n = 0; n < 12; n++)
        {
            ImGui::SameLine(chartbarwidth * (n + 1) - ImGui::CalcTextSize(monthlabels[n]).x / 2.0f);
            ImGui::Text(monthlabels[n]);
        }
    }

    void settemperaturechartcolour(float temperature)
    {
        float colour1 = 0.0f;
        float colour2 = 0.0f;
        float colour3 = 0.0f;

        if (temperature > 0.0f)
        {
            colour1 = 250.0f;
            colour2 = 250.0f - temperature * 3.0f;
            colour3 = 250.0f - temperature * 7.0f;
        }
        else
        {
            temperature = abs(temperature);
            colour1 = 250.0f - temperature * 7.0f;
            colour2 = 250.0f - temperature * 7.0f;
            colour3 = 250.0f;
        }

        colour1 = max(colour1, 0.0f) / 255.0f;
        colour2 = max(colour2, 0.0f) / 255.0f;
        colour3 = max(colour3, 0.0f) / 255.0f;

        ImGuiStyle& style = ImGui::GetStyle();
        style.Colors[ImGuiCol_PlotHistogram] = ImVec4(colour1, colour2, colour3, 1.00f);
        style.Colors[ImGuiCol_PlotHistogramHovered] = ImVec4(colour1, colour2, colour3, 1.00f);
    }

    void setrainfallchartcolour(float rainfall)
    {
        float colour1 = 255.0f - rainfall / 2.0f;
        float colour2 = 255.0f - rainfall / 2.0f;
        const float colour3 = 255.0f;

        colour1 = max(colour1, 0.0f) / 255.0f;
        colour2 = max(colour2, 0.0f) / 255.0f;

        ImGuiStyle& style = ImGui::GetStyle();
        style.Colors[ImGuiCol_PlotHistogram] = ImVec4(colour1, colour2, colour3 / 255.0f, 1.00f);
        style.Colors[ImGuiCol_PlotHistogramHovered] = ImVec4(colour1, colour2, colour3 / 255.0f, 1.00f);
    }

    void drawtemperaturechart(const char* title, const ImVec2& position, const array<float, 12>& sourcevalues, ImGuiWindowFlags window_flags)
    {
        array<float, 12> displayvalues = sourcevalues;
        float lowest = displayvalues[0];
        float highest = displayvalues[0];

        for (int n = 1; n < 12; n++)
        {
            lowest = min(lowest, displayvalues[n]);
            highest = max(highest, displayvalues[n]);
        }

        float subzero = 0.0f;

        if (lowest < 0.0f)
        {
            subzero = -lowest;

            for (float& value : displayvalues)
                value += subzero;

            lowest = displayvalues[0];
            highest = displayvalues[0];

            for (int n = 1; n < 12; n++)
            {
                lowest = min(lowest, displayvalues[n]);
                highest = max(highest, displayvalues[n]);
            }
        }

        const float plotmax = highest > 0.0f ? highest : 1.0f;

        ImGui::SetNextWindowPos(position, ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(518, 139), ImGuiCond_FirstUseEver);
        ImGui::Begin(title, NULL, window_flags);

        for (int n = 0; n < 13; n++)
        {
            ImGui::SameLine(chartbarwidth * (float)n + chartbarwidth / 2.0f);

            const float value = n < 12 ? displayvalues[n] : 0.0f;
            const float chartvalue[] = { value / plotmax };

            settemperaturechartcolour((n < 12 ? displayvalues[n] : 0.0f) - subzero);
            ImGui::PlotHistogram(" ", chartvalue, IM_ARRAYSIZE(chartvalue), 0, NULL, 0.0f, 1.0f, ImVec2(0.0f, chartbarwidth));
        }

        array<float, 12> labelvalues = displayvalues;
        for (float& value : labelvalues)
            value -= subzero;

        drawchartlabels(labelvalues);
        ImGui::End();
    }

    void drawrainfallchart(const char* title, const ImVec2& position, const array<float, 12>& values, ImGuiWindowFlags window_flags)
    {
        float highest = values[0];

        for (int n = 1; n < 12; n++)
            highest = max(highest, values[n]);

        const float plotmax = highest > 0.0f ? highest : 1.0f;

        ImGui::SetNextWindowPos(position, ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(518, 139), ImGuiCond_FirstUseEver);
        ImGui::Begin(title, NULL, window_flags);

        for (int n = 0; n < 13; n++)
        {
            ImGui::SameLine(chartbarwidth * (float)n + chartbarwidth / 2.0f);

            const float value = n < 12 ? values[n] : 0.0f;
            const float chartvalue[] = { value / plotmax };

            if (n < 12)
                setrainfallchartcolour(values[n]);
            else
                setrainfallchartcolour(0.0f);

            ImGui::PlotHistogram(" ", chartvalue, IM_ARRAYSIZE(chartvalue), 0, NULL, 0.0f, 1.0f, ImVec2(0.0f, chartbarwidth));
        }

        drawchartlabels(values);
        ImGui::End();
    }
}

void drawglobaltemperaturechartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, int poix, int poiy)
{
    const array<float, 12> values = buildquarterlyseries((float)world.jantemp(poix, poiy), (float)world.aprtemp(poix, poiy), (float)world.jultemp(poix, poiy), (float)world.octtemp(poix, poiy));
    drawtemperaturechart("Temperature", ImVec2(main_viewport->WorkPos.x + 680, main_viewport->WorkPos.y + 385), values, window_flags);
}

void drawglobalrainfallchartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, int poix, int poiy)
{
    const array<float, 12> values = buildquarterlyseries((float)world.janrain(poix, poiy), (float)world.aprrain(poix, poiy), (float)world.julrain(poix, poiy), (float)world.octrain(poix, poiy));
    drawrainfallchart("Precipitation", ImVec2(main_viewport->WorkPos.x + 680, main_viewport->WorkPos.y + 246), values, window_flags);
}

void drawregionaltemperaturechartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, region& region, int poix, int poiy)
{
    const int yy = region.lefty() + poiy / 16;
    const array<float, 12> values = buildquarterlyseries(
        (float)region.jantemp(poix, poiy),
        (float)region.aprtemp(poix, poiy, yy, world.height(), world.tilt(), world.eccentricity(), world.perihelion()),
        (float)region.jultemp(poix, poiy),
        (float)region.aprtemp(poix, poiy, yy, world.height(), world.tilt(), world.eccentricity(), world.perihelion()));

    drawtemperaturechart("Temperature##regional", ImVec2(main_viewport->WorkPos.x + 680, main_viewport->WorkPos.y + 385), values, window_flags);
}

void drawregionalrainfallchartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, region& region, int poix, int poiy)
{
    const array<float, 12> values = buildquarterlyseries((float)region.janrain(poix, poiy), (float)region.aprrain(poix, poiy), (float)region.julrain(poix, poiy), (float)region.octrain(poix, poiy));
    drawrainfallchart("Precipitation##regional", ImVec2(main_viewport->WorkPos.x + 680, main_viewport->WorkPos.y + 246), values, window_flags);
}
