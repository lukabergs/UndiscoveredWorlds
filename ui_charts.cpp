#include <array>
#include <cmath>
#include <algorithm>
#include <string>

#include "ui_charts.hpp"
#include "functions.hpp"
#include "map_appearance.hpp"

using namespace std;

namespace
{
    constexpr float chartbarwidth = 40.0f;
    constexpr float widgetpadding = 8.0f;
    constexpr float handlegrabradius = 7.0f;
    constexpr float curvesegments = 48.0f;
    const array<const char*, 12> monthlabels = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

    ImVec2 operator+(const ImVec2& lhs, const ImVec2& rhs)
    {
        return ImVec2(lhs.x + rhs.x, lhs.y + rhs.y);
    }

    ImVec2 operator-(const ImVec2& lhs, const ImVec2& rhs)
    {
        return ImVec2(lhs.x - rhs.x, lhs.y - rhs.y);
    }

    ImVec2 operator*(const ImVec2& lhs, float rhs)
    {
        return ImVec2(lhs.x * rhs, lhs.y * rhs);
    }

    ImVec2 clamp01(const ImVec2& value)
    {
        return ImVec2(clamp(value.x, 0.0f, 1.0f), clamp(value.y, 0.0f, 1.0f));
    }

    float normalizedvalue(float value, float minimumvalue, float maximumvalue)
    {
        if (maximumvalue <= minimumvalue)
            return 0.0f;

        return clamp((value - minimumvalue) / (maximumvalue - minimumvalue), 0.0f, 1.0f);
    }

    ImVec2 toscreen(const ImVec2& origin, const ImVec2& size, const ImVec2& normalized)
    {
        return ImVec2(origin.x + normalized.x * size.x, origin.y + (1.0f - normalized.y) * size.y);
    }

    ImVec2 tonormalized(const ImVec2& origin, const ImVec2& size, const ImVec2& screen)
    {
        const float width = max(size.x, 1.0f);
        const float height = max(size.y, 1.0f);
        return ImVec2((screen.x - origin.x) / width, 1.0f - ((screen.y - origin.y) / height));
    }

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

    void drawchartframe(const char* label, const ImVec2& size)
    {
        ImGui::TextUnformatted(label);
        ImVec2 canvas_size = size;
        if (canvas_size.x <= 0.0f)
            canvas_size.x = max(240.0f, ImGui::GetContentRegionAvail().x);
        if (canvas_size.y <= 0.0f)
            canvas_size.y = 220.0f;
        ImGui::PushID(label);
        ImGui::InvisibleButton("##chartframe", canvas_size);
        ImGui::PopID();
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

    void drawcurvehandle(ImDrawList* drawlist, const ImVec2& center, ImU32 color)
    {
        drawlist->AddCircleFilled(center, handlegrabradius, color, 16);
        drawlist->AddCircle(center, handlegrabradius, IM_COL32(0, 0, 0, 180), 16, 2.0f);
    }

    float clampcurveweight(float value)
    {
        return clamp(value, 0.1f, 8.0f);
    }

    int curvedegree(size_t controlpointcount)
    {
        if (controlpointcount <= 1)
            return 0;

        return min(3, static_cast<int>(controlpointcount) - 1);
    }

    vector<float> buildcurveknots(size_t controlpointcount)
    {
        const int degree = curvedegree(controlpointcount);
        const size_t knotcount = controlpointcount + static_cast<size_t>(degree) + 1;
        vector<float> knots(knotcount, 0.0f);

        if (controlpointcount == 0)
            return knots;

        const size_t interiorcount = knotcount - static_cast<size_t>((degree + 1) * 2);

        for (size_t index = 0; index < knotcount; index++)
        {
            if (index <= static_cast<size_t>(degree))
                knots[index] = 0.0f;
            else if (index >= controlpointcount)
                knots[index] = 1.0f;
            else
                knots[index] = static_cast<float>(index - degree) / static_cast<float>(interiorcount + 1);
        }

        return knots;
    }

    float evaluatecurvebasis(int index, int degree, float t, const vector<float>& knots)
    {
        if (degree == 0)
            return (t >= knots[index] && t < knots[index + 1]) || (t >= 1.0f && knots[index + 1] >= 1.0f) ? 1.0f : 0.0f;

        float left = 0.0f;
        float right = 0.0f;
        const float leftdenominator = knots[index + degree] - knots[index];
        const float rightdenominator = knots[index + degree + 1] - knots[index + 1];

        if (leftdenominator > 0.0001f)
            left = ((t - knots[index]) / leftdenominator) * evaluatecurvebasis(index, degree - 1, t, knots);

        if (rightdenominator > 0.0001f)
            right = ((knots[index + degree + 1] - t) / rightdenominator) * evaluatecurvebasis(index + 1, degree - 1, t, knots);

        return left + right;
    }

    ImVec2 evaluatenurbscurve(const vector<ImVec2>& controlpoints, const vector<float>& weights, float t)
    {
        if (controlpoints.empty())
            return ImVec2(0.0f, 0.0f);

        const int degree = curvedegree(controlpoints.size());
        const vector<float> knots = buildcurveknots(controlpoints.size());
        const float clampedt = clamp(t, 0.0f, 1.0f);
        float totalweight = 0.0f;
        ImVec2 point(0.0f, 0.0f);

        for (size_t index = 0; index < controlpoints.size(); index++)
        {
            const float basis = evaluatecurvebasis(static_cast<int>(index), degree, clampedt, knots);
            const float weight = clampcurveweight(index < weights.size() ? weights[index] : 1.0f);
            const float scaledbasis = basis * weight;
            totalweight += scaledbasis;
            point.x += controlpoints[index].x * scaledbasis;
            point.y += controlpoints[index].y * scaledbasis;
        }

        if (totalweight <= 0.0001f)
            return clamp01(controlpoints.back());

        return ImVec2(point.x / totalweight, point.y / totalweight);
    }

    vector<ImVec2> sampleweightedcurve(const vector<ImVec2>& controlpoints, const vector<float>& weights)
    {
        vector<ImVec2> points;
        points.reserve(static_cast<size_t>(curvesegments) + 1);

        for (int index = 0; index <= static_cast<int>(curvesegments); index++)
            points.push_back(evaluatenurbscurve(controlpoints, weights, static_cast<float>(index) / curvesegments));

        return points;
    }

    string formatmeasurement(float value, const char* unittext)
    {
        return formatnumber(static_cast<int>(roundf(value))) + (unittext != nullptr ? unittext : "");
    }

    void drawgradientbackground(ImDrawList* drawlist, const ImVec2& minpos, const ImVec2& maxpos, float minimumvalue, float maximumvalue, const MapGradientSettings* gradient)
    {
        const int steps = max(1, static_cast<int>(roundf(maxpos.x - minpos.x)));

        for (int step = 0; step < steps; step++)
        {
            const float factor0 = static_cast<float>(step) / static_cast<float>(steps);
            const float factor1 = static_cast<float>(step + 1) / static_cast<float>(steps);
            const float value = minimumvalue + ((factor0 + factor1) * 0.5f) * (maximumvalue - minimumvalue);
            const ImVec4 colour = gradient != nullptr
                ? samplegradientcolour(*gradient, static_cast<int>(roundf(value)))
                : ImVec4(0.16f, 0.18f, 0.22f, 1.0f);
            const ImU32 fill = ImGui::ColorConvertFloat4ToU32(ImVec4(colour.x, colour.y, colour.z, 0.28f));

            drawlist->AddRectFilled(
                ImVec2(minpos.x + factor0 * (maxpos.x - minpos.x), minpos.y),
                ImVec2(minpos.x + factor1 * (maxpos.x - minpos.x), maxpos.y),
                fill);
        }
    }

    void drawcontinuoushistogram(ImDrawList* drawlist, const ImVec2& minpos, const ImVec2& maxpos, const vector<float>& values, float minimumvalue, float maximumvalue, float sealevelvalue, const MapGradientSettings* gradient)
    {
        const int count = static_cast<int>(values.size());
        if (count <= 0)
            return;

        const float width = maxpos.x - minpos.x;
        const float height = maxpos.y - minpos.y;
        float maxcount = 0.0f;

        for (float value : values)
            maxcount = max(maxcount, value);

        if (maxcount <= 0.0f)
            maxcount = 1.0f;

        drawgradientbackground(drawlist, minpos, maxpos, minimumvalue, maximumvalue, gradient);

        vector<ImVec2> curvepoints;
        curvepoints.reserve(static_cast<size_t>(count));

        for (int index = 0; index < count; index++)
        {
            const float x = minpos.x + (static_cast<float>(index) / static_cast<float>(max(1, count - 1))) * width;
            const float y = maxpos.y - clamp(values[index] / maxcount, 0.0f, 1.0f) * height;
            curvepoints.emplace_back(x, y);
            drawlist->AddLine(ImVec2(x, maxpos.y), ImVec2(x, y), IM_COL32(255, 255, 255, 24), 1.0f);
        }

        for (int index = 1; index < count; index++)
        {
            const float samplevalue = minimumvalue + (static_cast<float>(index) / static_cast<float>(max(1, count - 1))) * (maximumvalue - minimumvalue);
            const ImVec4 colour = gradient != nullptr
                ? samplegradientcolour(*gradient, static_cast<int>(roundf(samplevalue)))
                : ImVec4(0.72f, 0.78f, 0.88f, 1.0f);
            drawlist->AddLine(curvepoints[index - 1], curvepoints[index], ImGui::ColorConvertFloat4ToU32(ImVec4(colour.x, colour.y, colour.z, 0.92f)), 2.0f);
        }

        const float sealevelx = minpos.x + normalizedvalue(sealevelvalue, minimumvalue, maximumvalue) * width;
        drawlist->AddLine(ImVec2(sealevelx, minpos.y), ImVec2(sealevelx, maxpos.y), IM_COL32(100, 160, 255, 210), 2.0f);
        drawcurvehandle(drawlist, ImVec2(sealevelx, maxpos.y - 10.0f), IM_COL32(100, 160, 255, 255));
        drawlist->AddRect(minpos, maxpos, IM_COL32(255, 255, 255, 60), 4.0f);
    }
}

bool drawheighthistogramwidget(const char* label, const vector<float>& values, float minimumvalue, float maximumvalue, float& sealevelvalue, const MapGradientSettings* gradient, const char* unittext, const ImVec2& size)
{
    bool changed = false;
    ImGui::BeginGroup();
    drawchartframe(label, size);

    const ImVec2 canvasmin = ImGui::GetItemRectMin();
    const ImVec2 canvasmax = ImGui::GetItemRectMax();
    const ImVec2 origin = canvasmin + ImVec2(widgetpadding, widgetpadding);
    const ImVec2 extent = canvasmax - canvasmin - ImVec2(widgetpadding * 2.0f, widgetpadding * 2.0f);
    const ImVec2 mouse = ImGui::GetIO().MousePos;
    ImDrawList* drawlist = ImGui::GetWindowDrawList();
    ImGuiStorage* storage = ImGui::GetStateStorage();
    const ImGuiID widgetid = ImGui::GetID(label);
    int dragmode = storage->GetInt(widgetid, -1);

    drawlist->AddRectFilled(canvasmin, canvasmax, IM_COL32(24, 28, 34, 255), 4.0f);
    drawlist->AddRect(canvasmin, canvasmax, IM_COL32(255, 255, 255, 35), 4.0f);
    drawcontinuoushistogram(drawlist, origin, origin + extent, values, minimumvalue, maximumvalue, sealevelvalue, gradient);

    if (ImGui::IsItemHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
    {
        const float sealevelx = origin.x + normalizedvalue(sealevelvalue, minimumvalue, maximumvalue) * extent.x;
        dragmode = fabsf(mouse.x - sealevelx) <= handlegrabradius * 2.0f ? 0 : -1;
        storage->SetInt(widgetid, dragmode);
    }

    if (dragmode == 0 && ImGui::IsMouseDown(ImGuiMouseButton_Left))
    {
        const float factor = clamp((mouse.x - origin.x) / max(extent.x, 1.0f), 0.0f, 1.0f);
        sealevelvalue = minimumvalue + factor * (maximumvalue - minimumvalue);
        changed = true;
    }

    if (dragmode != -1 && ImGui::IsMouseReleased(ImGuiMouseButton_Left))
        storage->SetInt(widgetid, -1);

    const string mintext = formatmeasurement(minimumvalue, unittext);
    const string maxtext = formatmeasurement(maximumvalue, unittext);
    const float baseline = ImGui::GetCursorPosX();
    const float availwidth = ImGui::GetContentRegionAvail().x;

    ImGui::TextUnformatted(mintext.c_str());
    ImGui::SameLine();
    ImGui::SetCursorPosX(baseline + max(0.0f, availwidth - ImGui::CalcTextSize(maxtext.c_str()).x));
    ImGui::TextUnformatted(maxtext.c_str());

    ImGui::EndGroup();
    return changed;
}

bool drawhypsometriccurvewidget(const char* label, const vector<ImVec2>& worldcurve, vector<ImVec2>& controlpoints, vector<float>& weights, int& selectedcontrolpoint, float& sealevelvalue, const ImVec2& size, const vector<ImVec2>* previewcurve, bool editable)
{
    bool changed = false;

    if (controlpoints.size() < 2)
        controlpoints = { ImVec2(0.0f, 0.0f), ImVec2(1.0f, 1.0f) };

    if (weights.size() != controlpoints.size())
        weights.assign(controlpoints.size(), 1.0f);

    controlpoints.front() = ImVec2(0.0f, 0.0f);
    controlpoints.back() = ImVec2(1.0f, 1.0f);
    weights.front() = 1.0f;
    weights.back() = 1.0f;
    selectedcontrolpoint = controlpoints.size() > 2
        ? clamp(selectedcontrolpoint, 1, static_cast<int>(controlpoints.size()) - 2)
        : -1;

    ImGui::BeginGroup();
    drawchartframe(label, size);

    const ImVec2 canvasmin = ImGui::GetItemRectMin();
    const ImVec2 canvasmax = ImGui::GetItemRectMax();
    const ImVec2 origin = canvasmin + ImVec2(widgetpadding, widgetpadding);
    const ImVec2 extent = canvasmax - canvasmin - ImVec2(widgetpadding * 2.0f, widgetpadding * 2.0f);
    const ImVec2 mouse = ImGui::GetIO().MousePos;
    ImDrawList* drawlist = ImGui::GetWindowDrawList();
    ImGuiStorage* storage = ImGui::GetStateStorage();
    const ImGuiID widgetid = ImGui::GetID(label);
    int dragmode = storage->GetInt(widgetid, -1);

    drawlist->AddRectFilled(canvasmin, canvasmax, IM_COL32(24, 28, 34, 255), 4.0f);
    drawlist->AddRect(canvasmin, canvasmax, IM_COL32(255, 255, 255, 35), 4.0f);

    const vector<ImVec2> curve = previewcurve != nullptr ? *previewcurve : sampleweightedcurve(controlpoints, weights);

    drawlist->AddLine(ImVec2(origin.x, origin.y + extent.y), ImVec2(origin.x + extent.x, origin.y + extent.y), IM_COL32(120, 120, 120, 120), 1.0f);
    drawlist->AddLine(ImVec2(origin.x, origin.y), ImVec2(origin.x, origin.y + extent.y), IM_COL32(120, 120, 120, 120), 1.0f);

    if (worldcurve.size() >= 2)
    {
        for (size_t index = 1; index < worldcurve.size(); index++)
        {
            const ImVec2 a = toscreen(origin, extent, clamp01(worldcurve[index - 1]));
            const ImVec2 b = toscreen(origin, extent, clamp01(worldcurve[index]));
            drawlist->AddLine(a, b, IM_COL32(120, 190, 255, 200), 2.0f);
        }
    }

    for (size_t index = 1; index < curve.size(); index++)
        drawlist->AddLine(toscreen(origin, extent, curve[index - 1]), toscreen(origin, extent, curve[index]), IM_COL32(240, 185, 95, 220), 2.0f);

    const float sealevely = origin.y + (1.0f - clamp(sealevelvalue, 0.0f, 1.0f)) * extent.y;
    drawlist->AddLine(ImVec2(origin.x, sealevely), ImVec2(origin.x + extent.x, sealevely), IM_COL32(100, 160, 255, 200), 2.0f);
    drawcurvehandle(drawlist, ImVec2(origin.x + 10.0f, sealevely), IM_COL32(100, 160, 255, 255));

    if (editable)
    {
        for (size_t index = 1; index + 1 < controlpoints.size(); index++)
        {
            const ImU32 colour = static_cast<int>(index) == selectedcontrolpoint ? IM_COL32(255, 218, 95, 255) : IM_COL32(240, 185, 95, 255);
            drawcurvehandle(drawlist, toscreen(origin, extent, clamp01(controlpoints[index])), colour);
        }
    }

    if (editable && ImGui::IsItemHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
    {
        dragmode = -1;
        float bestdistance = 999999.0f;

        for (size_t index = 1; index + 1 < controlpoints.size(); index++)
        {
            const ImVec2 handle = toscreen(origin, extent, clamp01(controlpoints[index]));
            const ImVec2 delta = mouse - handle;
            const float distance = delta.x * delta.x + delta.y * delta.y;

            if (distance < bestdistance && distance <= handlegrabradius * handlegrabradius * 4.0f)
            {
                bestdistance = distance;
                selectedcontrolpoint = static_cast<int>(index);
                dragmode = static_cast<int>(index) + 1;
            }
        }

        if (dragmode == -1 && fabsf(mouse.y - sealevely) <= handlegrabradius * 2.0f)
            dragmode = 0;
        else if (dragmode == -1 && controlpoints.size() < MAPGRADIENTMAXSTOPS && mouse.x >= origin.x && mouse.x <= origin.x + extent.x && mouse.y >= origin.y && mouse.y <= origin.y + extent.y)
        {
            const ImVec2 local = clamp01(tonormalized(origin, extent, mouse));
            size_t insertindex = 1;

            while (insertindex < controlpoints.size() && controlpoints[insertindex].x < local.x)
                insertindex++;

            const float minx = controlpoints[insertindex - 1].x;
            const float maxx = controlpoints[insertindex].x;
            const float miny = controlpoints[insertindex - 1].y;
            const float maxy = controlpoints[insertindex].y;
            controlpoints.insert(controlpoints.begin() + insertindex, ImVec2(clamp(local.x, minx, maxx), clamp(local.y, miny, maxy)));
            weights.insert(weights.begin() + insertindex, 1.0f);
            selectedcontrolpoint = static_cast<int>(insertindex);
            dragmode = selectedcontrolpoint + 1;
            changed = true;
        }

        storage->SetInt(widgetid, dragmode);
    }

    if (editable && dragmode != -1 && ImGui::IsMouseDown(ImGuiMouseButton_Left))
    {
        const ImVec2 local = tonormalized(origin, extent, mouse);

        if (dragmode == 0)
        {
            sealevelvalue = clamp(local.y, 0.0f, 1.0f);
            changed = true;
        }
        else
        {
            const int pointindex = dragmode - 1;

            if (pointindex > 0 && pointindex < static_cast<int>(controlpoints.size()) - 1)
            {
                const float minx = controlpoints[pointindex - 1].x;
                const float maxx = controlpoints[pointindex + 1].x;
                const float miny = controlpoints[pointindex - 1].y;
                const float maxy = controlpoints[pointindex + 1].y;
                controlpoints[pointindex] = ImVec2(clamp(local.x, minx, maxx), clamp(local.y, miny, maxy));
                changed = true;
            }
        }
    }

    if (editable && dragmode != -1 && ImGui::IsMouseReleased(ImGuiMouseButton_Left))
        storage->SetInt(widgetid, -1);

    ImGui::EndGroup();
    return changed;
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
