#include "app/application.hpp"

void Application::updateglobalmapscreen()
{
    ImGuiIO& io = ImGui::GetIO();

    if (screenmode == globalmapscreen)
    {
        showworldeditproperties = 0;
        showregionaltemperaturechart = 0;
        showregionalrainfallchart = 0;

        areafromregional = 0;
        brandnew = 0;

        // Main controls.

        string title;

        if (world->seed() >= 0)
            title = "Seed: " + to_string(world->seed());
       else
            title = "Custom";

        ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 10, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(161, 549), ImGuiCond_FirstUseEver);

        ImGui::Begin(title.c_str(), NULL, window_flags);

        ImVec2 pos = ImGui::GetWindowPos();

        ImGui::Text("World controls:");

        if (standardbutton("New world"))
        {
            brandnew = 0;
            seedentry = 0;

            screenmode = createworldscreen;
        }

        if (standardbutton("Load world"))
        {
            openFileDialog(kWorldLoadFilter);

            loadingworld = 1;
        }

        if (standardbutton("Save world"))
        {
            openFileDialog(".uww");

            savingworld = 1;
        }

        if (standardbutton("Custom world"))
        {
            showsetsize = 1;
        }

        ImGui::Dummy(ImVec2(0.0f, linespace));

        ImGui::SetNextItemWidth(0);

        ImGui::Text("Export options:");

        ImGui::PushItemWidth(100.0f);

        if (standardbutton("World maps"))
        {
            openFileDialog(".png");

            exportingworldmaps = 1;
        }

        if (standardbutton("Area maps"))
        {
            redrawGlobalRelief();

            screenmode = exportareascreen;
            areafromregional = 0;
        }

        ImGui::Dummy(ImVec2(0.0f, linespace));

        ImGui::SetNextItemWidth(0);

        ImGui::Text("Display map type:");

        ImGui::PushItemWidth(100.0f);

        drawmapviewbuttons([&](mapviewenum view) { showGlobalMapView(view); });

        ImGui::Dummy(ImVec2(0.0f, linespace));

        ImGui::SetNextItemWidth(0);

        ImGui::Text("Other controls:");

        ImGui::PushItemWidth(100.0f);

        if (standardbutton("Properties"))
            toggle(showworldproperties);

        if (standardbutton("Appearance"))
            toggle(showcolouroptions);

        if (standardbutton("Zoom"))
        {
            if (focused == 1)
            {
                redrawGlobalRelief(true);
                newx = poix;
                newy = poiy;

                screenmode = generatingregionscreen;
            }
        }

        ImGui::End();

        drawglobalinfowindow(main_viewport, window_flags, *world, focused == 1, poix, poiy, degree, cube, showglobaltemperaturechart, showglobalrainfallchart, newworld == 1, showabout);

        // Now check to see if the map has been clicked on.

        if (window.hasFocus() && sf::Mouse::isButtonPressed(sf::Mouse::Left) && io.WantCaptureMouse == 0)
        {
            sf::Vector2i mousepos = sf::Mouse::getPosition(window);

            float mult = ((float)world->width()+1.0f) / (float)DISPLAYMAPSIZEX;

            float fpoix = (float)(mousepos.x - globalmapxpos) * mult;
            float fpoiy = (float)(mousepos.y - globalmapypos) * mult;

            poix = (int)fpoix;
            poiy = (int)fpoiy;

            if (poix >= 0 && poix <= world->width() && poiy >= 0 && poiy <= world->height())
            {
                focused = 1;
                newworld = 0;

                highlight.setPosition(sf::Vector2f((float)mousepos.x, (float)mousepos.y));
            }
            else
            {
                if (focused == 1)
                {
                    focused = 0;
                    poix = -1;
                    poiy = -1;
                }
            }
        }
    }
}
