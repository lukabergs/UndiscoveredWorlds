#include "app/application.hpp"

void Application::updategeneratingregionscreen()
{
    if (screenmode == generatingregionscreen)
    {
        showcolouroptions = 0;
        showworldproperties = 0;
        showglobalrainfallchart = 0;
        showglobaltemperaturechart = 0;
        showregionalrainfallchart = 0;
        showregionaltemperaturechart = 0;

        if (generatingregionpass < 10) // This one has a non-functioning copy of the regional map screen controls, plus the "Please wait" message.
        {
            string title;

            if (world->seed() >= 0)
                title = "Seed: " + to_string(world->seed());
            else
                title = "Custom";

            title = title + "##regional";

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 10, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(161, 449), ImGuiCond_FirstUseEver);

            ImGui::Begin(title.c_str(), NULL, window_flags);

            ImGui::Text("World controls:");

            ImGui::PushItemWidth(100.0f);

            standardbutton("World map");

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Export options:");

            ImGui::PushItemWidth(100.0f);

            standardbutton("Regional maps");

            standardbutton("Area maps");

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Display map type:");

            ImGui::PushItemWidth(100.0f);

            for (const mapviewdefinition& definition : allmapviewdefinitions)
                standardbutton(definition.label);

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Other controls:");

            ImGui::PushItemWidth(100.0f);

            standardbutton("Properties");

            standardbutton("Appearance");

            ImGui::End();

            // Now the text box.

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 180, main_viewport->WorkPos.y + 542), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(1023, 140), ImGuiCond_FirstUseEver);

            title = "               ";

            ImGui::Begin(title.c_str(), NULL, window_flags);

            ImGui::End();

            // Now the additional element.

            showdeferredworkwindow(generatingregionpass, "Please wait...##generateregion ", "Generating region", ImVec2(main_viewport->WorkPos.x + 849, main_viewport->WorkPos.y + 364), ImVec2(200, 60));
        }
        else
        {
            newx = newx / 32;
            newy = newy / 32;

            newx = newx * 32;
            newy = newy * 32;

            newx = newx + 16;
            newy = newy + 16;

            region->setcentrex(newx);
            region->setcentrey(newy);

            mapview = relief;

            resetmapcache(regionalmaps);

            float progressstep = 1.0f / REGIONALCREATIONSTEPS;

            // Blank the regional map image first

            regionalreliefimage.create(regionalimagewidth, regionalimageheight, sf::Color::Black);

            updateTextureFromImage(regionalmaptexture, regionalreliefimage);
            regionalmap.setTexture(regionalmaptexture);

            // Now generate the regional map

            generateregionalmap(*world, *region, smalllake, island, peaks, riftblob, riftblobsize, 0, smudge, smallsmudge, squareroot);

            // Now draw the regional map image

            showRegionalMapView(mapview);

            // Sort out the minimap

            updateTextureFromImage(globalmaptexture, displayglobalreliefimage);
            minimap.setTexture(globalmaptexture);

            focused = 0;
            generatingregionpass = 0;

            screenmode = regionalmapscreen;
        }
    }
}
