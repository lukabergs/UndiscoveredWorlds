#include "app/application.hpp"

void Application::updateexportareascreen()
{
    ImGuiIO& io = ImGui::GetIO();

    if (screenmode == exportareascreen)
    {
        showglobaltemperaturechart = 0;
        showglobalrainfallchart = 0;
        showregionaltemperaturechart = 0;
        showregionalrainfallchart = 0;
        showworldproperties = 0;
        showworldeditproperties = 0;
        showcolouroptions = 0;
        newworld = 0;

        showsetsize = 0;

        // Work out the size of the currently selected area.

        int totalregions = 0;
        int maxtotalregions = 100; // Query area maps larger than this.

        if (areanex != -1)
        {
            float regiontilewidth = (float)REGIONALTILEWIDTH; //30;
            float regiontileheight = (float)REGIONALTILEHEIGHT; //30; // The width and height of the visible regional map, in tiles.

            int regionwidth = (int)regiontilewidth * 16;
            int regionheight = (int)regiontileheight * 16; // The width and height of the visible regional map, in pixels.

            int newareanwx = areanwx; // This is because the regions we'll be making will start to the north and west of the defined area.
            int newareanwy = areanwy;

            int newareanex = areanex;
            int newareaney = areaney;

            int newareaswx = areaswx;
            int newareaswy = areaswy;

            int newareasey = areasey;

            newareanwx = newareanwx / (int)regiontilewidth;
            newareanwx = newareanwx * (int)regiontilewidth;

            newareanwy = newareanwy / (int)regiontileheight;
            newareanwy = newareanwy * (int)regiontileheight;

            newareaswx = newareanwx;
            newareaney = newareanwy;

            float areatilewidth = (float)(newareanex - newareanwx);
            float areatileheight = (float)(newareasey - newareaney);

            int areawidth = (int)(areatilewidth * 16.0f);
            int areaheight = (int)(areatileheight * 16.0f);

            float fregionswide = areatilewidth / regiontilewidth;
            float fregionshigh = areatileheight / regiontileheight;

            int regionswide = (int)fregionswide;
            int regionshigh = (int)fregionshigh;

            if (regionswide != fregionswide)
                regionswide++;

            if (regionshigh != fregionshigh)
                regionshigh++;

            totalregions = regionswide * regionshigh; // This is how many regional maps we're going to have to do.
        }

        // Main controls.

        ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 10, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(161, 138), ImGuiCond_FirstUseEver);

        ImGui::Begin("Export custom area", NULL, window_flags);

        ImGui::PushItemWidth(100.0f);

        if (standardbutton("Clear selection"))
        {
            areanex = -1;
            areaney = -1;
            areasex = -1;
            areasey = -1;
            areaswx = -1;
            areaswy = -1;
            areanwx = -1;
            areanwy = -1;
        }

        if (standardbutton("Export maps"))
        {
            if (areanex != -1)
            {
                if (totalregions <= maxtotalregions)
                {
                    openFileDialog(".png");

                    exportingareamaps = 1;
                }
                else
                    showareawarning = 1;
            }
        }

        ImGui::Text("  ");

        if (standardbutton("Cancel"))
        {
            areanex = -1;
            areaney = -1;
            areasex = -1;
            areasey = -1;
            areaswx = -1;
            areaswy = -1;
            areanwx = -1;
            areanwy = -1;

            if (areafromregional == 1)
                screenmode = regionalmapscreen;
            else
                screenmode = globalmapscreen;
        }

        ImGui::End();

        // Now the text box.

        ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 180, main_viewport->WorkPos.y + 542), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(1023, 151), ImGuiCond_FirstUseEver);
        string title = "            ";
        string areatext = "Total regions: " + to_string(totalregions); // 117 OK; 143 not OK; 132 not OK; 121 not OK; 112 OK; 116 OK; 110

        ImGui::Begin(title.c_str(), NULL, window_flags);
        ImGui::PushItemWidth((float)world->width() / 2.0f);
        areatext="This screen allows you to export maps at the same scale as the regional map, but of larger areas.\n\nClick on the map to pick the corners of the area you want to export. You can re-select corners to fine-tune the area.\n\nWhen you are done, click on 'export maps'. The program will ask you to specify the filename under which to save the maps, and then create them.";

        ImGui::Text(areatext.c_str(), world->width() / 2);
        ImGui::End();

        // Now check to see if the map has been clicked on.

        if (window.hasFocus() && sf::Mouse::isButtonPressed(sf::Mouse::Left) && io.WantCaptureMouse == 0)
        {
            sf::Vector2i mousepos = sf::Mouse::getPosition(window);

            float mult = ((float)world->width() + 1.0f) / (float)DISPLAYMAPSIZEX;

            float xpos = ((float)mousepos.x - (float)globalmapxpos) * mult;
            float ypos = ((float)mousepos.y - (float)globalmapypos) * mult;

            poix = (int)xpos;
            poiy = (int)ypos;

            if (poix >= 0 && poix <= world->width() && poiy >= 0 && poiy <= world->height())
            {
                if (areaswx == -1) // If we don't have any corners yet
                {
                    areanex = poix + 2;
                    areaney = poiy;
                    areasex = poix + 2;
                    areasey = poiy + 2;
                    areaswx = poix;
                    areaswy = poiy + 2;
                    areanwx = poix;
                    areanwy = poiy;
                }
                else // If we do have all corners
                {
                    int nedistx = areanex - poix;
                    int nedisty = areaney - poiy;

                    int sedistx = areasex - poix;
                    int sedisty = areasey - poiy;

                    int swdistx = areaswx - poix;
                    int swdisty = areaswy - poiy;

                    int nwdistx = areanwx - poix;
                    int nwdisty = areanwy - poiy;

                    float nedist = (float)sqrt(nedistx * nedistx + nedisty * nedisty);
                    float sedist = (float)sqrt(sedistx * sedistx + sedisty * sedisty);
                    float swdist = (float)sqrt(swdistx * swdistx + swdisty * swdisty);
                    float nwdist = (float)sqrt(nwdistx * nwdistx + nwdisty * nwdisty);

                    short id = 2; // This will identify which corner this is going to be.
                    float mindist = nedist;

                    if (sedist < mindist)
                    {
                        id = 4;
                        mindist = sedist;
                    }

                    if (swdist < mindist)
                    {
                        id = 6;
                        mindist = swdist;
                    }

                    if (nwdist < mindist)
                        id = 8;

                    if (id == 2)
                    {
                        areanex = poix;
                        areaney = poiy;
                        areasex = poix;
                        areanwy = poiy;
                    }

                    if (id == 4)
                    {
                        areasex = poix;
                        areasey = poiy;
                        areanex = poix;
                        areaswy = poiy;
                    }

                    if (id == 6)
                    {
                        areaswx = poix;
                        areaswy = poiy;
                        areanwx = poix;
                        areasey = poiy;
                    }

                    if (id == 8)
                    {
                        areanwx = poix;
                        areanwy = poiy;
                        areaswx = poix;
                        areaney = poiy;
                    }

                    if (areanwx > areanex)
                        swap(areanwx, areanex);

                    if (areaswx > areasex)
                        swap(areaswx, areasex);

                    if (areanwy > areaswy)
                        swap(areanwy, areaswy);

                    if (areaney > areasey)
                        swap(areaney, areasey);
                }
            }
        }
    }
}
