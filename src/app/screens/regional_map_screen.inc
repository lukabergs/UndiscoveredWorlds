// Part of main() in app/main.cpp: included inside its frame loop, and uses its locals.

        if (screenmode == regionalmapscreen)
        {
            showglobaltemperaturechart = 0;
            showglobalrainfallchart = 0;
            showworldeditproperties = 0;

            areafromregional = 1;

            showsetsize = 0;
            newworld = 0;

            // Main controls.

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

            if (standardbutton("World map"))
            {
                focused = 0;
                poix = -1;
                poiy = -1;

                regionalreliefimage.create(regionalimagewidth, regionalimageheight, sf::Color::Black);

                redrawRegionalRelief();

                screenmode = globalmapscreen;

                showcolouroptions = 0;
                showworldproperties = 0;
                showregionalrainfallchart = 0;
                showregionaltemperaturechart = 0;
            }

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Export options:");

            ImGui::PushItemWidth(100.0f);

            if (standardbutton("Regional maps"))
            {
                openFileDialog(".png");

                exportingregionalmaps = 1;
            }

            if (standardbutton("Area maps"))
            {
                screenmode = exportareascreen;
                areafromregional = 1;
            }

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Display map type:");

            ImGui::PushItemWidth(100.0f);

            drawmapviewbuttons(showRegionalMapView);

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Other controls:");

            ImGui::PushItemWidth(100.0f);

            if (standardbutton("Properties"))
                toggle(showworldproperties);

            if (standardbutton("Appearance"))
                toggle(showcolouroptions);

            ImGui::End();

            drawregionalinfowindow(main_viewport, window_flags, *world, *region, focused == 1, poix, poiy, degree, cube, showregionaltemperaturechart, showregionalrainfallchart, showabout);

            if (window.hasFocus())
            {
                // Now check to see if the map has been clicked on.

                if (sf::Mouse::isButtonPressed(sf::Mouse::Left) && io.WantCaptureMouse == 0)
                {
                    sf::Vector2i mousepos = sf::Mouse::getPosition(window);

                    poix = mousepos.x - (int)regionalmapxpos;
                    poiy = mousepos.y - (int)regionalmapypos;

                    if (poix >= 0 && poix < regionalmapimagewidth && poiy >= 0 && poiy < regionalmapimageheight)
                    {
                        poix = poix + region->regwidthbegin();
                        poiy = poiy + region->regheightbegin();

                        focused = 1;

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

                    float mult = ((float)world->width() + 1.0f) / (float)DISPLAYMAPSIZEX;

                    float fpoix = (float)(mousepos.x - minimapxpos) * mult * 2.0f;
                    float fpoiy = (float)(mousepos.y - minimapypos) * mult * 2.0f;

                    int minipoix = (int)fpoix;
                    int minipoiy = (int)fpoiy;

                    if (minipoix >= 0 && minipoiy <= world->width() && minipoiy >= 0 && minipoiy <= world->height())// If the minimap has been clicked on.
                    {
                        screenmode = generatingregionscreen;
                        newx = minipoix;
                        newy = minipoiy;
                    }
                }

                // Check to see if the cursor keys have been pressed.

                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left))
                {
                    int regionalmapmove = REGIONALTILEWIDTH;

                    newx = region->centrex() - regionalmapmove;

                    if (newx > 0)
                    {
                        newx = newx / 32;
                        newx = newx * 32;
                        newx = newx + 16;
                    }

                    if (newx < regionalmapmove * 2)
                        newx = wrap(newx, world->width());

                    newy = region->centrey();

                    screenmode = generatingregionscreen;
                }

                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
                {
                    int regionalmapmove = REGIONALTILEWIDTH;

                    newx = region->centrex() + regionalmapmove;

                    newx = newx / 32;
                    newx = newx * 32;
                    newx = newx + 16;

                    if (newx > world->width())
                        newx = wrap(newx, world->width());

                    newy = region->centrey();

                    screenmode = generatingregionscreen;
                }

                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
                {
                    int regionalmapmove = REGIONALTILEHEIGHT;

                    newy = region->centrey() - regionalmapmove;

                    newy = newy / 32;
                    newy = newy * 32;
                    newy = newy + 16;

                    if (newy > regionalmapmove)
                    {
                        newx = region->centrex();
                        screenmode = generatingregionscreen;
                    }
                }

                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
                {
                    int regionalmapmove = REGIONALTILEHEIGHT;

                    newy = region->centrey() + regionalmapmove;

                    newy = newy / 32;
                    newy = newy * 32;
                    newy = newy + 16;

                    if (newy < world->height())
                    {
                        newx = region->centrex();
                        screenmode = generatingregionscreen;
                    }
                }
            }
        }
