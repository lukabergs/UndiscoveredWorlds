// Part of main() in app/main.cpp: included inside its frame loop, and uses its locals.

        if (showcolouroptions && screenmode != settingsloadfailure)
        {
            const int colouralign = 390;
            const int otheralign = 360;

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 300, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(780, 640), ImGuiCond_FirstUseEver);

            ImGui::Begin("Map appearance", NULL, window_flags);

            ImGui::Checkbox("Show landmass outline", &showmapoutline);
            ImGui::SameLine(260.0f);
            ImGui::ColorEdit3("Outline", (float*)&outlinecolour);
            ImGui::SameLine(520.0f);
            ImGui::TextUnformatted("Applies to every map view.");

            if (ImGui::BeginTabBar("MapAppearanceTabs"))
            {
                for (const mapviewdefinition& definition : allmapviewdefinitions)
                {
                    if (!ImGui::BeginTabItem(definition.label))
                        continue;

                    drawmapviewappearancetab(definition, *world, appearance, selectedgradientstops, colouralign, otheralign);

                    ImGui::EndTabItem();
                }

                ImGui::EndTabBar();
            }

            ImGui::SetCursorPos(ImVec2(140.0f, 600.0f));

            if (ImGui::Button("Save", ImVec2(120.0f, 0.0f)))
            {
                openFileDialog(".uws");

                savingsettings = 1;
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Save map appearance settings.");

            ImGui::SameLine();
            if (ImGui::Button("Load", ImVec2(120.0f, 0.0f)))
            {
                openFileDialog(".uws");

                loadingsettings = 1;
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Load map appearance settings.");

            ImGui::SameLine();
            if (ImGui::Button("Default", ImVec2(120.0f, 0.0f)))
            {
                initialisemapcolours(*world);
                syncappearancesettings(*world, appearance);
                refreshhighlightdisplay(*world, highlightdisplay);

                colourschanged = 1;
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Restore the default map appearance settings.");

            ImGui::SameLine();
            if (ImGui::Button("Close", ImVec2(120.0f, 0.0f)))
            {
                showcolouroptions = 0;
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Close the map appearance panel.");

            ImGui::End();
        }
