#include "app/application.hpp"

void Application::updateimportscreen()
{
    ImGuiIO& io = ImGui::GetIO();

    if (screenmode == importscreen)
    {
        showglobaltemperaturechart = 0;
        showglobalrainfallchart = 0;
        showregionaltemperaturechart = 0;
        showregionalrainfallchart = 0;
        showworldproperties = 0;

        showcolouroptions = 0;
        showworldproperties = 0;

        if (generationworkbench.active)
        {
            const ImportedClimateMaps* workbenchclimate = (importedclimatemaps.hasTemperature || importedclimatemaps.hasPrecipitation) ? &importedclimatemaps : nullptr;
            const planet& displayedworkbenchworld = getworkbenchdisplayworld(generationworkbench, *world);
            const GenerationWorkbenchPanelResult workbenchpanel = drawgenerationworkbenchpanel(main_viewport, window_flags, displayedworkbenchworld, generationworkbench);

            GenerationExecutionContext workbenchcontext;
            workbenchcontext.dodeltas = currentdeltas;
            workbenchcontext.socialoptions = socialoptionsfromworldcontrols();
            workbenchcontext.landshape = landshape;
            workbenchcontext.chainland = chainland;
            workbenchcontext.smalllake = smalllake;
            workbenchcontext.largelake = largelake;
            workbenchcontext.okmountains = &OKmountains;
            workbenchcontext.importedClimate = workbenchclimate;

            auto refreshworkbenchdisplay = [&]()
            {
                invalidateGlobalMaps();

                if (workbenchfinished(generationworkbench))
                    showGlobalMapView(relief);
                else
                    showGlobalMapView(preferredworkbenchmapview(getcurrentgenerationstage(generationworkbench).id));
            };

            auto previewworkbenchstage = [&](std::string& errormessage)
            {
                if (generationworkbench.active == false || workbenchfinished(generationworkbench))
                    return;

                runworkbenchaction([&]()
                {
                    previewcurrentgenerationstage(generationworkbench, *world, workbenchcontext, &errormessage);
                });
            };

            if (getcurrentgenerationstage(generationworkbench).id == GenerationStageId::inland_seas
                && window.hasFocus()
                && ImGui::IsMouseClicked(ImGuiMouseButton_Left)
                && io.WantCaptureMouse == 0)
            {
                sf::Vector2i mousepos = sf::Mouse::getPosition(window);
                const float mult = static_cast<float>(displayedworkbenchworld.width() + 1) / static_cast<float>(DISPLAYMAPSIZEX);
                const int clickedx = static_cast<int>((static_cast<float>(mousepos.x) - globalmapxpos) * mult);
                const int clickedy = static_cast<int>((static_cast<float>(mousepos.y) - globalmapypos) * mult);

                if (clickedx >= 0 && clickedx <= displayedworkbenchworld.width() && clickedy >= 0 && clickedy <= displayedworkbenchworld.height())
                {
                    if (cycleinlandseacomponentpolicyat(generationworkbench, displayedworkbenchworld.grid().cellof(clickedx, clickedy)))
                        refreshworkbenchdisplay();
                }
            }

            string workbencherror;

            if (workbenchpanel.displayChanged)
            {
                refreshworkbenchdisplay();
            }

            if (workbenchpanel.discardPreviewRequested)
            {
                clearworkbenchpreview(generationworkbench);
                refreshworkbenchdisplay();
            }

            if (workbenchpanel.controlsChanged)
            {
                clearworkbenchpreview(generationworkbench);
                previewworkbenchstage(workbencherror);
                refreshworkbenchdisplay();
            }
            else if (workbenchpanel.recomputeRequested)
            {
                clearworkbenchpreview(generationworkbench);
                previewworkbenchstage(workbencherror);
                refreshworkbenchdisplay();
            }
            else if (workbenchpanel.previewModeChanged)
            {
                if (generationworkbench.ui.previewEnabled
                    && generationworkbench.previewAvailable == false
                    && getcurrentgenerationstage(generationworkbench).id != GenerationStageId::plate_tectonics)
                {
                    previewworkbenchstage(workbencherror);
                }

                refreshworkbenchdisplay();
            }

            if (workbenchpanel.applyRequested)
            {
                runworkbenchaction([&]()
                {
                    applycurrentgenerationstage(generationworkbench, *world, workbenchcontext, &workbencherror);
                });

                if (generationworkbench.active && workbenchfinished(generationworkbench) == false)
                    previewworkbenchstage(workbencherror);

                refreshworkbenchdisplay();
            }

            if (workbenchpanel.skipRequested)
            {
                skipcurrentgenerationstage(generationworkbench, *world);

                if (generationworkbench.active && workbenchfinished(generationworkbench) == false)
                    previewworkbenchstage(workbencherror);

                refreshworkbenchdisplay();
            }

            if (workbenchpanel.resetPreviewRequested)
            {
                clearworkbenchpreview(generationworkbench);
                refreshworkbenchdisplay();
            }

            if (workbenchpanel.backRequested)
            {
                if (stepbackworkbenchsession(generationworkbench, *world))
                {
                    if (generationworkbench.ui.previewEnabled)
                        previewworkbenchstage(workbencherror);

                    refreshworkbenchdisplay();
                }
            }

            if (workbenchpanel.abortRequested)
            {
                seedentry = generationworkbench.ui.tectonicSeed;
                const bool proceduralsession = generationworkbench.procedural;
                resetworkbenchsession(generationworkbench);

                if (proceduralsession)
                {
                    showworldgenerationoptions = 0;
                    screenmode = createworldscreen;
                }
                else
                    refreshworkbenchdisplay();
            }

            if (workbencherror.empty() == false)
                updatereport(workbencherror);

            if (workbenchfinished(generationworkbench))
            {
                const bool proceduralsession = generationworkbench.procedural;
                resetworkbenchsession(generationworkbench);

                if (proceduralsession)
                {
                    invalidateGlobalMaps();
                    showGlobalMapView(relief);
                    screenmode = movingtoglobalmapscreen;
                }
            }
        }
        else
        {
        // Main controls.

        ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 10, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(161, 735), ImGuiCond_FirstUseEver);

        ImGui::Begin("Custom##custom", NULL, window_flags);

        ImGui::SetNextItemWidth(0);
        ImGui::Text("Value mode:");

        int importvaluemode = static_cast<int>(mapimportsettings.valueMode);
        const char* importmodeitems[] = { "Red", "Strip" };

        ImGui::PushItemWidth(100.0f);

        if (ImGui::Combo("##importvaluemode", &importvaluemode, importmodeitems, IM_ARRAYSIZE(importmodeitems)))
            mapimportsettings.valueMode = static_cast<MapImportValueMode>(importvaluemode);

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("PNG imports can read the red channel directly or decode colours through a 1-pixel-tall gradient strip. TIFF imports always use the stored uint16 sample values.");

        if (mapimportsettings.valueMode == MapImportValueMode::gradientStrip)
        {
            if (standardbutton("Scale strip"))
            {
                openFileDialog(kGradientStripFilter);
                importinggradientstrip = 1;
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Load a 1-pixel-tall PNG strip. The leftmost colour maps to the minimum value, and each pixel to the right adds the chosen increment.");

            ImGui::InputFloat("Minimum", &mapimportsettings.gradientMinimum, 1.0f, 10.0f, "%.3f");

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Value assigned to the leftmost colour in the gradient strip.");

            ImGui::InputFloat("Step", &mapimportsettings.gradientIncrement, 0.1f, 1.0f, "%.3f");

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Amount added for each colour step to the right in the gradient strip.");

            const char* stripstatus = mapimportsettings.gradientStripPath.empty() ? "Strip: none" : "Strip: loaded";
            ImGui::TextUnformatted(stripstatus);

            if (!mapimportsettings.gradientStripPath.empty() && ImGui::IsItemHovered())
                ImGui::SetTooltip(mapimportsettings.gradientStripPath.c_str());
        }

        ImGui::Dummy(ImVec2(0.0f, linespace));

        ImGui::Text("Import:");

        if (standardbutton("Land map"))
        {
            openFileDialog(kRasterImportFilter);

            importinglandmap = 1;
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("TIFF mode: uint16 values are metres above sea level, with 0 for sea. PNG red mode: red 0 is sea and higher values are elevation above sea level in 10-metre steps. PNG strip mode: decoded values are treated as elevation above sea level.");

        if (standardbutton("Sea map"))
        {
            openFileDialog(kRasterImportFilter);

            importingseamap = 1;
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("TIFF mode: uint16 values are metres below sea level, with 0 for land. PNG red mode: red 0 is land and higher values are depth below sea level in 50-metre steps. PNG strip mode: decoded values are treated as depth below sea level.");

        if (standardbutton("Mountains"))
        {
            openFileDialog(kRasterImportFilter);

            importingmountainsmap = 1;
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("TIFF mode: uint16 values are imported mountain height directly. PNG red mode: red values scale the imported mountain height. PNG strip mode: decoded values are used directly.");

        if (standardbutton("Volcanoes"))
        {
            openFileDialog(kVolcanoImportFilter);

            importingvolcanoesmap = 1;
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("PNG only. Red or strip decoding sets volcano height. Green=0 for shield volcano, or higher for stratovolcano. Blue=0 for extinct, or higher for active.");

        if (standardbutton("Temperature"))
        {
            openFileDialog(kRasterImportFilter);

            importingtemperaturemap = 1;
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("TIFF mode: uint16 values are imported directly as mean annual temperature. PNG red mode maps 0..255 linearly to -60 C..+60 C. PNG strip mode uses the decoded value directly.");

        if (standardbutton("Precipitation"))
        {
            openFileDialog(kRasterImportFilter);

            importingprecipitationmap = 1;
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("TIFF mode: uint16 values are imported directly as mean precipitation. PNG red mode maps 0..255 linearly to 0..1020. PNG strip mode uses the decoded value directly.");

        ImGui::Dummy(ImVec2(0.0f, linespace));

        ImGui::SetNextItemWidth(0);

        ImGui::Text("Generate terrain:");

        ImGui::PushItemWidth(100.0f);

        if (standardbutton("Terrain"))
            toggle(showterrainchooser);

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Create global terrain with the tectonic pipeline and continental shelves.");

        ImGui::Dummy(ImVec2(0.0f, linespace));

        ImGui::SetNextItemWidth(0);

        ImGui::Text("Generate elements:");

        ImGui::PushItemWidth(100.0f);

        if (standardbutton("Coastlines"))
        {
            runprocess("terrain.coastline_refinement");

            // Now redraw the map.

            resetmapcache(globalmaps);

            mapview = relief;
            showGlobalMapView(mapview);
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Disrupt straight edges on coastlines.");

        if (standardbutton("Shelves"))
        {
            generateimportedshelves(*world);

            // Now redraw the map.

            resetmapcache(globalmaps);

            mapview = relief;
            showGlobalMapView(mapview);
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Generate continental shelves.");

        if (standardbutton("Oceanic ridges"))
        {
            generateimportedoceanridges(*world);

            // Now redraw the map.

            resetmapcache(globalmaps);

            mapview = relief;
            showGlobalMapView(mapview);
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Generate mid-ocean ridges.");

        if (standardbutton("Sea bed"))
        {
            generateimportedseabed(*world);

            // Now redraw the map.

            resetmapcache(globalmaps);

            mapview = relief;
            showGlobalMapView(mapview);
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Create random depth variation across the oceans.");

        if (standardbutton("Land elevation"))
        {
            generateimportedlandelevation(*world);

            // Now redraw the map.

            resetmapcache(globalmaps);

            mapview = relief;
            showGlobalMapView(mapview);
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Generate random elevation variation across the land.");

        if (standardbutton("Mountains##2"))
        {
            generateimportedmountainchains(*world, landshape, chainland);

            // Now redraw the map.

            resetmapcache(globalmaps);

            mapview = relief;
            showGlobalMapView(mapview);
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Generate random mountain ranges.");

        if (standardbutton("Hills"))
        {
            generateimportedmountainchains(*world, landshape, chainland);

            // Now redraw the map.

            resetmapcache(globalmaps);

            mapview = relief;
            showGlobalMapView(mapview);
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Generate random ranges of hills.");

        if (standardbutton("Craters"))
        {
            generateimportedcraters(*world, squareroot);

            // Now redraw the map.

            resetmapcache(globalmaps);

            mapview = relief;
            showGlobalMapView(mapview);
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Generate random craters.");

        ImGui::Dummy(ImVec2(0.0f, linespace));

        ImGui::SetNextItemWidth(0);

        ImGui::Text("Other controls:");

        ImGui::PushItemWidth(100.0f);

        if (standardbutton("Properties"))
            toggle(showworldeditproperties);

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Adjust the world's properties.");

        if (standardbutton("Workbench"))
        {
            initializeimportedworkbenchsession(generationworkbench, *world, importedclimatemaps.hasTemperature || importedclimatemaps.hasPrecipitation);
            previewcurrentgenerationstage(generationworkbench, *world, GenerationExecutionContext{ currentdeltas, socialoptionsfromworldcontrols(), landshape, chainland, smalllake, largelake, &OKmountains, (importedclimatemaps.hasTemperature || importedclimatemaps.hasPrecipitation) ? &importedclimatemaps : nullptr }, nullptr);
            invalidateGlobalMaps();
            showGlobalMapView(preferredworkbenchmapview(getcurrentgenerationstage(generationworkbench).id));
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Open the staged generation workbench for imported or hand-edited terrain.");

        if (standardbutton("Finish"))
        {
            screenmode = completingimportscreen;
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Calculate climates, lakes, and rivers, and finish the world.");


        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Registers this world as a climate benchmark run with its diagnostics and maps.");

        if (standardbutton("Cancel"))
        {
            brandnew = 1;
            seedentry = 0;

            screenmode = createworldscreen;
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Return to the world creation screen.");

        ImGui::End();

        // Now the text box.

        ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 180, main_viewport->WorkPos.y + 542), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(1023, 155), ImGuiCond_FirstUseEver);

        string title = "            ";

        string importtext = "You can use the 'import' buttons to load in your own maps. These must be " + formatnumber(world->width() + 1) + " x " + formatnumber(world->height() + 1) + " pixels.\nLand, sea, mountains, temperature, and precipitation imports accept grayscale uint16 .tif/.tiff files and also support the older PNG modes. Volcanoes and gradient strips still use PNG because they rely on colour channels.\nUse 'Terrain' to open the procedural terrain workbench, or 'Workbench' to step through terrain and climate stages on imported or hand-edited maps.\nYou can use the other buttons to tweak imported maps directly, and the 'Properties' panel to change world settings such as temperatures or rainfall.\nWhen you are done, click 'Finish' to run the default completion flow.";

        ImGui::Begin(title.c_str(), NULL, window_flags);
        ImGui::PushItemWidth((float)(world->width() / 2));
        ImGui::Text(importtext.c_str(), world->width() / 2);
        ImGui::End();
        }
    }
}
