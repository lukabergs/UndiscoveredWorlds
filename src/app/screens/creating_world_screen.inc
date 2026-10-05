// Part of main() in app/main.cpp: included inside its frame loop, and uses its locals.

        if (screenmode == creatingworldscreen)
        {
            if (showdeferredworkwindow(creatingworldpass, "Please wait...##createworld ", "Generating world...", ImVec2(main_viewport->WorkPos.x + 506, main_viewport->WorkPos.y + 171), ImVec2(173, 68)))
            {
                updatereport("Generating world from seed: " + to_string(world->seed()) + ":");
                updatereport("");
                beginworldgendebugrun(world->seed(), &worldgenerationdebug.options);

                if (worldgenerationdebug.options.visualizePlateTectonicsRealtime)
                {
                    setworldgenvisualizationcallback([&]()
                    {
                        invalidateGlobalMaps();
                        showGlobalMapView(relief);
                        window.clear();
                        window.draw(globalmap);
                        window.display();
                    });
                }
                else
                    clearworldgenvisualizationcallback();

                begintimedreporting();

                initialiseworld(*world);
                world->clear();

                changeworldproperties(*world);
                applyworlddimensions(*world, worldpropertycontrols);

                fast_srand(world->seed());

                resizeglobaldisplayforworld(*world, globaldisplay, highlightdisplay);

                invalidateGlobalMaps();

                GenerationScratch scratch;
                initializegenerationscratch(scratch);

                GenerationWorkbenchUiState ui;
                applytectonicsoptionsfromdebug(ui);
                ui.seaLevel = world->sealevel();

                GenerationExecutionContext context;
                context.dorivers = true;
                context.dolakes = true;
                context.dodeltas = true;
                context.socialoptions = socialoptionsfromdebug();
                context.landshape = landshape;
                context.chainland = chainland;
                context.smalllake = smalllake;
                context.largelake = largelake;

                string errormessage;
                runremaininggenerationstages(*world, scratch, ui, context, getgenerationstageindex(GenerationStageId::geosphere_initialization), &errormessage);

                if (errormessage.empty() == false)
                    updatereport(errormessage);

                // Now draw a new map

                mapview = relief;

                showGlobalMapView(mapview);

                endtimedreporting();
                endworldgendebugrun();
                updatereport("");
                updatereport("World generation completed.");
                updatereport("");

                focused = 0;

                creatingworldpass = 0;
                screenmode = movingtoglobalmapscreen;
                newworld = 1;
            }
        }
