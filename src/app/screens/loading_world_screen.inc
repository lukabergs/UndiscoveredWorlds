// Part of main() in app/main.cpp: included inside its frame loop, and uses its locals.

        if (screenmode == loadingworldscreen)
        {
            if (showdeferredworkwindow(loadingworldpass, "Please wait...##loadworld ", "Loading world...", ImVec2(main_viewport->WorkPos.x + 507, main_viewport->WorkPos.y + 173), ImVec2(173, 68)))
            {
                bool success=world->loadworld(filepathname);

                if (success == 0) // Failed to load
                {
                    screenmode = loadfailure;
                    loadingworld = 0;
                    loadingworldpass = 0;
                }
                else
                {
                    syncworldpropertycontrols(*world, worldpropertycontrols);
                    syncappearancesettings(*world, appearance);

                    // Now draw the new world.

                    resizeglobaldisplayforworld(*world, globaldisplay, highlightdisplay);

                    resetmapcache(globalmaps);

                    mapview = relief;
                    showGlobalMapView(mapview);

                    focused = 0;
                    newworld = 1;

                    filepathname = "";
                    filepath = "";

                    screenmode = globalmapscreen;
                    loadingworld = 0;
                    loadingworldpass = 0;
                }
            }
        }
