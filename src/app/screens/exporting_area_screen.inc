// Part of main() in app/main.cpp: included inside its frame loop, and uses its locals.

        if (screenmode == exportingareascreen)
        {
            if (showdeferredworkwindow(exportingareapass, "Please wait...##exportarea ", "Generating area maps...", ImVec2(main_viewport->WorkPos.x + 300, main_viewport->WorkPos.y + 200), ImVec2(173, 68)))
            {
                int oldregionalcentrex = region->centrex();
                int oldregionalcentrey = region->centrey();

                mapviewenum oldmapview = mapview;

                float regiontilewidth = (float)REGIONALTILEWIDTH; //30;
                float regiontileheight = (float)REGIONALTILEHEIGHT; //30; // The width and height of the visible regional map, in tiles.

                int regionwidth = (int)(regiontilewidth * 16.0f);
                int regionheight = (int)(regiontileheight * 16.0f); // The width and height of the visible regional map, in pixels.

                int origareanwx = areanwx; // This is because the regions we'll be making will start to the north and west of the defined area.
                int origareanwy = areanwy;

                int origareanex = areanex;
                int origareaney = areaney;

                int origareaswx = areaswx;
                int origareaswy = areaswy;

                areanwx = areanwx / (int)regiontilewidth;
                areanwx = areanwx * (int)regiontilewidth;

                areanwy = areanwy / (int)regiontileheight;
                areanwy = areanwy * (int)regiontileheight;

                areaswx = areanwx;
                areaney = areanwy;

                int woffset = (origareanwx - areanwx) * 16;
                int noffset = (origareanwy - areanwy) * 16;

                float areatilewidth = (float)(areanex - areanwx);
                float areatileheight = (float)(areasey - areaney);

                int areawidth = (int)areatilewidth * 16;
                int areaheight = (int)areatileheight * 16;

                float imageareatilewidth = (float)(origareanex - origareanwx);
                float imageareatileheight = (float)(areasey - origareaney);

                int imageareawidth = (int)imageareatilewidth * 16;
                int imageareaheight = (int)imageareatileheight * 16;

                float fregionswide = areatilewidth / regiontilewidth;
                float fregionshigh = areatileheight / regiontileheight;

                int regionswide = (int)fregionswide;
                int regionshigh = (int)fregionshigh;

                if (regionswide != fregionswide)
                    regionswide++;

                if (regionshigh != fregionshigh)
                    regionshigh++;

                int totalregions = regionswide * regionshigh; // This is how many regional maps we're going to have to do.

                if (areafromregional == 1)
                    totalregions++; // Because we'll have to redo the regional map we came from.

                initialiseregion(*world, *region); // We'll do all this using the same region object as usual. We could create a new region object for it, but that seems to lead to inexplicable crashes, so we won't.

                // Now we need to prepare the images that we're going to copy the regional maps onto.

                sf::Image areareliefimage;
                sf::Image areaelevationimage;
                sf::Image areatemperatureimage;
                sf::Image areaprecipitationimage;
                sf::Image areaclimateimage;
                sf::Image areabiomeimage;
                sf::Image areariversimage;
                sf::Image areageologyimage;
                sf::Image areabasinsimage;
                sf::Image areaerosionimage;
                sf::Image areadepositionimage;
                sf::Image areafertilityimage;
                sf::Image arearesourcesimage;
                sf::Image areasuitabilityimage;
                sf::Image areasettlementsimage;
                sf::Image areapopulationimage;
                sf::Image areainfrastructureimage;
                sf::Image areapolitiesimage;
                sf::Image areatradeimage;
                sf::Image arearecentconflictimage;

                auto getAreaExportImage = [&](mapviewenum view) -> sf::Image&
                {
                    switch (view)
                    {
                    case relief: return areareliefimage;
                    case elevation: return areaelevationimage;
                    case temperature: return areatemperatureimage;
                    case precipitation: return areaprecipitationimage;
                    case climate: return areaclimateimage;
                    case biomes: return areabiomeimage;
                    case rivers: return areariversimage;
                    case geology: return areageologyimage;
                    case basins: return areabasinsimage;
                    case erosion: return areaerosionimage;
                    case deposition: return areadepositionimage;
                    case fertility: return areafertilityimage;
                    case resources: return arearesourcesimage;
                    case suitability: return areasuitabilityimage;
                    case settlements_map: return areasettlementsimage;
                    case population: return areapopulationimage;
                    case infrastructure_map: return areainfrastructureimage;
                    case polities_map: return areapolitiesimage;
                    case trade_map: return areatradeimage;
                    case recent_conflict_map: return arearecentconflictimage;
                    }

                    return areareliefimage;
                };

                for (const mapviewdefinition& definition : allmapviewdefinitions)
                    getAreaExportImage(definition.view).create(imageareawidth + 1, imageareaheight + 1, sf::Color::Black);

                // Now it's time to make the regions, one by one, generate their maps, and copy those maps over onto the area images.

                for (int i = 0; i < regionswide; i++)
                {
                    int centrex = i * (int)regiontilewidth + ((int)regiontilewidth / 2) + areanwx;

                    for (int j = 0; j < regionshigh; j++)
                    {
                        int centrey = j * (int)regiontileheight + ((int)regiontileheight / 2) + areanwy;

                        // First, create the new region.

                        region->setcentrex(centrex);
                        region->setcentrey(centrey);

                        generateregionalmap(*world, *region, smalllake, island, peaks, riftblob, riftblobsize, 0, smudge, smallsmudge,squareroot);

                        // Now generate the maps.

                        resetmapcache(regionalmaps);

                        rebuildAllRegionalMaps();

                        // Now copy those maps into the images that will be exported.

                        const int destx = i * regionwidth - woffset;
                        const int desty = j * regionheight - noffset;

                        auto copyRegionalLayer = [&](sf::Image& source, sf::Image& destination)
                        {
                            const int copyx0 = max(1, destx);
                            const int copyy0 = max(1, desty);
                            const int copyx1 = min(imageareawidth - 1, destx + regionalimagewidth - 1);
                            const int copyy1 = min(imageareaheight - 1, desty + regionalimageheight - 1);

                            if (copyx0 > copyx1 || copyy0 > copyy1)
                                return;

                            destination.copy(source, copyx0, copyy0, sf::IntRect(copyx0 - destx, copyy0 - desty, copyx1 - copyx0 + 1, copyy1 - copyy0 + 1), false);
                        };

                        for (const mapviewdefinition& definition : allmapviewdefinitions)
                            copyRegionalLayer(getmapimage(regionalmaps, definition.view), getAreaExportImage(definition.view));
                    }
                }

                region->setcentrex(oldregionalcentrex); // Move the region back to where it started.
                region->setcentrey(oldregionalcentrey);

                if (areafromregional == 1) // If we're going to go back to the regional map, we need to redo it.
                    generateregionalmap(*world, *region, smalllake, island, peaks, riftblob, riftblobsize, 0, smudge, smallsmudge, squareroot);

                // Now just save the images.

                filepathname.resize(filepathname.size() - 4);

                for (const mapviewdefinition& definition : allmapviewdefinitions)
                    getAreaExportImage(definition.view).saveToFile(filepathname + " " + definition.exportstem + ".png");

                // Clean up.

                filepathname = "";
                filepath = "";

                areanex = -1;
                areaney = -1;
                areasex = -1;
                areasey = -1;
                areaswx = -1;
                areaswy = -1;
                areanwx = -1;
                areanwy = -1;

                mapview = oldmapview;

                resetmapcache(regionalmaps);

                exportingareapass = 0;

                if (areafromregional == 1)
                    screenmode = regionalmapscreen;
                else
                    screenmode = globalmapscreen;
            }
        }
