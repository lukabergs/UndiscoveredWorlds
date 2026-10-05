// Part of main() in app/main.cpp: included inside its frame loop, and uses its locals.

        if (showworldproperties)
        {
            int rightalign = 290;

            int topleftfigures = 105;
            int bottomleftfigures = 225;
            int toprightfigures = 410;
            int bottomrightfigures = 460;

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 265, main_viewport->WorkPos.y + 60), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(640, 500), ImGuiCond_FirstUseEver);

            ImGui::Begin("World properties", NULL, window_flags & ~ImGuiWindowFlags_NoResize);

            std::array<int, 32> climateareas = {};
            int totallandarea = 0;

            for (int i = 0; i <= world->width(); i++)
            {
                for (int j = 0; j <= world->height(); j++)
                {
                    short thisclimate = world->climate(i, j);

                    if (thisclimate > 0 && thisclimate < static_cast<short>(climateareas.size()))
                    {
                        climateareas[thisclimate]++;
                        totallandarea++;
                    }
                }
            }

            if (ImGui::BeginTabBar("WorldPropertiesTabs"))
            {
                if (ImGui::BeginTabItem("Overview"))
                {

            string sizeinfo = "Size:";

            string sizevalue = "";

            if (world->size() == 0)
                sizevalue = "small";

            if (world->size() == 1)
                sizevalue = "medium";

            if (world->size() == 2)
                sizevalue = "large";

            stringstream ss3;
            ss3 << fixed << setprecision(5) << world->eccentricity();

            string eccentricityinfo = "Eccentricity:";
            string eccentricityvalue = ss3.str();

            stringstream ss4;
            ss4 << fixed << setprecision(2) << world->gravity();

            string gravityinfo = "Surface gravity: ";
            string gravityvalue = ss4.str() + "g";

            string perihelioninfo = "Perihelion:";
            string perihelionvalue = "";

            if (world->perihelion() == 0)
                perihelionvalue = "January";

            if (world->perihelion() == 1)
                perihelionvalue = "July";

            stringstream ss5;
            ss5 << fixed << setprecision(2) << world->lunar();

            string lunarinfo = "Lunar pull:";
            string lunarvalue = ss5.str();

            string typeinfo = "Category: ";
            string typevalue = "tectonic";

            string rotationinfo = "Rotation:";
            string rotationvalue = "";

            if (world->rotation())
                rotationvalue = "west to east";
            else
                rotationvalue = "east to west";

            stringstream ss;
            ss << fixed << setprecision(2) << world->tilt();

            string tiltinfo = "Obliquity:";
            string tiltvalue = ss.str() + degree;

            stringstream ss2;
            ss2 << fixed << setprecision(2) << world->tempdecrease();

            string tempdecreaseinfo = "Heat decrease per vertical km:";
            string tempdecreasevalue = ss2.str() + degree;

            string northpoleinfo = "North pole adjustment:";
            string northpolevalue = to_string(world->northpolaradjust()) + degree;

            string southpoleinfo = "South pole adjustment:";
            string southpolevalue = to_string(world->southpolaradjust()) + degree;

            string averageinfo = "Average global temperature:";
            string averagevalue = to_string(world->averagetemp()) + degree;

            stringstream ss6;
            ss6 << fixed << setprecision(2) << world->waterpickup();

            string moistureinfo = "Moisture pickup rate:";
            string moisturevalue = ss6.str();

            string glacialinfo = "Glaciation temperature:";
            string glacialvalue=to_string(world->glacialtemp())+degree;

            ImGui::Text(sizeinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Size of planet. (Earth: large; Mars: medium; Moon: small)");

            ImGui::SameLine((float)topleftfigures);
            ImGui::Text(sizevalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(gravityinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Affects mountain and valley sizes. (Earth: 1.00g)");

            ImGui::SameLine((float)toprightfigures);
            ImGui::Text(gravityvalue.c_str());

            ImGui::Text(typeinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Unified terrain generator in use.");

            ImGui::SameLine((float)topleftfigures);
            ImGui::Text(typevalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(lunarinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Affects tides and coastal regions. (Earth: 1.00)");

            ImGui::SameLine((float)toprightfigures);
            ImGui::Text(lunarvalue.c_str());

            ImGui::Text(rotationinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Affects weather patterns. (Earth: west to east)");

            ImGui::SameLine((float)topleftfigures);
            ImGui::Text(rotationvalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(tiltinfo.c_str());

            string tilttip = "Affects seasonal variation in temperature. (Earth: 22.5" + degree + ")";

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(tilttip.c_str());

            ImGui::SameLine((float)toprightfigures);
            ImGui::Text(tiltvalue.c_str());

            ImGui::Text(eccentricityinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("How elliptical the orbit is. (Earth: 0.0167)");

            ImGui::SameLine((float)topleftfigures);
            ImGui::Text(eccentricityvalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(perihelioninfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("When the planet is closest to the sun. (Earth: January)");

            ImGui::SameLine((float)toprightfigures);
            ImGui::Text(perihelionvalue.c_str());

            ImGui::Text("   ");

            ImGui::Text(tempdecreaseinfo.c_str());

            string tempdecreasetip = "Affects how much colder it gets higher up. (Earth: 6.5" + degree + ")";

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(tempdecreasetip.c_str());

            ImGui::SameLine((float)bottomleftfigures);
            ImGui::Text(tempdecreasevalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(moistureinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Affects how much moisture wind picks up from the ocean. (Earth: 1.0)");

            ImGui::SameLine((float)bottomrightfigures);
            ImGui::Text(moisturevalue.c_str());

            ImGui::Text("   ");

            ImGui::Text(averageinfo.c_str());

            string avetip = "Earth: 14" + degree;

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(avetip.c_str());

            ImGui::SameLine((float)bottomleftfigures);
            ImGui::Text(averagevalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(glacialinfo.c_str());

            string glacialtip = "Areas below this average temperature may show signs of past glaciation. (Earth: 4" + degree + ")";

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(glacialtip.c_str());

            ImGui::SameLine((float)bottomrightfigures);
            ImGui::Text(glacialvalue.c_str());

            ImGui::Text(northpoleinfo.c_str());

            string northtip = "Adjustment to north pole temperature. (Earth: +3" + degree + ")";

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(northtip.c_str());

            ImGui::SameLine((float)bottomleftfigures);
            ImGui::Text(northpolevalue.c_str());

            ImGui::Text(southpoleinfo.c_str());

            string southtip = "Adjustment to south pole temperature. (Earth: -3" + degree + ")";

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(southtip.c_str());

            ImGui::SameLine((float)bottomleftfigures);
            ImGui::Text(southpolevalue.c_str());
                    ImGui::EndTabItem();
                }

                if (ImGui::BeginTabItem("Climate"))
                {
                    ImGui::Text("Percent of land area");
                    ImGui::Text(" ");

                    if (ImGui::BeginTable("ClimateAreas", 3, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY, ImVec2(0.0f, 360.0f)))
                    {
                        ImGui::TableSetupColumn("Code", ImGuiTableColumnFlags_WidthFixed, 55.0f);
                        ImGui::TableSetupColumn("Climate");
                        ImGui::TableSetupColumn("Area %", ImGuiTableColumnFlags_WidthFixed, 80.0f);
                        ImGui::TableHeadersRow();

                        for (int climate = 1; climate <= world->climatenumber() && climate < static_cast<int>(climateareas.size()); climate++)
                        {
                            if (climateareas[climate] == 0)
                                continue;

                            const float percentage = totallandarea > 0 ? (static_cast<float>(climateareas[climate]) * 100.0f) / static_cast<float>(totallandarea) : 0.0f;

                            ImGui::TableNextRow();
                            ImGui::TableSetColumnIndex(0);
                            ImGui::Text("%s", getclimatecode((short)climate).c_str());
                            ImGui::TableSetColumnIndex(1);
                            ImGui::Text("%s", getclimatename((short)climate).c_str());
                            ImGui::TableSetColumnIndex(2);
                            ImGui::Text("%.2f", percentage);
                        }

                        ImGui::EndTable();
                    }

                    ImGui::EndTabItem();
                }

                if (ImGui::BeginTabItem("Tectonics"))
                {
                    const planet& inspectworld = currentglobalworld();

                    ImGui::Text("World metadata");
                    ImGui::Separator();
                    ImGui::Text("time_myr: %.3f", inspectworld.tectonictimemyr());
                    ImGui::Text("delta_time_myr: %.6f", inspectworld.tectonicdeltatimemyr());

                    ImGui::Text("plate_count: %d", inspectworld.tectonicplatecount());
                    ImGui::Text("sea_level_m: %d", inspectworld.tectonicsealevelm());
                    ImGui::Text("boundary_segments: %d", inspectworld.tectonicboundarysegmentcount());
                    ImGui::Text("deforming_regions: %d", inspectworld.tectonicdeformingregioncount());

                    ImGui::Spacing();
                    ImGui::Text("Point inspector");
                    ImGui::Separator();

                    int inspectx = 0;
                    int inspecty = 0;
                    const bool hasinspectionpoint = gettectonicinspectionpoint(currentglobalworld(), *region, screenmode, focused == 1, poix, poiy, inspectx, inspecty);

                    if (hasinspectionpoint == false)
                    {
                        ImGui::TextUnformatted("Select a point on the global or regional map to inspect tectonics.");
                    }
                    else
                    {
                        const planet& inspectworld = currentglobalworld();
                        const float crustagemyr = inspectworld.tectoniccrustagemyr(inspectx, inspecty);
                        const float crustthickness = inspectworld.tectoniccrustthickness(inspectx, inspecty);
                        const CrustClass crustclass = inspectworld.tectoniccrustclass(inspectx, inspecty);
                        const float uplift = inspectworld.tectonicuplifttendency(inspectx, inspecty);
                        const float subsidence = inspectworld.tectonicsubsidencetendency(inspectx, inspecty);
                        const float strain = inspectworld.tectonicaccumulatedstrain(inspectx, inspecty);
                        const BoundaryType boundarytype = inspectworld.tectonicboundarytype(inspectx, inspecty);
                        const int boundarydistance = inspectworld.tectonicboundarydistance(inspectx, inspecty);
                        const int boundarysegmentid = inspectworld.tectonicboundarysegmentid(inspectx, inspecty);
                        const int nearestboundaryid = inspectworld.tectonicnearestboundaryid(inspectx, inspecty);
                        const float boundaryhistory = inspectworld.tectonicboundaryhistory(inspectx, inspecty);
                        const int deformingregionid = inspectworld.tectonicdeformingregionid(inspectx, inspecty);
                        const DeformingRegionType deformingregiontype = inspectworld.tectonicdeformingregiontype(inspectx, inspecty);
                        const float deformationrate = inspectworld.tectonicdeformationrate(inspectx, inspecty);
                        const float deformationvelocityx = inspectworld.tectonicdeformationvelocityx(inspectx, inspecty);
                        const float deformationvelocityy = inspectworld.tectonicdeformationvelocityy(inspectx, inspecty);

                        ImGui::Text("global_x: %d", inspectx);
                        ImGui::Text("global_y: %d", inspecty);
                        ImGui::Text("crust_age_myr: %.3f", crustagemyr);
                        ImGui::Text("crust_thickness: %.3f", crustthickness);
                        ImGui::Text("crust_class: %s", crustclasslabel(crustclass));
                        ImGui::Text("uplift_tendency: %.4f", uplift);
                        ImGui::Text("subsidence_tendency: %.4f", subsidence);
                        ImGui::Text("accumulated_strain: %.4f", strain);
                        ImGui::Text("boundary_type: %s", boundarytypelabel(boundarytype));
                        ImGui::Text("boundary_distance: %d", boundarydistance);
                        ImGui::Text("boundary_segment_id: %d", boundarysegmentid);
                        ImGui::Text("nearest_boundary_id: %d", nearestboundaryid);
                        ImGui::Text("boundary_history (app-derived): %.4f", boundaryhistory);
                        ImGui::Text("deforming_region_id: %d", deformingregionid);
                        ImGui::Text("deforming_region_type: %s", deformingregiontypelabel(deformingregiontype));
                        ImGui::Text("deformation_rate: %.4f", deformationrate);
                        ImGui::Text("deformation_velocity_x: %.4f", deformationvelocityx);
                        ImGui::Text("deformation_velocity_y: %.4f", deformationvelocityy);

                        if (boundarysegmentid > 0 || nearestboundaryid > 0 || deformingregionid > 0)
                        {
                            const TectonicBoundarySegment* boundarysegment = inspectworld.findtectonicboundarysegment(boundarysegmentid);
                            const TectonicBoundarySegment* nearestboundary = inspectworld.findtectonicboundarysegment(nearestboundaryid);
                            const TectonicDeformingRegion* deformingregion = inspectworld.findtectonicdeformingregion(deformingregionid);

                            ImGui::Spacing();
                            ImGui::Text("Linked summaries");
                            ImGui::Separator();

                            drawboundarysegmentsummary("Boundary segment", boundarysegmentid, boundarysegment);

                            if (nearestboundaryid > 0)
                            {
                                if (nearestboundaryid == boundarysegmentid)
                                    ImGui::Text("Nearest boundary: id %d (same as boundary segment)", nearestboundaryid);
                                else
                                    drawboundarysegmentsummary("Nearest boundary", nearestboundaryid, nearestboundary);
                            }

                            drawdeformingregionsummary("Deforming region", deformingregionid, deformingregion);
                        }

                        drawmineraldepositsummary(inspectworld, inspectx, inspecty);
                    }

                    ImGui::EndTabItem();
                }

                ImGui::EndTabBar();
            }

            if (ImGui::Button("Close", ImVec2(120.0f, 0.0f)))
                showworldproperties = 0;

            ImGui::End();
        }
