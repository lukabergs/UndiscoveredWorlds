#include "app/application.hpp"

const char* crustclasslabel(CrustClass crustclass)
{
    switch (crustclass)
    {
    case CrustClass::oceanic:
        return "oceanic";
    case CrustClass::transitional:
        return "transitional";
    case CrustClass::continental:
        return "continental";
    case CrustClass::none:
    default:
        return "none";
    }
}

const char* boundarytypelabel(BoundaryType boundarytype)
{
    switch (boundarytype)
    {
    case BoundaryType::convergent:
        return "convergent";
    case BoundaryType::divergent:
        return "divergent";
    case BoundaryType::transform:
        return "transform";
    case BoundaryType::passive_margin:
        return "passive_margin";
    case BoundaryType::none:
    default:
        return "none";
    }
}

const char* deformingregiontypelabel(DeformingRegionType regiontype)
{
    switch (regiontype)
    {
    case DeformingRegionType::continental_rift:
        return "continental_rift";
    case DeformingRegionType::diffuse_collision:
        return "diffuse_collision";
    case DeformingRegionType::none:
    default:
        return "none";
    }
}

const char* geologicregimelabel(GeologicRegime regime)
{
    switch (regime)
    {
    case GeologicRegime::convergent_arc:
        return "convergent_arc";
    case GeologicRegime::continent_collision:
        return "continent_collision";
    case GeologicRegime::divergent_rift:
        return "divergent_rift";
    case GeologicRegime::transform:
        return "transform";
    case GeologicRegime::passive_margin:
        return "passive_margin";
    case GeologicRegime::mid_ocean_ridge:
        return "mid_ocean_ridge";
    case GeologicRegime::trench_adjacent:
        return "trench_adjacent";
    case GeologicRegime::stable:
    default:
        return "stable";
    }
}

void drawboundarysegmentsummary(const char* label, int requestedid, const TectonicBoundarySegment* segment)
{
    if (requestedid <= 0)
        return;

    if (segment == nullptr)
    {
        ImGui::Text("%s: id %d (not found)", label, requestedid);
        return;
    }

    ImGui::Text("%s: id %d | type %s | plates %d/%d | cells %d | age %.2f Myr",
        label, segment->id, boundarytypelabel(segment->boundaryType), segment->leftPlateId, segment->rightPlateId,
        segment->cellCount, static_cast<float>(segment->ageMyr));
    ImGui::Text("  centroid(%.1f E, %.1f N) | normal %.3f | shear %.3f | len %.0f km | regime %s",
        segment->centroidLongitudeDeg, segment->centroidLatitudeDeg, segment->averageNormalMotion, segment->averageShearMotion,
        segment->lengthKm, geologicregimelabel(segment->geologicRegime));
}

void drawdeformingregionsummary(const char* label, int requestedid, const TectonicDeformingRegion* region)
{
    if (requestedid <= 0)
        return;

    if (region == nullptr)
    {
        ImGui::Text("%s: id %d (not found)", label, requestedid);
        return;
    }

    ImGui::Text("%s: id %d | type %s | plates %d/%d | cells %d | age %.2f Myr",
        label, region->id, deformingregiontypelabel(region->type), region->primaryPlateId, region->secondaryPlateId,
        region->cellCount, static_cast<float>(region->ageMyr));
    ImGui::Text("  centroid(%.1f E, %.1f N) | rate %.3f | vel(%.3f E, %.3f N) | normal %.3f | shear %.3f",
        region->centroidLongitudeDeg, region->centroidLatitudeDeg, region->averageDeformationRate,
        region->averageInterpolatedVelocityX, region->averageInterpolatedVelocityY,
        region->averageNormalMotion, region->averageShearMotion);
}

// Deposits in the inspected cell, with the history that explains each one.
void drawmineraldepositsummary(const planet& world, int x, int y)
{
    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    int shown = 0;
    for (const mineralsystems::Deposit& deposit : world.mineraldeposits())
    {
        const int depositx = std::clamp(static_cast<int>(deposit.u * columns), 0, columns - 1);
        const int deposity = std::clamp(static_cast<int>(deposit.v * rows), 0, rows - 1);
        if (depositx != x || deposity != y)
            continue;
        if (shown++ == 0)
        {
            ImGui::Spacing();
            ImGui::Text("Mineral deposits");
            ImGui::Separator();
        }
        ImGui::TextWrapped("%s", mineralsystems::describedeposit(deposit, &world.mineraldeposits()).c_str());
    }
}

bool gettectonicinspectionpoint(const planet& world, const region& region, screenmodeenum screenmode, bool focused, int poix, int poiy, int& x, int& y)
{
    if (focused == false)
        return false;

    if (screenmode == globalmapscreen || screenmode == exportareascreen || screenmode == importscreen)
    {
        x = wrap(poix, world.width());
        y = std::clamp(poiy, 0, world.height());
        return true;
    }

    if (screenmode == regionalmapscreen || screenmode == generatingregionscreen)
    {
        x = wrap(region.leftx() + poix / 16, world.width());
        y = std::clamp(region.lefty() + poiy / 16, 0, world.height());
        return true;
    }

    return false;
}

int createrandomseednumber()
{
    int seed = random(0, 9);

    for (int n = 1; n <= 7; n++)
    {
        if (n == 7)
            seed = seed + (random(1, 9) * (int)pow(10, n));
        else
            seed = seed + (random(0, 9) * (int)pow(10, n));
    }

    return seed;
}

void initglobalmapdisplay(planet& world, GlobalMapDisplay& display)
{
    display.textureSize.x = world.width() + 1;
    display.textureSize.y = world.height() + 1;

    for (mapviewenum thismapview : allmapviews)
        getmapimage(display.maps, thismapview).create(display.textureSize.x, display.textureSize.y, sf::Color::Black);

    display.displayTextureSize.x = DISPLAYMAPSIZEX;
    display.displayTextureSize.y = DISPLAYMAPSIZEY;

    for (mapviewenum thismapview : allmapviews)
        getdisplaymapimage(display.maps, thismapview).create(display.displayTextureSize.x, display.displayTextureSize.y, sf::Color::Black);

    updateTextureFromImage(display.texture, getdisplaymapimage(display.maps, relief));

    display.sprite.setTexture(display.texture);
    display.sprite.setPosition(sf::Vector2f(display.x, display.y));

    display.minimap.setTexture(display.texture);
    display.minimap.setScale(0.5f, 0.5f);
    display.minimap.setPosition(sf::Vector2f(display.minimapX, display.minimapY));
}

void inithighlightdisplay(planet& world, HighlightDisplay& display)
{
    display.image.create(display.size, display.size, sf::Color::Transparent);
    display.miniImage.create(display.miniSize, display.miniSize, sf::Color::Transparent);

    drawhighlightobjects(world, display.image, display.size, display.miniImage, display.miniSize);

    updateTextureFromImage(display.texture, display.image);
    display.sprite.setTexture(display.texture);
    display.sprite.setOrigin((float)(display.size / 2), (float)(display.size / 2));

    updateTextureFromImage(display.miniTexture, display.miniImage);
    display.miniSprite.setTexture(display.miniTexture);
    display.miniSprite.setOrigin((float)(display.miniSize / 2), (float)(display.miniSize / 2));
}

void initregionalmapdisplay(RegionalMapDisplay& display)
{
    for (mapviewenum thismapview : allmapviews)
        getmapimage(display.maps, thismapview).create(display.imageWidth, display.imageHeight, sf::Color::Black);

    updateTextureFromImage(display.texture, getmapimage(display.maps, relief));

    display.sprite.setTexture(display.texture);
    display.sprite.setPosition(sf::Vector2f(display.x, display.y));
}

void refreshminihighlightdisplay(HighlightDisplay& display)
{
    updateTextureFromImage(display.miniTexture, display.miniImage);
    display.miniSprite.setTexture(display.miniTexture);
    display.miniSprite.setOrigin((float)display.miniSize / 2.0f, (float)display.miniSize / 2.0f);
}

void refreshhighlightdisplay(planet& world, HighlightDisplay& display)
{
    drawhighlightobjects(world, display.image, display.size, display.miniImage, display.miniSize);

    updateTextureFromImage(display.texture, display.image);
    display.sprite.setTexture(display.texture);

    refreshminihighlightdisplay(display);
}

void resizeglobaldisplayforworld(planet& world, GlobalMapDisplay& globaldisplay, HighlightDisplay& highlightdisplay)
{
    adjustforsize(world, globaldisplay.textureSize, globaldisplay.maps, highlightdisplay.image, highlightdisplay.size, highlightdisplay.miniImage, highlightdisplay.miniSize);
    refreshminihighlightdisplay(highlightdisplay);
}

void configureimguistyle()
{
    ImGuiStyle& style = ImGui::GetStyle();
    const float highlight1 = 0.40f;
    const float highlight2 = 0.40f;
    const float highlight3 = 0.40f;

    style.WindowMinSize = ImVec2(160, 20);
    style.FramePadding = ImVec2(4, 2);
    style.ItemSpacing = ImVec2(6, 4);
    style.ItemInnerSpacing = ImVec2(2, 4);
    style.Alpha = 0.95f;
    style.WindowRounding = 4.0f;
    style.FrameRounding = 2.0f;
    style.IndentSpacing = 6.0f;
    style.ColumnsMinSpacing = 50.0f;
    style.GrabMinSize = 14.0f;
    style.GrabRounding = 16.0f;
    style.ScrollbarSize = 12.0f;
    style.ScrollbarRounding = 16.0f;

    style.Colors[ImGuiCol_Text] = ImVec4(0.86f, 0.93f, 0.89f, 0.78f);
    style.Colors[ImGuiCol_TextDisabled] = ImVec4(0.86f, 0.93f, 0.89f, 0.28f);
    style.Colors[ImGuiCol_WindowBg] = ImVec4(0.13f, 0.14f, 0.17f, 1.00f);
    style.Colors[ImGuiCol_ChildBg] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_PopupBg] = ImVec4(0.20f, 0.22f, 0.27f, 0.9f);
    style.Colors[ImGuiCol_Border] = ImVec4(0.31f, 0.31f, 1.00f, 0.00f);
    style.Colors[ImGuiCol_BorderShadow] = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    style.Colors[ImGuiCol_FrameBg] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_FrameBgHovered] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_FrameBgActive] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_TitleBg] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_TitleBgCollapsed] = ImVec4(0.20f, 0.22f, 0.27f, 0.75f);
    style.Colors[ImGuiCol_TitleBgActive] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_MenuBarBg] = ImVec4(0.20f, 0.22f, 0.27f, 0.47f);
    style.Colors[ImGuiCol_ScrollbarBg] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_ScrollbarGrab] = ImVec4(0.09f, 0.15f, 0.16f, 1.00f);
    style.Colors[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.29f, 0.18f, 0.92f, 0.78f);
    style.Colors[ImGuiCol_ScrollbarGrabActive] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_CheckMark] = ImVec4(0.86f, 0.93f, 0.89f, 0.78f);
    style.Colors[ImGuiCol_SliderGrab] = ImVec4(0.47f, 0.77f, 0.83f, 0.14f);
    style.Colors[ImGuiCol_SliderGrabActive] = ImVec4(0.28f, 0.18f, 0.92f, 1.00f);
    style.Colors[ImGuiCol_Button] = ImVec4(0.47f, 0.77f, 0.83f, 0.14f);
    style.Colors[ImGuiCol_ButtonHovered] = ImVec4(highlight1, highlight2, highlight3, 0.86f);
    style.Colors[ImGuiCol_ButtonActive] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_Header] = ImVec4(highlight1, highlight2, highlight3, 0.76f);
    style.Colors[ImGuiCol_HeaderHovered] = ImVec4(highlight1, highlight2, highlight3, 0.86f);
    style.Colors[ImGuiCol_HeaderActive] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_Separator] = ImVec4(0.14f, 0.16f, 0.19f, 1.00f);
    style.Colors[ImGuiCol_SeparatorHovered] = ImVec4(highlight1, highlight2, highlight3, 0.78f);
    style.Colors[ImGuiCol_SeparatorActive] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_ResizeGrip] = ImVec4(0.47f, 0.77f, 0.83f, 0.04f);
    style.Colors[ImGuiCol_ResizeGripHovered] = ImVec4(highlight1, highlight2, highlight3, 0.78f);
    style.Colors[ImGuiCol_ResizeGripActive] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_PlotLines] = ImVec4(0.86f, 0.93f, 0.89f, 0.63f);
    style.Colors[ImGuiCol_PlotLinesHovered] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_PlotHistogram] = ImVec4(0.86f, 0.93f, 0.89f, 0.63f);
    style.Colors[ImGuiCol_PlotHistogramHovered] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_TextSelectedBg] = ImVec4(highlight1, highlight2, highlight3, 0.43f);
    style.Colors[ImGuiCol_Tab] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_TabHovered] = ImVec4(highlight1, highlight2, highlight3, 0.86f);
    style.Colors[ImGuiCol_TabActive] = ImVec4(highlight1, highlight2, highlight3, 0.86f);
    style.Colors[ImGuiCol_TabUnfocused] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_TabUnfocusedActive] = ImVec4(highlight1, highlight2, highlight3, 0.86f);
}

Application::Application()
{
    string currentversionstring = "Current version: " + to_string(currentversion);

    updatereport(currentversionstring.c_str());

    window.setFramerateLimit(60);
    ImGui::SFML::Init(window);
    ImGui::GetIO().IniFilename = "configs/imgui.ini";

    main_viewport = ImGui::GetMainViewport();

    // Setting up Dear ImGUI style
    // Adapted from Cinder-ImGui by Simon Geilfus - https://github.com/simongeilfus/Cinder-ImGui

    configureimguistyle();

    ImGuiIO& io = ImGui::GetIO();
    font1 = io.Fonts->AddFontDefault();
    io.Fonts->Build();
    IM_ASSERT(font1 != NULL);
    ImGui::SFML::UpdateFontTexture();

    initialiseworld(*world);
    initialisemapcolours(*world);

    width = world->width();
    height = world->height();

    initialiseregion(*world, *region);

    regionalmapimagewidth = region->regwidthend() - region->regwidthbegin() + 1;
    regionalmapimageheight = region->regheightend() - region->regheightbegin() + 1;

    initglobalmapdisplay(*world, globaldisplay);
    inithighlightdisplay(*world, highlightdisplay);
    initregionalmapdisplay(regionaldisplay);

    arearectangle.setOutlineThickness(2);
    arearectangle.setFillColor(sf::Color::Transparent);

    createlandshapetemplates(landshape);
    createchainlandtemplates(chainland);
    createsmalllaketemplates(smalllake);
    createlargelaketemplates(largelake);
    createislandtemplates(island);
    createsmudgetemplates(smudge);
    createsmallsmudgetemplates(smallsmudge);
    loadpeaktemplates(peaks);
    createriftblob(riftblob, riftblobsize / 2);

    fast_srand((long)time(0));

    for (int n = 0; n <= ((MAXCRATERRADIUS * MAXCRATERRADIUS + MAXCRATERRADIUS) * 24); n++)
        squareroot[n] = (int)sqrt(n);

    resetmapcache(globalmaps);
    resetmapcache(regionalmaps);

    worldpropertycontrols = makeworldpropertycontrols(*world);
    appearance = makeappearancesettings(*world);
}

void Application::openFileDialog(const char* filter)
{
    IGFD::FileDialogConfig dialogConfig;
    const AppEnvironmentConfig& appenv = getappenvironment();
    string dialogpath = appenv.defaultWorldDirectory.string();

    if (strcmp(filter, ".uws") == 0)
        dialogpath = appenv.defaultAppearanceDirectory.string();
    else if (strstr(filter, ".png") != nullptr || strstr(filter, ".tif") != nullptr)
        dialogpath = appenv.defaultImageDirectory.string();

    if (dialogpath.empty())
        dialogpath = ".";

    dialogConfig.path = dialogpath;
    ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey", "Choose File", filter, dialogConfig);
}

// Runs one registered process on the world, as a recipe step does (its own seeded random stream).
bool Application::runprocess(const char* id)
{
    const WorldProcessDefinition* process = findworldprocess(id);
    if (process == nullptr)
        return false;

    GenerationScratch scratch;
    WorldStageParameters parameters;
    GenerationExecutionContext context;
    context.dodeltas = currentdeltas;
    context.socialoptions = socialoptionsfromworldcontrols();
    context.landshape = landshape;
    context.chainland = chainland;
    context.smalllake = smalllake;
    context.largelake = largelake;
    context.okmountains = &OKmountains;

    std::string error;
    const bool done = executegenerationstage(process->stage, *world, scratch, parameters, context, &error);
    if (!done)
        updatereport(error);
    return done;
}

void Application::invalidateGlobalMaps()
{
    resetmapcache(globalmaps);
}

void Application::invalidateRegionalMaps()
{
    resetmapcache(regionalmaps);
}

planet& Application::currentglobalworld()
{
    if (screenmode == importscreen && generationworkbench.active && generationworkbench.ui.previewEnabled && generationworkbench.previewAvailable)
        return generationworkbench.previewWorld;

    return *world;
}

sf::Color Application::blendcolour(const sf::Color& low, const sf::Color& high, float factor)
{
    factor = std::clamp(factor, 0.0f, 1.0f);

    auto lerpchannel = [&](sf::Uint8 left, sf::Uint8 right)
    {
        return static_cast<sf::Uint8>(std::clamp(static_cast<int>(std::round(static_cast<float>(left) + (static_cast<float>(right) - static_cast<float>(left)) * factor)), 0, 255));
    };

    return sf::Color(lerpchannel(low.r, high.r), lerpchannel(low.g, high.g), lerpchannel(low.b, high.b));
}

sf::Color Application::tosfcolour(const ImVec4& colour, float alpha)
{
    return sf::Color(
        static_cast<sf::Uint8>(std::clamp(static_cast<int>(std::round(colour.x * 255.0f)), 0, 255)),
        static_cast<sf::Uint8>(std::clamp(static_cast<int>(std::round(colour.y * 255.0f)), 0, 255)),
        static_cast<sf::Uint8>(std::clamp(static_cast<int>(std::round(colour.z * 255.0f)), 0, 255)),
        static_cast<sf::Uint8>(std::clamp(static_cast<int>(std::round(alpha * 255.0f)), 0, 255)));
}

sf::Color Application::terrainrendercolour(const planet& renderworld, int worldx, int worldy)
{
    const int altitude = renderworld.map(worldx, worldy) - renderworld.sealevel();
    return tosfcolour(samplegradientcolour(generationworkbench.ui.terrainMapGradient, altitude));
}

sf::Color Application::tectonicboundaryoverlaycolour(const planet& renderworld, int worldx, int worldy, bool& visible)
{
    visible = false;
    const BoundaryType boundarytype = renderworld.tectonicboundarytype(worldx, worldy);
    const DeformingRegionType regiontype = renderworld.tectonicdeformingregiontype(worldx, worldy);

    if (boundarytype == BoundaryType::none && regiontype == DeformingRegionType::none)
        return sf::Color::Transparent;

    visible = true;
    sf::Color target(194, 177, 127);

    switch (boundarytype)
    {
    case BoundaryType::convergent:
        target = sf::Color(186, 86, 66);
        break;
    case BoundaryType::divergent:
        target = sf::Color(216, 166, 78);
        break;
    case BoundaryType::transform:
        target = sf::Color(57, 133, 142);
        break;
    case BoundaryType::passive_margin:
        target = sf::Color(194, 177, 127);
        break;
    case BoundaryType::none:
    default:
        break;
    }

    if (regiontype == DeformingRegionType::continental_rift)
        target = blendcolour(target, sf::Color(235, 193, 88), 0.45f);
    else if (regiontype == DeformingRegionType::diffuse_collision)
        target = blendcolour(target, sf::Color(160, 92, 72), 0.55f);

    return target;
}

void Application::styleworkbenchgloballayer(mapviewenum view, planet& renderworld)
{
    if (screenmode != importscreen || generationworkbench.active == false || workbenchfinished(generationworkbench))
        return;

    const GenerationStageId stageid = getcurrentgenerationstage(generationworkbench).id;
    const GenerationScratch& visiblescratch = (generationworkbench.previewAvailable && generationworkbench.ui.previewEnabled)
        ? generationworkbench.previewScratch
        : generationworkbench.committedScratch;
    const bool terrainrenderactive = isterraingenerationstage(stageid);
    const bool showplateboundaries = generationworkbench.ui.showPlateBoundaries;
    const bool showfeatureoverlay = generationworkbench.ui.showFeatureOverlay && visiblescratch.stageFeatureCellCount > 0;
    const bool inlandseaselection = stageid == GenerationStageId::inland_seas;

    if (!terrainrenderactive && !showplateboundaries && !showfeatureoverlay && !inlandseaselection)
        return;

    auto applystyle = [&](sf::Image& image)
    {
        const int imagewidth = static_cast<int>(image.getSize().x);
        const int imageheight = static_cast<int>(image.getSize().y);
        const int worldwidth = renderworld.width() + 1;
        const int worldheight = renderworld.height() + 1;

        for (int imagey = 0; imagey < imageheight; imagey++)
        {
            const int starty = (imagey * worldheight) / imageheight;
            const int endy = std::max(starty, ((imagey + 1) * worldheight - 1) / imageheight);

            for (int imagex = 0; imagex < imagewidth; imagex++)
            {
                const int startx = (imagex * worldwidth) / imagewidth;
                const int endx = std::max(startx, ((imagex + 1) * worldwidth - 1) / imagewidth);
                const int samplex = std::clamp((startx + endx) / 2, 0, renderworld.width());
                const int sampley = std::clamp((starty + endy) / 2, 0, renderworld.height());

                if (terrainrenderactive)
                    image.setPixel(imagex, imagey, terrainrendercolour(renderworld, samplex, sampley));

                if (inlandseaselection)
                {
                    bool overlayvisible = false;
                    sf::Color overlaycolour = sf::Color::Transparent;

                    for (int worldy = starty; worldy <= endy && !overlayvisible; worldy++)
                    {
                        for (int worldx = startx; worldx <= endx && !overlayvisible; worldx++)
                        {
                            const int wrappedx = wrap(worldx, renderworld.width());
                            const int clampedy = std::clamp(worldy, 0, renderworld.height());
                            const mesh::CellId cell = renderworld.grid().cellof(wrappedx, clampedy);
                            const int componentid = cell < visiblescratch.inlandSeaComponentIds.size() ? visiblescratch.inlandSeaComponentIds[cell] : 0;
                            if (componentid <= 0)
                                continue;

                            for (const GenerationWaterComponent& component : visiblescratch.inlandSeaComponents)
                            {
                                if (component.id != componentid)
                                    continue;

                                overlayvisible = true;
                                switch (component.policy)
                                {
                                case GenerationComponentPolicy::keep:
                                    overlaycolour = sf::Color(255, 0, 0);
                                    break;
                                case GenerationComponentPolicy::drain:
                                    overlaycolour = sf::Color(0, 255, 0);
                                    break;
                                case GenerationComponentPolicy::fill:
                                    overlaycolour = sf::Color(0, 0, 255);
                                    break;
                                }
                                break;
                            }
                        }
                    }

                    if (overlayvisible)
                        image.setPixel(imagex, imagey, overlaycolour);
                }

                if (showfeatureoverlay)
                {
                    bool overlayvisible = false;

                    for (int worldy = starty; worldy <= endy && !overlayvisible; worldy++)
                    {
                        for (int worldx = startx; worldx <= endx && !overlayvisible; worldx++)
                        {
                            const int wrappedx = wrap(worldx, renderworld.width());
                            const int clampedy = std::clamp(worldy, 0, renderworld.height());
                            const mesh::CellId cell = renderworld.grid().cellof(wrappedx, clampedy);
                            if (cell < visiblescratch.stageFeatureMask.size() && visiblescratch.stageFeatureMask[cell] != 0)
                            {
                                overlayvisible = true;
                                break;
                            }
                        }
                    }

                    if (overlayvisible)
                        image.setPixel(imagex, imagey, blendcolour(image.getPixel(imagex, imagey), sf::Color(255, 0, 0), generationworkbench.ui.featureOverlayOpacity));
                }

                if (showplateboundaries)
                {
                    bool overlayvisible = false;
                    sf::Color overlaycolour = sf::Color::Transparent;

                    for (int worldy = starty; worldy <= endy && overlayvisible == false; worldy++)
                    {
                        for (int worldx = startx; worldx <= endx; worldx++)
                        {
                            overlaycolour = tectonicboundaryoverlaycolour(renderworld, wrap(worldx, renderworld.width()), std::clamp(worldy, 0, renderworld.height()), overlayvisible);
                            if (overlayvisible)
                                break;
                        }
                    }

                    if (overlayvisible)
                        image.setPixel(imagex, imagey, blendcolour(image.getPixel(imagex, imagey), overlaycolour, 0.72f));
                }
            }
        }
    };

    maplayer& layer = getmaplayer(globalmaps, view);
    applystyle(layer.image);
    applystyle(layer.displayimage);
}

void Application::showGlobalMapView(mapviewenum view, bool updateMiniMap)
{
    mapview = view;
    planet& renderworld = currentglobalworld();
    withworkbenchglobalrenderoverrides(renderworld, [&]()
    {
        drawglobalmapimage(view, renderworld, globalmaps);
        styleworkbenchgloballayer(view, renderworld);
        updateTextureFromImage(globalmaptexture, getdisplaymapimage(globalmaps, view));
        globalmap.setTexture(globalmaptexture);

        if (updateMiniMap)
            minimap.setTexture(globalmaptexture);
    });
}

void Application::redrawGlobalRelief(bool updateMiniMap)
{
    showGlobalMapView(relief, updateMiniMap);
}

void Application::showRegionalMapView(mapviewenum view)
{
    mapview = view;
    applyregionalmapview(view, *world, *region, regionalmaps, regionalmaptexture, regionalmap);
}

void Application::redrawRegionalRelief()
{
    showRegionalMapView(relief);
}

void Application::rebuildAllGlobalMaps()
{
    planet& renderworld = currentglobalworld();
    withworkbenchglobalrenderoverrides(renderworld, [&]()
    {
        drawallglobalmapimages(renderworld, globalmaps);
    });
}

void Application::rebuildAllRegionalMaps()
{
    drawallregionalmapimages(*world, *region, regionalmaps);
}

SocialGenerationOptions Application::socialoptionsfromdebug()
{
    return makesocialoptions(
        worldgenerationdebug.options.socialEnabled,
        worldgenerationdebug.options.socialMode,
        worldgenerationdebug.options.usePrehistory,
        worldgenerationdebug.options.historyYears);
}

void Application::applytectonicsoptionsfromdebug(GenerationWorkbenchUiState& ui)
{
    ui.tectonicDurationMyr = worldgenerationdebug.options.plateTectonicsDurationMyr;
    ui.tectonicPlateCount = worldgenerationdebug.options.plateTectonicsPlateCount;
    ui.tectonicDeltaTimeMyr = worldgenerationdebug.options.plateTectonicsDeltaTimeMyr;
}

SocialGenerationOptions Application::socialoptionsfromworldcontrols()
{
    return makesocialoptions(
        currentsocialenabled,
        currentsocialmode,
        currentuseprehistory,
        currenthistoryyears);
}

void Application::startproceduralworkbenchsession()
{
    initialiseworld(*world);
    world->clear();
    resetworkbenchsession(generationworkbench);
    clearimportedclimatemaps(importedclimatemaps);
    mapimportsettings = MapImportSettings{};

    changeworldproperties(*world);
    applyworlddimensions(*world, worldpropertycontrols);
    world->setseed(seedentry);
    syncworldpropertycontrols(*world, worldpropertycontrols);
    currentdeltas = true;

    resizeglobaldisplayforworld(*world, globaldisplay, highlightdisplay);
    fast_srand(world->seed());

    initializeproceduralworkbenchsession(generationworkbench, *world, worldgenerationdebug.options.plateTectonicsDurationMyr, worldgenerationdebug.options.plateTectonicsPlateCount);
    applytectonicsoptionsfromdebug(generationworkbench.ui);

    screenmode = importscreen;
    brandnew = 0;

    GenerationExecutionContext context;
    context.dodeltas = currentdeltas;
    context.socialoptions = socialoptionsfromworldcontrols();
    context.landshape = landshape;
    context.chainland = chainland;
    context.smalllake = smalllake;
    context.largelake = largelake;
    context.okmountains = &OKmountains;

    runworkbenchaction([&]()
    {
        previewcurrentgenerationstage(generationworkbench, *world, context, nullptr);
    });

    invalidateGlobalMaps();
    showGlobalMapView(preferredworkbenchmapview(getcurrentgenerationstage(generationworkbench).id));
}

bool Application::showdeferredworkwindow(short& pass, const char* windowTitle, const char* message, ImVec2 position, ImVec2 size)
{
    if (pass >= 10)
        return true;

    ImGui::SetNextWindowPos(position, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(size, ImGuiCond_FirstUseEver);
    ImGui::Begin(windowTitle, NULL, window_flags);
    ImGui::Text(message);
    ImGui::End();

    pass++;
    return false;
}

void Application::run()
{
    sf::Clock deltaClock;
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            ImGui::SFML::ProcessEvent(window, event);

            if (event.type == sf::Event::Closed)
                window.close();
        }

        ImGuiIO& io = ImGui::GetIO();

        ImGui::SFML::Update(window, deltaClock.restart());

        if (screenmode == movingtoglobalmapscreen)
            screenmode = globalmapscreen;

        ImGui::PushFont(font1);

        //ImGui::ShowDemoWindow();

        // First, draw the GUI.

        updatecreateworldscreen();
        updatecreatingworldscreen();
        updateglobalmapscreen();
        updateregionalmapscreen();
        updateexportareascreen();
        updateimportscreen();

        // These screens all display a "Please wait" message ten times (for some reason doing it once or twice doesn't actually display it) and then do something time-consuming.
        updatecompletingimportscreen();
        updateloadingworldscreen();
        updatesavingworldscreen();
        updateexportingareascreen();
        updategeneratingregionscreen();
        updateloadfailurescreen();
        updatesettingsloadfailurescreen();

        window.clear();

        updateappearancepanel();
        updateworldpropertiespanel();
        drawchartwindows();
        drawmaps();
        handlefiledialogs();
        drawdialogwindows();

        ImGui::PopFont();

        ImGui::SFML::Render(window);
        window.display();

        applyappearancechanges();
    }

    ImGui::SFML::Shutdown();
}

void Application::drawchartwindows()
{
    // World edit properties screen, if being shown

    drawworldeditpropertieswindow(main_viewport, window_flags, worldpropertycontrols, degree, showworldeditproperties);

    // Global temperature chart, if being shown

    if (showglobaltemperaturechart)
        drawglobaltemperaturechartwindow(main_viewport, window_flags, *world, poix, poiy);

    // Global precipitation chart, if being shown

    if (showglobalrainfallchart)
        drawglobalrainfallchartwindow(main_viewport, window_flags, *world, poix, poiy);

    // Regional temperature chart, if being shown

    if (showregionaltemperaturechart)
        drawregionaltemperaturechartwindow(main_viewport, window_flags, *world, *region, poix, poiy);

    // Regional precipitation chart, if being shown

    if (showregionalrainfallchart)
        drawregionalrainfallchartwindow(main_viewport, window_flags, *region, poix, poiy);
}

void Application::drawmaps()
{
    // Now draw the graphical elements.

    if (screenmode == globalmapscreen || screenmode == exportareascreen || screenmode == importscreen)
        window.draw(globalmap);

    if (screenmode == exportareascreen) // Area selection rectangle.
    {
        if (areaswx != -1)
        {
            float mult = ((float)world->width() + 1.0f) / (float)DISPLAYMAPSIZEX;

            float sizex = ((float)areanex - (float)areanwx) / mult;
            float sizey = ((float)areaswy - (float)areanwy) / mult;

            arearectangle.setSize(sf::Vector2f(sizex, sizey));
            arearectangle.setOutlineColor(sf::Color(world->highlight1(), world->highlight2(), world->highlight3()));

            sf::Vector2f position;

            position.x = (float)areanwx / mult + (float)globalmapxpos;
            position.y = (float)areanwy / mult + (float)globalmapypos;

            arearectangle.setPosition(position);

            window.draw(arearectangle);
        }
    }

    if (screenmode == regionalmapscreen || screenmode == generatingregionscreen)
    {
        window.draw(regionalmap);

        window.draw(minimap);

        float mult = (float)DISPLAYMAPSIZEX / ((float)world->width() + 1.0f);
        mult = mult * 0.5f;

        float posx = (float)minimapxpos + (float)region->centrex() * mult;
        float posy = (float)minimapypos + (float)region->centrey() * mult;

        minihighlight.setPosition(sf::Vector2f(posx, posy));

        if (screenmode != generatingregionscreen)
        {
            window.draw(minihighlight);
        }
    }

    if (screenmode == globalmapscreen || screenmode == regionalmapscreen)
    {
        if (focused == 1)
            window.draw(highlight);
    }
}

void Application::handlefiledialogs()
{
    // Now draw file dialogues if needed

    handlefiledialog(loadingworld, filepathname, filepath, [&](const std::string&, const std::string&)
    {
        resetworkbenchsession(generationworkbench);
        screenmode = loadingworldscreen;
    });

    handlefiledialog(savingworld, filepathname, filepath, [&](const std::string&, const std::string&)
    {
        screenmode = savingworldscreen;
    });

    handlefiledialog(loadingsettings, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
    {
        bool found = loadsettings(*world, selectedpath);

        if (found == 1)
        {
            syncappearancesettings(*world, appearance);
            refreshhighlightdisplay(*world, highlightdisplay);

            colourschanged = 1;
        }
        else
        {
            oldscreenmode = screenmode;
            screenmode = settingsloadfailure;
        }
    });

    handlefiledialog(savingsettings, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
    {
        savesettings(*world, selectedpath);
    });

    handlefiledialog(exportingworldmaps, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
    {
        rebuildAllGlobalMaps();

        std::string exportpath = selectedpath;
        exportpath.resize(exportpath.size() - 4);

        for (const mapviewdefinition& definition : allmapviewdefinitions)
            getmapimage(globalmaps, definition.view).saveToFile(exportpath + " " + definition.exportstem + ".png");
    });

    handlefiledialog(exportingregionalmaps, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
    {
        rebuildAllRegionalMaps();

        std::string exportpath = selectedpath;
        exportpath.resize(exportpath.size() - 4);

        for (const mapviewdefinition& definition : allmapviewdefinitions)
            getmapimage(regionalmaps, definition.view).saveToFile(exportpath + " " + definition.exportstem + ".png");
    });

    handlefiledialog(exportingareamaps, filepathname, filepath, [&](const std::string&, const std::string&)
    {
        screenmode = exportingareascreen;
    });

    // These sections are for loading in maps in the import screen.

    handlefiledialog(importinglandmap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
    {
        if (importlandheightmap(*world, selectedpath, &mapimportsettings))
        {
            resetworkbenchsession(generationworkbench);
            invalidateGlobalMaps();
            showGlobalMapView(relief);
        }
    });

    handlefiledialog(importingseamap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
    {
        if (importseadepthmap(*world, selectedpath, &mapimportsettings))
        {
            resetworkbenchsession(generationworkbench);
            invalidateGlobalMaps();
            showGlobalMapView(relief);
        }
    });

    handlefiledialog(importingmountainsmap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
    {
        if (importmountainmap(*world, selectedpath, OKmountains, &mapimportsettings))
        {
            resetworkbenchsession(generationworkbench);
            invalidateGlobalMaps();
            showGlobalMapView(relief);
        }
    });

    handlefiledialog(importingvolcanoesmap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
    {
        if (importvolcanomap(*world, selectedpath, &mapimportsettings))
        {
            resetworkbenchsession(generationworkbench);
            invalidateGlobalMaps();
            showGlobalMapView(relief);
        }
    });

    handlefiledialog(importingtemperaturemap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
    {
        if (importtemperaturemap(*world, selectedpath, importedclimatemaps, &mapimportsettings))
        {
            resetworkbenchsession(generationworkbench);
            invalidateGlobalMaps();
            showGlobalMapView(temperature);
        }
    });

    handlefiledialog(importingprecipitationmap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
    {
        if (importprecipitationmap(*world, selectedpath, importedclimatemaps, &mapimportsettings))
        {
            resetworkbenchsession(generationworkbench);
            invalidateGlobalMaps();
            showGlobalMapView(precipitation);
        }
    });

    handlefiledialog(importinggradientstrip, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
    {
        mapimportsettings.gradientStripPath = selectedpath;
    });
}

void Application::drawdialogwindows()
{
    // Set size window, for custom worlds.

    if (drawcustomworldsizewindow(main_viewport, window_flags, showsetsize, brandnew, worldpropertycontrols))
    {
        initialiseworld(*world);
        world->clear();
        resetworkbenchsession(generationworkbench);
        clearimportedclimatemaps(importedclimatemaps);
        mapimportsettings = MapImportSettings{};

        world->setsize(currentsize);
        applyworlddimensions(*world, worldpropertycontrols);

        resizeglobaldisplayforworld(*world, globaldisplay, highlightdisplay);
        world->setgravity(getdefaultgravityforsize(currentsize));

        syncworldpropertycontrols(*world, worldpropertycontrols);
        currentdeltas = true;

        int width = world->width();
        int height = world->height();
        int val = world->sealevel() - 5000;

        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
                world->setnom(i, j, val);
        }

        mapview = relief;

        resetmapcache(globalmaps);

        showGlobalMapView(mapview);

        world->setseed(0 - createrandomseednumber());

        screenmode = importscreen;
    }

    // Window for creating custom worlds.

    if (drawterrainchooserwindow(main_viewport, window_flags, showterrainchooser))
    {
        resetworkbenchsession(generationworkbench);
        initializeproceduralworkbenchsession(generationworkbench, *world, worldgenerationdebug.options.plateTectonicsDurationMyr, worldgenerationdebug.options.plateTectonicsPlateCount);
        applytectonicsoptionsfromdebug(generationworkbench.ui);
        previewcurrentgenerationstage(generationworkbench, *world, GenerationExecutionContext{}, nullptr);
        invalidateGlobalMaps();
        showGlobalMapView(preferredworkbenchmapview(getcurrentgenerationstage(generationworkbench).id));
    }

    // Warning window for over-large area maps.

    if (drawareawarningwindow(main_viewport, window_flags, showareawarning))
    {
        openFileDialog(".png");
        exportingareamaps = 1;
    }

    // Window for displaying information about the program.

    drawaboutwindow(main_viewport, window_flags, showabout, currentversion);
}

void Application::applyappearancechanges()
{
    // Now update the colours if necessary.

    const AppearanceChangeFlags appearanceChanges = getappearancechanges(*world, appearance);
    const bool shouldApplyAppearance = showcolouroptions == 1 && (colourschanged == 1 || appearanceChanges.any());
    const bool shouldRedrawAppearanceMaps = colourschanged == 1 || appearanceChanges.mapAppearanceChanged;

    if (shouldApplyAppearance)
    {
        applyappearancesettings(*world, appearance);

        if (appearanceChanges.highlightChanged)
            refreshhighlightdisplay(*world, highlightdisplay);

        if (shouldRedrawAppearanceMaps)
        {
            invalidateGlobalMaps();
            showGlobalMapView(mapview, true);

            if (screenmode == regionalmapscreen)
            {
                invalidateRegionalMaps();
                showRegionalMapView(mapview);
            }
        }

        colourschanged = 0;
    }
}
