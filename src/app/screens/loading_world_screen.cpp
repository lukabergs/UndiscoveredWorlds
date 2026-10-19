#include "app/application.hpp"
#include "pipeline/export_routing.hpp"

// A .uww save, or a run bundle chosen by its manifest.json (a cube bundle gives its lat-lon export). A bundle carries the world's fields:
// appearance settings are the defaults, and the tables (settlements, polities, routes, deposits) are not loaded.
bool Application::loadworldfile(const std::string& filename)
{
    const std::filesystem::path path(filename);
    if (path.extension() != ".json")
        return world->loadworld(filename);

    try
    {
        initialiseworld(*world);
        initialisemapcolours(*world);
        world->clear();
        loadplanetbundle(*world, path.parent_path());
        // The desktop runs lat-lon worlds: a cube bundle's world is replaced by its export.
        if (world->grid().kind() == mesh::GridKind::cube)
            *world = *exportworld(*world);
        return true;
    }
    catch (const std::exception& error)
    {
        updatereport(error.what());
        return false;
    }
}

void Application::updateloadingworldscreen()
{
    if (screenmode == loadingworldscreen)
    {
        if (showdeferredworkwindow(loadingworldpass, "Please wait...##loadworld ", "Loading world...", ImVec2(main_viewport->WorkPos.x + 507, main_viewport->WorkPos.y + 173), ImVec2(173, 68)))
        {
            bool success = loadworldfile(filepathname);

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
}
