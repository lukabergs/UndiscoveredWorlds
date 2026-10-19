#include "app/application.hpp"
#include "pipeline/export_routing.hpp"

void Application::updatesavingworldscreen()
{
    if (screenmode == savingworldscreen)
    {
        if (showdeferredworkwindow(savingworldpass, "Please wait...##saveworld ", "Saving world...", ImVec2(main_viewport->WorkPos.x + 507, main_viewport->WorkPos.y + 173), ImVec2(173, 68)))
        {
            saveworldfile(*world, filepathname);   // a cube world saves its lat-lon export

            filepathname = "";
            filepath = "";

            savingworld = 0;

            savingworldpass = 0;

            screenmode = movingtoglobalmapscreen;
        }
    }
}
