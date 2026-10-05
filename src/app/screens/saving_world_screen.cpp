#include "app/application.hpp"

void Application::updatesavingworldscreen()
{
    if (screenmode == savingworldscreen)
    {
        if (showdeferredworkwindow(savingworldpass, "Please wait...##saveworld ", "Saving world...", ImVec2(main_viewport->WorkPos.x + 507, main_viewport->WorkPos.y + 173), ImVec2(173, 68)))
        {
            world->saveworld(filepathname);

            filepathname = "";
            filepath = "";

            savingworld = 0;

            savingworldpass = 0;

            screenmode = movingtoglobalmapscreen;
        }
    }
}
