#include "app/application.hpp"

void Application::updatecompletingimportscreen()
{
    if (screenmode == completingimportscreen)
    {
        if (showdeferredworkwindow(completingimportpass, "Please wait...##completeimport ", "Finishing world...", ImVec2(main_viewport->WorkPos.x + 507, main_viewport->WorkPos.y + 173), ImVec2(173, 68)))
        {
            updatereport("Generating custom world:");
            updatereport("");
            begintimedreporting();

            // Plug in the world settings.

            applyworldpropertycontrols(*world, worldpropertycontrols);
            completeimportedworldgeneration(*world, currentdeltas, socialoptionsfromworldcontrols(), smalllake, largelake, landshape, OKmountains, &importedclimatemaps);

            // Now draw a new map

            invalidateGlobalMaps();

            mapview = relief;

            showGlobalMapView(mapview);

            endtimedreporting();
            updatereport("");
            updatereport("World generation completed.");
            updatereport("");

            focused = 0;

            completingimportpass = 0;
            screenmode = movingtoglobalmapscreen;
            newworld = 1;
        }
    }
}
