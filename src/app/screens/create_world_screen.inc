// Part of main() in app/main.cpp: included inside its frame loop, and uses its locals.

        if (screenmode == createworldscreen)
        {
            showglobaltemperaturechart = 0;
            showglobalrainfallchart = 0;
            showregionaltemperaturechart = 0;
            showregionalrainfallchart = 0;
            showworldproperties = 0;
            showworldeditproperties = 0;
            newworld = 0;

            bool openedplatetectonicsmenu = false;

            if (!showworldgenerationoptions)
            {
                const CreateWorldScreenActions actions = drawcreateworldscreen(main_viewport, window_flags);

                if (actions.openLoadDialog)
                {
                    openFileDialog(".uww");
                    loadingworld = 1;
                }

                if (actions.openCustomWorld)
                    showsetsize = 1;

                if (actions.openPlateTectonicsMenu)
                {
                    showworldgenerationoptions = 1;
                    openedplatetectonicsmenu = true;
                }
            }

            if (showworldgenerationoptions && openedplatetectonicsmenu == false &&
                drawworldgenerationoptionswindow(main_viewport, window_flags, showworldgenerationoptions, worldpropertycontrols))
            {
                startproceduralworkbenchsession();
            }
        }
