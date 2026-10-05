// Part of main() in app/main.cpp: included inside its frame loop, and uses its locals.

        if (screenmode == loadfailure)
        {
            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 446, main_viewport->WorkPos.y + 174), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(296, 129), ImGuiCond_FirstUseEver);

            ImGui::Begin("Loading unsuccessful!", NULL, window_flags);

            ImGui::Text("This world file is not compatible with the");
            ImGui::Text("current version of Undiscovered Worlds.");

            ImGui::Text(" ");
            ImGui::Text(" ");

            ImGui::SameLine((float)135);

            if (ImGui::Button("OK"))
            {
                if (brandnew)
                    screenmode = createworldscreen;
                else
                    screenmode = globalmapscreen;
            }

            ImGui::End();
        }
