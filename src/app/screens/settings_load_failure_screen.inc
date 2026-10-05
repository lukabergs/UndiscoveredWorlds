// Part of main() in app/main.cpp: included inside its frame loop, and uses its locals.

        if (screenmode == settingsloadfailure)
        {
            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 446, main_viewport->WorkPos.y + 174), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(310, 129), ImGuiCond_FirstUseEver);

            ImGui::Begin("Loading unsuccessful!##settings ", NULL, window_flags);

            ImGui::Text("This settings file is not compatible with the");
            ImGui::Text("current version of Undiscovered Worlds.");

            ImGui::Text(" ");
            ImGui::Text(" ");

            ImGui::SameLine((float)135);

            if (ImGui::Button("OK"))
                screenmode = oldscreenmode;

            ImGui::End();
        }
