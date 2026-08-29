#pragma once

#include "imgui.h"

struct CreateWorldScreenActions
{
    bool openLoadDialog = false;
    bool openCustomWorld = false;
    bool openPlateTectonicsMenu = false;
};

CreateWorldScreenActions drawcreateworldscreen(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, float currentversion, float latestversion);
