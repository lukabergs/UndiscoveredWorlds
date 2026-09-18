#pragma once

#include <string>

#include "ImGuiFileDialog.h"

struct FileDialogState
{
    bool loadingWorld = false;
    bool savingWorld = false;
    bool loadingSettings = false;
    bool savingSettings = false;
    bool exportingWorldMaps = false;
    bool exportingRegionalMaps = false;
    bool exportingAreaMaps = false;
    bool importingLandMap = false;
    bool importingSeaMap = false;
    bool importingMountainsMap = false;
    bool importingVolcanoesMap = false;
    bool importingTemperatureMap = false;
    bool importingPrecipitationMap = false;
    bool importingGradientStrip = false;
    std::string filepathname;
    std::string filepath;
};

template <typename Callback>
void handlefiledialog(bool& active, std::string& filepathname, std::string& filepath, Callback&& onAccept)
{
    if (!active || !ImGuiFileDialog::Instance()->Display("ChooseFileDlgKey"))
        return;

    filepathname.clear();
    filepath.clear();

    if (ImGuiFileDialog::Instance()->IsOk())
    {
        filepathname = ImGuiFileDialog::Instance()->GetFilePathName();
        filepath = ImGuiFileDialog::Instance()->GetCurrentPath();

        if (!filepathname.empty())
            onAccept(filepathname, filepath);
    }

    ImGuiFileDialog::Instance()->Close();
    active = false;
}
