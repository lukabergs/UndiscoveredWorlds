#include <string>

#include "ui_panels.hpp"
#include "functions.hpp"

using namespace std;

namespace
{
    struct infopanellayout
    {
        float pos2 = 100.0f;
        float pos3 = 220.0f;
        float pos4 = 350.0f;
        float pos5 = 570.0f;
        float pos6 = 690.0f;
        float pos7 = 850.0f;
        float aboutx = 993.0f;
        float abouty = 122.0f;
    };

    string formatcoordinate(int degrees, int minutes, int seconds, bool negative, const string& degree)
    {
        string text = to_string(degrees) + degree + " " + to_string(minutes) + "' " + to_string(seconds);

        if (negative)
            text = "-" + text;

        return text;
    }

    string formatseaicetext(int seaice)
    {
        if (seaice == 2)
            return "permanent";

        if (seaice == 1)
            return "seasonal";

        return "none";
    }

    string getwindtext(int wind)
    {
        if (wind == 0 || wind > 50)
            return "none";

        if (wind > 0)
            return "westerly";

        return "easterly";
    }

    string getglobalspecialtext(planet& world, int x, int y, bool sea)
    {
        string specialtext;

        switch (world.special(x, y))
        {
        case 110:
            specialtext = "Salt pan";
            break;

        case 120:
            specialtext = "Dunes";
            break;

        case 130:
            specialtext = "Wetlands";
            break;

        case 131:
            specialtext = "Brackish wetlands";
            break;

        case 132:
            specialtext = "Salt wetlands";
            break;
        }

        if (world.volcano(x, y) > 0)
            specialtext = sea ? "Submarine volcano" : "Volcano";

        return specialtext;
    }

    string getregionalspecialtext(region& region, int x, int y, bool sea, bool river)
    {
        string specialtext;

        if (region.barrierisland(x, y) == 0 && (region.rivervalley(x, y) == 1 || river))
            specialtext = "River valley";

        switch (region.special(x, y))
        {
        case 110:
            specialtext = "Salt pan";
            break;

        case 120:
            specialtext = "Dunes";
            break;

        case 130:
            specialtext = "Wetlands";
            break;

        case 131:
            specialtext = "Brackish wetlands";
            break;

        case 132:
            specialtext = "Salt wetlands";
            break;
        }

        if (region.mud(x, y))
            specialtext = "Mud flats";

        if (region.sand(x, y))
            specialtext = "Sandy beach";

        if (region.shingle(x, y))
            specialtext = "Shingle beach";

        if (region.sand(x, y) && region.shingle(x, y))
            specialtext = "Sandy/shingle beach";

        if (region.mud(x, y) && region.sand(x, y))
            specialtext = "Muddy/sandy beach";

        if (region.mud(x, y) && region.shingle(x, y))
            specialtext = "Muddy/shingle beach";

        if (region.mud(x, y) && region.sand(x, y) && region.shingle(x, y))
            specialtext = "Muddy/sandy/shingle beach";

        if (region.volcano(x, y))
            specialtext = sea ? "Submarine volcano" : "Volcano";

        if (region.barrierisland(x, y))
        {
            if (region.sand(x, y))
                specialtext = "Sandbar";
            else if (specialtext == "")
                specialtext = "Barrier island";
            else
                specialtext = "Barrier island. " + specialtext;
        }

        if (region.mangrove(x, y))
            specialtext = specialtext + ". Mangrove";

        return specialtext;
    }

    void drawinfowindowfooter(const infopanellayout& layout, bool& showabout)
    {
        ImGui::SetCursorPosX(layout.aboutx);
        ImGui::SetCursorPosY(layout.abouty);

        if (ImGui::Button("?"))
            toggle(showabout);
    }
}

void drawglobalinfowindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, bool focused, int poix, int poiy, const string& degree, const string& cube, bool& showtemperaturechart, bool& showrainfallchart, bool newworld, bool& showabout)
{
    const infopanellayout layout;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 180, main_viewport->WorkPos.y + 539), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(1023, 155), ImGuiCond_FirstUseEver);
    ImGui::Begin("Information", NULL, window_flags);

    if (focused)
    {
        const int sealevel = world.sealevel();
        const bool lake = world.truelake(poix, poiy) != 0;
        const bool sea = world.sea(poix, poiy) == 1;
        const bool river = !sea && world.riveraveflow(poix, poiy) > 0 && !lake;

        int longdegrees = 0;
        int longminutes = 0;
        int longseconds = 0;
        bool longneg = false;
        world.longitude(poix, longdegrees, longminutes, longseconds, longneg);
        const string longtext = formatcoordinate(longdegrees, longminutes, longseconds, longneg, degree);

        int latdegrees = 0;
        int latminutes = 0;
        int latseconds = 0;
        bool latneg = false;
        world.latitude(poiy, latdegrees, latminutes, latseconds, latneg);
        const string lattext = formatcoordinate(latdegrees, latminutes, latseconds, latneg, degree);

        const string windtext = getwindtext(world.wind(poix, poiy));

        string elevationtext;
        int pointelevation = world.map(poix, poiy);

        if (world.lakesurface(poix, poiy) != 0)
            pointelevation = world.lakesurface(poix, poiy) - 1;

        if (!sea)
        {
            if (pointelevation > sealevel)
            {
                elevationtext = formatnumber(pointelevation - sealevel) + " metres";

                if (world.seatotal() != 0)
                    elevationtext = elevationtext + " above sea level";
            }
            else if (lake)
            {
                elevationtext = formatnumber(sealevel - pointelevation) + " metres below sea level";
            }
            else if (world.seatotal() != 0)
            {
                elevationtext = formatnumber(world.lakesurface(poix, poiy) - sealevel) + " metres above sea level";
            }
            else
            {
                elevationtext = formatnumber(world.lakesurface(poix, poiy) - sealevel) + " metres";
            }
        }
        else
        {
            elevationtext = formatnumber(sealevel - pointelevation) + " metres below sea level";
        }

        const short climatetype = world.climate(poix, poiy);
        const string climatetext = getclimatename(climatetype) + " (" + getclimatecode(climatetype) + ")";
        const string specialstext = getglobalspecialtext(world, poix, poiy, sea);

        string glacialtext;
        if (!sea && (world.jantemp(poix, poiy) + world.jultemp(poix, poiy)) / 2 <= world.glacialtemp())
            glacialtext = "Ancient glacial region";

        string janflowtext;
        string julflowtext;
        string flowdirtext;

        if (river)
        {
            flowdirtext = getdirstring(world.riverdir(poix, poiy));
            janflowtext = formatnumber(world.riverjan(poix, poiy)) + " m" + cube + "/s";
            julflowtext = formatnumber(world.riverjul(poix, poiy)) + " m" + cube + "/s";
        }

        string lakedepthtext;
        if (lake)
        {
            elevationtext = formatnumber(world.lakesurface(poix, poiy) - sealevel) + " metres";

            if (world.seatotal() != 0)
                elevationtext = elevationtext + " above sea level";

            const int depth = world.lakesurface(poix, poiy) - world.nom(poix, poiy);
            const string salt = world.special(poix, poiy) == 100 ? " (salty)" : "";
            lakedepthtext = formatnumber(depth) + " metres" + salt;
        }

        const string seaicetext = formatseaicetext(world.seaice(poix, poiy));

        ImGui::Text("Longitude:");
        ImGui::SameLine(layout.pos2);
        ImGui::Text(longtext.c_str());

        ImGui::SameLine(layout.pos3);
        ImGui::Text("Elevation:");
        ImGui::SameLine(layout.pos4);
        ImGui::Text(elevationtext.c_str());

        if (lake)
        {
            ImGui::SameLine(layout.pos5);
            ImGui::Text(("Lake depth:  " + lakedepthtext).c_str());
        }

        if (river)
        {
            ImGui::SameLine(layout.pos5);
            ImGui::Text("River direction:");
            ImGui::SameLine(layout.pos6);
            ImGui::Text(flowdirtext.c_str());
        }

        ImGui::SameLine(layout.pos7);
        ImGui::Text("Charts:");

        ImGui::Text("Latitude:");
        ImGui::SameLine(layout.pos2);
        ImGui::Text(lattext.c_str());

        ImGui::SameLine(layout.pos3);
        ImGui::Text("Prevailing wind:");
        ImGui::SameLine(layout.pos4);
        ImGui::Text(windtext.c_str());

        if (river)
        {
            ImGui::SameLine(layout.pos5);
            ImGui::Text("January flow:");
            ImGui::SameLine(layout.pos6);
            ImGui::Text(janflowtext.c_str());
        }

        ImGui::SameLine(layout.pos7 + 20.0f);
        if (ImGui::Button("Temperature", ImVec2(120.0f, 0.0f)))
            toggle(showtemperaturechart);

        ImGui::Text(" ");

        if (river)
        {
            ImGui::SameLine(layout.pos5);
            ImGui::SetCursorPosY(ImGui::GetCursorPosY() - 4.0f);
            ImGui::Text("July flow:");
            ImGui::SameLine(layout.pos6);
            ImGui::SetCursorPosY(ImGui::GetCursorPosY() - 4.0f);
            ImGui::Text(julflowtext.c_str());
        }

        ImGui::SameLine(layout.pos7 + 20.0f);
        if (ImGui::Button("Precipitation", ImVec2(120.0f, 0.0f)))
            toggle(showrainfallchart);

        if (sea)
        {
            ImGui::Text("Sea ice:");
            ImGui::SameLine(layout.pos2);
            ImGui::Text(seaicetext.c_str());
        }
        else
        {
            ImGui::Text("Climate:");
            ImGui::SameLine(layout.pos2);
            ImGui::Text(climatetext.c_str());
        }

        ImGui::Text(" ");
        ImGui::SameLine(layout.pos2);

        if (specialstext != "")
            ImGui::Text(specialstext.c_str());
        else
            ImGui::Text(glacialtext.c_str());
    }

    if (newworld)
        ImGui::Text("Welcome to a new world!");

    drawinfowindowfooter(layout, showabout);
    ImGui::End();
}

void drawregionalinfowindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, region& region, bool focused, int poix, int poiy, const string& degree, const string& cube, bool& showtemperaturechart, bool& showrainfallchart, bool& showabout)
{
    const infopanellayout layout;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 180, main_viewport->WorkPos.y + 539), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(1023, 155), ImGuiCond_FirstUseEver);
    ImGui::Begin("Information##regional ", NULL, window_flags);

    if (focused)
    {
        const int sealevel = world.sealevel();
        const int width = world.width();
        const int height = world.height();
        const int xx = region.leftx() + poix / 16;
        const int yy = region.lefty() + poiy / 16;
        const bool lake = region.truelake(poix, poiy) != 0;
        const bool sea = region.sea(poix, poiy) == 1;
        const bool river = !lake && region.map(poix, poiy) > sealevel && (region.riveraveflow(poix, poiy) > 0 || region.fakeaveflow(poix, poiy) > 0);
        const bool beach = region.sand(poix, poiy) || region.shingle(poix, poiy) || region.mud(poix, poiy);

        int longdegrees = 0;
        int longminutes = 0;
        int longseconds = 0;
        bool longneg = false;
        region.longitude(poix, xx, poix / 16, width, longdegrees, longminutes, longseconds, longneg);
        const string longtext = formatcoordinate(longdegrees, longminutes, longseconds, longneg, degree);

        int latdegrees = 0;
        int latminutes = 0;
        int latseconds = 0;
        bool latneg = false;
        region.latitude(poiy, yy, poiy / 16, height, latdegrees, latminutes, latseconds, latneg);
        const string lattext = formatcoordinate(latdegrees, latminutes, latseconds, latneg, degree);

        const string windtext = getwindtext(world.wind(xx, yy));

        string elevationtext;
        int pointelevation = region.map(poix, poiy);

        if (region.lakesurface(poix, poiy) != 0)
            pointelevation = region.lakesurface(poix, poiy);

        if (!sea)
        {
            if (pointelevation > sealevel)
            {
                elevationtext = formatnumber(pointelevation - sealevel) + " metres";

                if (world.seatotal() != 0)
                    elevationtext = beach ? "Sea level" : elevationtext + " above sea level";
            }
            else if (lake)
            {
                elevationtext = formatnumber(sealevel - pointelevation) + " metres below sea level";
            }
            else if (world.seatotal() != 0)
            {
                elevationtext = formatnumber(region.lakesurface(poix, poiy) - sealevel) + " metres above sea level";
            }
            else
            {
                elevationtext = formatnumber(region.lakesurface(poix, poiy) - sealevel) + " metres";
            }
        }
        else
        {
            elevationtext = formatnumber(sealevel - pointelevation) + " metres below sea level";
        }

        const short climatetype = region.climate(poix, poiy);
        const string climatetext = getclimatename(climatetype) + " (" + getclimatecode(climatetype) + ")";
        const string specialstext = getregionalspecialtext(region, poix, poiy, sea, river);
        const string glacialtext;

        string janflowtext;
        string julflowtext;
        string flowdirtext;

        if (river)
        {
            flowdirtext = getdirstring(region.riverdir(poix, poiy));

            int janflow = region.riverjan(poix, poiy);
            int julflow = region.riverjul(poix, poiy);

            if (flowdirtext == "")
            {
                flowdirtext = getdirstring(region.fakedir(poix, poiy));
                janflow = region.fakejan(poix, poiy);
                julflow = region.fakejul(poix, poiy);
            }

            janflowtext = formatnumber(janflow) + " m" + cube + "/s";
            julflowtext = formatnumber(julflow) + " m" + cube + "/s";
        }

        string lakedepthtext;
        if (lake)
        {
            elevationtext = formatnumber(region.lakesurface(poix, poiy) - sealevel) + " metres";

            if (world.seatotal() != 0)
                elevationtext = elevationtext + " above sea level";

            const int depth = region.lakesurface(poix, poiy) - region.map(poix, poiy);
            const string salt = region.special(poix, poiy) == 100 ? " (salty)" : "";
            lakedepthtext = formatnumber(depth) + " metres" + salt;
        }

        const string seaicetext = formatseaicetext(region.seaice(poix, poiy));

        ImGui::Text("Longitude:");
        ImGui::SameLine(layout.pos2);
        ImGui::Text(longtext.c_str());

        ImGui::SameLine(layout.pos3);
        ImGui::Text("Elevation:");
        ImGui::SameLine(layout.pos4);
        ImGui::Text(elevationtext.c_str());

        if (lake)
        {
            ImGui::SameLine(layout.pos5);
            ImGui::Text(("Lake depth:  " + lakedepthtext).c_str());
        }

        if (river)
        {
            ImGui::SameLine(layout.pos5);
            ImGui::Text("River direction:");
            ImGui::SameLine(layout.pos6);
            ImGui::Text(flowdirtext.c_str());
        }

        ImGui::SameLine(layout.pos7);
        ImGui::Text("Charts:");

        ImGui::Text("Latitude:");
        ImGui::SameLine(layout.pos2);
        ImGui::Text(lattext.c_str());

        ImGui::SameLine(layout.pos3);
        ImGui::Text("Prevailing wind:");
        ImGui::SameLine(layout.pos4);
        ImGui::Text(windtext.c_str());

        if (river)
        {
            ImGui::SameLine(layout.pos5);
            ImGui::Text("January flow:");
            ImGui::SameLine(layout.pos6);
            ImGui::Text(janflowtext.c_str());
        }

        ImGui::SameLine(layout.pos7 + 20.0f);
        if (ImGui::Button("Temperature", ImVec2(120.0f, 0.0f)))
            toggle(showtemperaturechart);

        ImGui::Text(" ");

        if (river)
        {
            ImGui::SameLine(layout.pos5);
            ImGui::SetCursorPosY(ImGui::GetCursorPosY() - 4.0f);
            ImGui::Text("July flow:");
            ImGui::SameLine(layout.pos6);
            ImGui::SetCursorPosY(ImGui::GetCursorPosY() - 4.0f);
            ImGui::Text(julflowtext.c_str());
        }

        ImGui::SameLine(layout.pos7 + 20.0f);
        if (ImGui::Button("Precipitation", ImVec2(120.0f, 0.0f)))
            toggle(showrainfallchart);

        if (sea)
        {
            ImGui::Text("Sea ice:");
            ImGui::SameLine(layout.pos2);
            ImGui::Text(seaicetext.c_str());
        }
        else
        {
            ImGui::Text("Climate:");
            ImGui::SameLine(layout.pos2);
            ImGui::Text(climatetext.c_str());
        }

        ImGui::Text(" ");
        ImGui::SameLine(layout.pos2);

        if (specialstext != "")
            ImGui::Text(specialstext.c_str());
        else
            ImGui::Text(glacialtext.c_str());
    }

    drawinfowindowfooter(layout, showabout);
    ImGui::End();
}
