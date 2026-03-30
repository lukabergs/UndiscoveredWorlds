#include <string>
#include <vector>

#include <SFML/Graphics.hpp>

#include "map_imports.hpp"
#include "functions.hpp"

using namespace std;

namespace
{
    bool loadworldsizedimage(planet& world, const string& filepathname, sf::Image& importimage)
    {
        if (!importimage.loadFromFile(filepathname))
            return false;

        const sf::Vector2u imagesize = importimage.getSize();

        return imagesize.x == world.width() + 1 && imagesize.y == world.height() + 1;
    }
}

bool importlandheightmap(planet& world, const string& filepathname)
{
    sf::Image importimage;

    if (!loadworldsizedimage(world, filepathname, importimage))
        return false;

    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const sf::Color colour = importimage.getPixel(i, j);

            if (colour.r != 0)
                world.setnom(i, j, colour.r * 10 + sealevel);
            else
                world.setnom(i, j, sealevel - 5000);
        }
    }

    return true;
}

bool importseadepthmap(planet& world, const string& filepathname)
{
    sf::Image importimage;

    if (!loadworldsizedimage(world, filepathname, importimage))
        return false;

    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const sf::Color colour = importimage.getPixel(i, j);

            if (colour.r != 0)
            {
                int elev = sealevel - colour.r * 50;

                if (elev < 1)
                    elev = 1;

                world.setnom(i, j, elev);
            }
        }
    }

    return true;
}

bool importmountainmap(planet& world, const string& filepathname, vector<vector<bool>>& okmountains)
{
    sf::Image importimage;

    if (!loadworldsizedimage(world, filepathname, importimage))
        return false;

    const int width = world.width();
    const int height = world.height();
    vector<vector<int>> rawmountains(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const sf::Color colour = importimage.getPixel(i, j);
            rawmountains[i][j] = colour.r != 0 ? colour.r * 65 : 0;
        }
    }

    createmountainsfromraw(world, rawmountains, okmountains);
    return true;
}

bool importvolcanomap(planet& world, const string& filepathname)
{
    sf::Image importimage;

    if (!loadworldsizedimage(world, filepathname, importimage))
        return false;

    const int width = world.width();
    const int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const sf::Color colour = importimage.getPixel(i, j);

            if (colour.r != 0)
            {
                int elev = colour.r * 45;
                bool strato = colour.g > 0;

                if (!strato)
                    elev = elev / 2;

                if (colour.b == 0)
                    elev = 0 - elev;

                world.setvolcano(i, j, elev);
                world.setstrato(i, j, strato);
            }
        }
    }

    return true;
}
