#include "land_distance.hpp"
#include "planet.hpp"
#include "functions.hpp"
#include <queue>

using namespace std;

namespace terraindetail
{
bool buildlanddistancefield(planet& world, vector<vector<short>>& distances, int maximumdistance, const std::function<bool(int, int)>& issource)
{
    const int width = world.width();
    const int height = world.height();
    queue<pair<int, int>> frontier;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            distances[i][j] = static_cast<short>(maximumdistance + 1);

            if (world.sea(i, j) == 0 && issource(i, j))
            {
                distances[i][j] = 0;
                frontier.emplace(i, j);
            }
        }
    }

    if (frontier.empty())
        return false;

    const int directions[8][2] =
    {
        { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 },
        { 0, 1 }, { -1, 1 }, { -1, 0 }, { -1, -1 }
    };

    while (frontier.empty() == false)
    {
        const pair<int, int> current = frontier.front();
        frontier.pop();

        const short currentdistance = distances[current.first][current.second];

        if (currentdistance >= maximumdistance)
            continue;

        for (const auto& direction : directions)
        {
            const int nx = wrap(current.first + direction[0], width);
            const int ny = current.second + direction[1];

            if (ny < 0 || ny > height || world.sea(nx, ny) == 1)
                continue;

            const short nextdistance = static_cast<short>(currentdistance + 1);

            if (distances[nx][ny] <= nextdistance)
                continue;

            distances[nx][ny] = nextdistance;
            frontier.emplace(nx, ny);
        }
    }

    return true;
}
}
