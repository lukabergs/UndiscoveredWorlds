// Shapes on plain grids: lines, circles, boxes, flood fills, smoothing, curves, Voronoi cells.
// Moved unchanged from misc.cpp (Jonathan Hill, 2022).

#include "core/grid_shapes.hpp"

#include <cmath>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include "core/classes.hpp"
#include "core/deterministic_random.hpp"
#include "core/grid_coordinates.hpp"
#include "core/procedural_noise.hpp"

#define ARRAYWIDTH 2048
#define ARRAYHEIGHT 1024

using namespace std;

// This function creates a rift blob template.

void createriftblob(vector<vector<float>>& riftblob, int size)
{
    for (int i = 0; i <= size * 2; i++)
    {
        for (int j = 0; j <= size * 2; j++)
            riftblob[i][j] = 0;
    }

    float startvalue = 1.0f; // The starting value for the centre of the circle.
    float endvalue = 0.15f; // The end value for the edge of the circle.
    float step = (startvalue - endvalue) / size; // Reduce the value by this much at each step of the circle drawing.

    float currentvalue = startvalue;

    int centrex = size + 1;
    int centrey = size + 1;

    for (int radius = 1; radius <= size; radius++)
    {
        for (int i = -radius; i <= radius; i++)
        {
            int ii = centrex + i;

            for (int j = -radius; j <= radius; j++)
            {
                int jj = centrey + j;

                if (riftblob[ii][jj] == 0 && i * i + j * j < radius * radius + radius)
                    riftblob[ii][jj] = currentvalue;
            }
        }
        currentvalue = currentvalue - step;
    }
}

// This normalises one value to another, so they are as close as possible, given the wrapped maximum.

int normalise(int x, int y, int max)
{
    while (abs(x - y) > abs(x - (y - max)))
        y = y - max;

    while (abs(x - y) > abs(x - (y + max)))
        y = y + max;

    return y;
}

// This function smoothes a vector.

void smooth(vector<vector<int>>& arr, int width, int height, int maxelev, int amount, bool vary)
{
    int varyx = width / 2;
    int varyy = random(1, height - 1);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 1; j < height; j++)
        {
            bool goaheadx = 0;
            bool goaheady = 0;

            if (vary == 1) // This will vary the amount of smoothing over the map, to make it more interesting.
            {
                if (i < varyx)
                {
                    if (i > 10 && random(0, varyx) < i * 2)
                        goaheadx = 1;
                }
                else
                {
                    if (i<width - 10 && random(varyx, width)>i / 2)
                        goaheadx = 1;
                }

                if (j < varyy)
                {
                    if (random(0, varyy) < j * 2)
                        goaheady = 1;
                }
                else
                {
                    if (random(varyy, height) > j / 2)
                        goaheady = 1;
                }
                if (abs(varyx - i) < 100) // Force it to smooth in the centre of the map.
                {
                    //goaheadx=1;
                    //goaheady=1;
                }
            }
            else
            {
                goaheadx = 1;
                goaheady = 1;
            }

            if (goaheadx == 1 && goaheady == 1)
            {
                float crount = 0.0f;
                float ave = 0.0f;

                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk < 0)
                        kk = width;

                    if (kk > width)
                        kk = 0;

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        //ave=ave+nomwrap(k,l);
                        ave = ave + (float)arr[kk][l];
                        crount++;
                    }
                }

                ave = ave / crount;

                if (ave > 0.0f && ave < (float)maxelev)
                    arr[i][j] = (int)ave;
            }
        }
    }
}

// This function sees whether the given point on the given vector is on the edge of the filled area.
// Note, this wraps around the edges.

bool edge(vector<vector<bool>>& arr, int width, int height, int i, int j)
{
    if (arr[i][j] == 0)
        return 0;

    for (int k = i - 1; k <= i + 1; k++)
    {
        int kk = k;

        if (kk<0 || kk>width)
            kk = wrap(kk, width);

        for (int l = j - 1; l <= j + 1; l++)
        {
            if (l >= 0 && l <= height)
            {
                if (arr[kk][l] == 0)
                    return 1;
            }
        }
    }

    return 0;
}

// This draws a straight line on the supplied vector.
// Uses Bresenham's line algorithm - http://rosettacode.org/wiki/Bitmap/Bresenham%27s_line_algorithm#C.2B.2B

void drawline(vector<vector<bool>>& arr, int x1, int y1, int x2, int y2)
{
    const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
    if (steep)
    {
        std::swap(x1, y1);
        std::swap(x2, y2);
    }

    if (x1 > x2)
    {
        std::swap(x1, x2);
        std::swap(y1, y2);
    }

    const float dx = (float)(x2 - x1);
    const float dy = (float)fabs(y2 - y1);

    float error = dx / 2.0f;
    const int ystep = (y1 < y2) ? 1 : -1;
    int y = (int)y1;

    const int maxX = (int)x2;

    for (int x = (int)x1; x <= maxX; x++)
    {
        if (steep)
            arr[y][x] = 1;
        else
            arr[x][y] = 1;

        error -= dy;
        if (error < 0)
        {
            y += ystep;
            error += dx;
        }
    }
}

// This does the same thing, but allows lines to pass between the left and right sides.

void drawlinewrapped(vector<vector<bool>>& arr, int width, int height, int x1, int y1, int x2, int y2)
{
    if (x2 > x1)
    {
        int dist = x2 - x1;

        int altdist = (x1 + width) - x2;

        if (altdist < dist)
            x1 = x1 + width;
    }
    else
    {
        int dist = x1 - x2;

        int altdist = (x2 + width) - x1;

        if (altdist < dist)
            x2 = x2 + width;
    }
    
    const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
    if (steep)
    {
        std::swap(x1, y1);
        std::swap(x2, y2);
    }

    if (x1 > x2)
    {
        std::swap(x1, x2);
        std::swap(y1, y2);
    }

    const float dx = (float)(x2 - x1);
    const float dy = (float)fabs(y2 - y1);

    float error = dx / 2.0f;
    const int ystep = (y1 < y2) ? 1 : -1;
    int y = (int)y1;

    const int maxX = (int)x2;

    for (int x = (int)x1; x <= maxX; x++)
    {
        int a = x;
        int b = y;

        if (steep)
        {
            a = y;
            b = x;
        }

        if (a < 0)
            a = a + width;

        if (a > width)
            a = a - width;

        if (b < 0)
            b = b + width;

        if (b > width)
            b = b - width;
        
        arr[a][b] = 1;

        error -= dy;
        if (error < 0)
        {
            y += ystep;
            error += dx;
        }
    }
}

// This function draws a filled circle on the specified vector.
// Uses code by OBese87 and Libervurto - https://forum.thegamecreators.com/thread/197341

void drawcircle(vector<vector<int>>& arr, int x, int y, int col, int radius)
{
    float dots = radius * 6.28f;
    float t = 360.0f / dots;
    int quarter = (int)(dots / 4.0f);

    for (int i = 1; i <= quarter; i++)
    {
        int u = (int)((float)sin(i * t) * (float)radius);
        int v = (int)((float)cos(i * t) * (float)radius);

        box(arr, x - u, y - v, x + u, y + v, col);
    }
}

// This function draws a filled box on the specified vector.

void box(vector<vector<int>>& arr, int x1, int y1, int x2, int y2, int col)
{
    for (int i = x1; i <= x2; i++)
    {
        int ii = i;

        if (ii<0 || i>ARRAYWIDTH - 1)
            ii = wrap(ii, ARRAYWIDTH - 1);

        for (int j = y1; j <= y2; j++)
        {
            if (j >= 0 && j < ARRAYHEIGHT)
                arr[ii][j] = col;
        }
    }
}

// As above, but on a three-dimensional vector.

void drawcircle3d(vector<vector<vector<int>>>& arr, int x, int y, int col, int radius, int index)
{
    float dots = radius * 6.28f;
    float t = 360.0f / dots;
    int quarter = (int)(dots / 4.0f);

    for (int i = 1; i <= quarter; i++)
    {
        int u = (int)((float)sin(i * t) * (float)radius);
        int v = (int)((float)cos(i * t) * (float)radius);

        box3d(arr, x - u, y - v, x + u, y + v, col, index);
    }
}

// Again, on a three-dimensional vector.

void box3d(vector<vector<vector<int>>>& arr, int x1, int y1, int x2, int y2, int col, int index)
{
    for (int i = x1; i <= x2; i++)
    {
        for (int j = y1; j <= y2; j++)
        {
            if (i >= 0 && i <= 20 && j >= 0 && j <= 20)
                arr[index][i][j] = col;
        }
    }
}

// This function draws splines. Based on functions by Markus - https://forum.thegamecreators.com/thread/202580
// (This is easily the most important bit of borrowed code in this project. Thank you Markus!)

twofloats curvepos(twofloats p0, twofloats p1, twofloats p2, twofloats p3, float t)
{
    twofloats pt;

    pt.x = 0.5f * ((2.0f * p1.x) + (-p0.x + p2.x) * t + (2.0f * p0.x - 5.0f * p1.x + 4.0f * p2.x - p3.x) * (t * t) + (-p0.x + 3.0f * p1.x - 3.0f * p2.x + p3.x) * (t * t * t)); 
    pt.y = 0.5f * ((2.0f * p1.y) + (-p0.y + p2.y) * t + (2.0f * p0.y - 5.0f * p1.y + 4.0f * p2.y - p3.y) * (t * t) + (-p0.y + 3.0f * p1.y - 3.0f * p2.y + p3.y) * (t * t * t));

    return (pt);
}

// This function does a flood fill on a bool vector.

void fill(vector<vector<bool>>& arr, int width, int height, int startx, int starty, bool replacement)
{
    if (startx < 0 || startx >= width || starty < 0 || starty >= height)
        return;

    if (arr[startx][starty] == replacement)
        return;

    int row[] = { -1,0,0,1 };
    int col[] = { 0,-1,1,0 };

    bool target = arr[startx][starty]; // This is the "colour" that we're going to change into the replacement "colour".

    twointegers node;

    node.x = startx;
    node.y = starty;

    queue<twointegers> q; // Create a queue
    q.push(node); // Put our starting node into the queue

    while (q.empty() != 1) // Keep going while there's anything left in the queue
    {
        node = q.front(); // Take the node at the front of the queue
        q.pop(); // And then pop it out of the queue

        if (node.x >= 0 && node.x < width && node.y >= 0 && node.y < height && arr[node.x][node.y] == target)
        {
            arr[node.x][node.y] = replacement;

            for (int k = 0; k < 4; k++) // Look at the four neighbouring nodes in turn
            {
                twointegers nextnode;

                nextnode.x = node.x + row[k];
                nextnode.y = node.y + col[k];

                if (nextnode.x >= 0 && nextnode.x < width && nextnode.y >= 0 && nextnode.y < height)
                {
                    if (arr[nextnode.x][nextnode.y] == target) // If this node is the right "colour"
                        q.push(nextnode); // Put that node onto the queue
                }
            }
        }
    }
}

// This function does the same thing, but with a continental mask to avoid going into the sea.

void fillcontinent(vector<vector<bool>>& arr, vector<vector<short>>& mask, short maskcheck, int width, int height, int startx, int starty, bool replacement)
{
    if (startx < 0 || startx >= width || starty < 0 || starty >= height)
        return;

    if (arr[startx][starty] == replacement)
        return;

    if (mask[startx][starty] != maskcheck)
        return;

    int row[] = { -1,0,0,1 };
    int col[] = { 0,-1,1,0 };

    bool target = arr[startx][starty]; // This is the "colour" that we're going to change into the replacement "colour".

    twointegers node;

    node.x = startx;
    node.y = starty;

    queue<twointegers> q; // Create a queue
    q.push(node); // Put our starting node into the queue

    while (q.empty() != 1) // Keep going while there's anything left in the queue
    {
        node = q.front(); // Take the node at the front of the queue
        q.pop(); // And then pop it out of the queue

        if (node.x >= 0 && node.x < width && node.y >= 0 && node.y < height && arr[node.x][node.y] == target)
        {
            arr[node.x][node.y] = replacement;

            for (int k = 0; k < 4; k++) // Look at the four neighbouring nodes in turn
            {
                twointegers nextnode;

                nextnode.x = node.x + row[k];
                nextnode.y = node.y + col[k];

                if (nextnode.x >= 0 && nextnode.x < width && nextnode.y >= 0 && nextnode.y < height)
                {
                    if (arr[nextnode.x][nextnode.y] == target && mask[nextnode.x][nextnode.y] == maskcheck) // If this node is the right "colour"
                        q.push(nextnode); // Put that node onto the queue
                }
            }
        }
    }
}

// This function creates a voronoi map.

void makevoronoi(vector<vector<short>>& voronoi, int width, int height, int points)
{
    vector<vector<int>> pointslist(points, vector<int>(2, 0));

    int margin = points / 6;

    int widthmargin = width / 10;
    int widthfarmargin = width - widthmargin;
    int heightmargin = height / 10;
    int heightfarmargin = height - heightmargin;

    for (int n = 0; n <= margin; n = n + 3) // Force it to put some points near the edges.
    {
        if (random(1, 2) == 1)
        {
            pointslist[n][0] = random(0, widthmargin);
            pointslist[n][1] = random(0, height);
        }
        else
        {
            pointslist[n][0] = random(widthfarmargin, width);
            pointslist[n][1] = random(0, height);
        }

        pointslist[n + 1][0] = random(0, width);
        pointslist[n + 1][1] = random(0, heightmargin);

        pointslist[n + 2][0] = random(0, width);
        pointslist[n + 2][1] = random(heightfarmargin, height);

    }

    for (int n = margin + 1; n < points; n++)
    {
        pointslist[n][0] = random(0, width);
        pointslist[n][1] = random(0, height);
    }

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            int closest = 0;
            int dist = 10000000;

            for (int point = 0; point < points; point++)
            {
                int xdist = i - pointslist[point][0];
                int ydist = j - pointslist[point][1];
                int pointdist = xdist * xdist + ydist * ydist;

                if (pointdist < dist)
                {
                    closest = point;
                    dist = pointdist;
                }
            }

            voronoi[i][j] = closest;
        }
    }
}

std::vector<mesh::CellId> columnmajorcells(const mesh::CellGrid& grid)
{
    std::vector<mesh::CellId> cells;
    cells.reserve(grid.cellcount());
    if (grid.kind() == mesh::GridKind::latlon && !grid.issubgrid())
    {
        const int columns = grid.columns();
        const int rows = grid.rows();
        for (int x = 0; x < columns; x++)
            for (int y = 0; y < rows; y++)
                cells.push_back(static_cast<mesh::CellId>(y) * static_cast<mesh::CellId>(columns) + static_cast<mesh::CellId>(x));
    }
    else if (grid.kind() == mesh::GridKind::cube && !grid.issubgrid())
    {
        const mesh::CellId n = static_cast<mesh::CellId>(grid.facecells());
        for (mesh::CellId face = 0; face < 6; face++)
            for (mesh::CellId i = 0; i < n; i++)
                for (mesh::CellId j = 0; j < n; j++)
                    cells.push_back(face * n * n + j * n + i);
    }
    else
    {
        for (mesh::CellId c = 0; c < grid.cellcount(); c++)
            cells.push_back(c);
    }
    return cells;
}

void legacyslot(const mesh::CellGrid& grid, mesh::CellId c, int& x, int& y)
{
    if (grid.kind() == mesh::GridKind::latlon)
    {
        x = grid.column(c);
        y = grid.row(c);
        return;
    }
    x = static_cast<int>(c % ARRAYWIDTH);
    y = static_cast<int>(c / ARRAYWIDTH);
}
