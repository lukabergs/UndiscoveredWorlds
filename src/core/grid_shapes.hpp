#pragma once

// Shapes on plain grids: lines, circles, boxes, flood fills, smoothing, curves, Voronoi cells.

#include <array>
#include <string>
#include <vector>

#include "core/classes.hpp"
#include "core/mesh/cell_grid.hpp"

void createriftblob(std::vector<std::vector<float>>& riftblob, int size);
int normalise(int x, int y, int max);
void smooth(std::vector<std::vector<int>>& arr, int width, int height, int maxelev, int amount, bool vary);
bool edge(std::vector<std::vector<bool>>& arr, int width, int height, int i, int j);
void drawline(std::vector<std::vector<bool>>& arr, int x1, int y1, int x2, int y2);
void drawlinewrapped(std::vector<std::vector<bool>>& arr, int width, int height, int x1, int y1, int x2, int y2);
void drawcircle(std::vector<std::vector<int>>& arr, int x, int y, int col, int radius);
void box(std::vector<std::vector<int>>& arr, int x1, int y1, int x2, int y2, int col);
void drawcircle3d(std::vector<std::vector<std::vector<int>>>& arr, int x, int y, int col, int radius, int index);
void box3d(std::vector<std::vector<std::vector<int>>>& arr, int x1, int y1, int x2, int y2, int col, int index);
twofloats curvepos(twofloats p0, twofloats p1, twofloats p2, twofloats p3, float t);
void fill(std::vector<std::vector<bool>>& arr, int width, int height, int x, int y, bool replacement);
void fillcontinent(std::vector<std::vector<bool>>& arr, std::vector<std::vector<short>>& mask, short maskcheck, int width, int height, int startx, int starty, bool replacement);
void makevoronoi(std::vector<std::vector<short>>& voronoi, int width, int height, int points);

// Cell-grid forms of the legacy (x, y) loops (docs/architecture/sphere-mesh.md, sections 3 and 4).

// The cells in the order of the legacy loops "for x { for y }": on lat-lon column by column,
// each from north to south; on the cube face by face, each face column by column (i outer, j
// inner); a subgrid in ascending id order. The order of "for y { for x }" is the id order.
std::vector<mesh::CellId> columnmajorcells(const mesh::CellGrid& grid);

// The slot [x][y] of a cell in a legacy ARRAYWIDTH x ARRAYHEIGHT (2048 x 1024) array: its column
// and row on lat-lon; elsewhere the id folded into rows of 2048 (the arrays of the generation
// scratch and the workbench hold a cube's cells that way until they store per-cell vectors).
void legacyslot(const mesh::CellGrid& grid, mesh::CellId c, int& x, int& y);

// The 3 x 3 box around a cell as neighbour slots (0: the cell itself) in the order of the legacy
// loops "for (i = x - 1; i <= x + 1) for (j = y - 1; j <= y + 1)": NW, W, SW, N, the cell, S, NE, E, SE.
constexpr std::array<int, 9> boxslots = { 8, 7, 6, 1, 0, 5, 2, 3, 4 };
