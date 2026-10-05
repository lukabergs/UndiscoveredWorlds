#pragma once

// Shapes on plain grids: lines, circles, boxes, flood fills, smoothing, curves, Voronoi cells.

#include <string>
#include <vector>

#include "core/classes.hpp"

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
