#pragma once

// Seeded random numbers and procedural noise shared by terrain, climate and
// hydrology. No world state: callers pass arrays and sizes.

#include <vector>

// Shared legacy generator (seeded per stage).
void fast_srand(long seed);
int fast_rand(void);
int random(int a, int b);
int randomsign(int a);
float randomsign(float a);
// The C library generator, used by a few legacy paths.
int altrandom(int a, int b);
int altrandomsign(int a);

// Average of two values on a wrapped range, and linear blends by percentage.
int wrappedaverage(int x, int y, int max);
int tilt(int a, int b, int percentage);
float tilt(float a, float b, int percentage);

// Diamond-square fractal on [0..awidth] x [0..aheight].
void createfractal(std::vector<std::vector<int>>& arr, int awidth, int aheight, int grain, float valuemod, float valuemod2,
    int min, int max, bool extreme, bool wrapped);
void warp(std::vector<std::vector<int>>& map, int width, int height, int maxelev, int warpfactor, bool vary);
void monstrouswarp(std::vector<std::vector<int>>& map, int width, int height, int maxelev, int warpfactor);
void shift(std::vector<std::vector<int>>& map, int width, int height, int offset);
void flip(std::vector<std::vector<int>>& arr, int awidth, int aheight, bool vert, bool horiz);
