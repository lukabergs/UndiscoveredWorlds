#pragma once

// Climate-driven dunes, salt pans and surface roughness.

class planet;
class boolshapetemplate;

void createergs(planet& world, boolshapetemplate smalllake[], boolshapetemplate largelake[], boolshapetemplate shape[]);
void createsaltpans(planet& world, boolshapetemplate smalllake[], boolshapetemplate largelake[]);
void refineroughnessmap(planet& world);
