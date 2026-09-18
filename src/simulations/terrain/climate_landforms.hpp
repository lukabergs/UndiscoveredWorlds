#pragma once

// Climate-driven dunes, salt pans, river-adjacent terrain and surface roughness.

class planet;
class boolshapetemplate;

void broadenfastlemterrainfromrivers(planet& world);
void createergs(planet& world, boolshapetemplate smalllake[], boolshapetemplate largelake[], boolshapetemplate shape[]);
void createsaltpans(planet& world, boolshapetemplate smalllake[], boolshapetemplate largelake[]);
void refineroughnessmap(planet& world);
