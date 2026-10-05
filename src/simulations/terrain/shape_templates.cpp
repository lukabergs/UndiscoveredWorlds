// Shape templates of the terrain module: peaks, mountain ranges, lakes, islands,
// smudges and land shapes. The data is assets/terrain/shape-templates.txt, one
// template per line ("key: number of values: values"); this file decodes it.

#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/app_environment.hpp"
#include "core/classes.hpp"
#include "world/functions.hpp"

using namespace std;

namespace
{
// Every template of the asset file, read once.
const map<string, vector<int>>& templatedata()
{
    static const map<string, vector<int>> data = []
    {
        const auto path = getappenvironment().terrainTemplatePath;
        ifstream file(path);
        if (!file)
            throw runtime_error("Cannot read the terrain shape templates: " + path.string());
        map<string, vector<int>> result;
        string line;
        while (getline(file, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            const size_t keyend = line.find(':');
            const size_t countend = line.find(':', keyend + 1);
            if (keyend == string::npos || countend == string::npos)
                throw runtime_error("Malformed line in " + path.string() + ": " + line.substr(0, 40));
            const size_t count = stoul(line.substr(keyend + 1, countend - keyend - 1));
            vector<int>& values = result[line.substr(0, keyend)];
            values.reserve(count);
            istringstream stream(line.substr(countend + 1));
            for (int value; stream >> value;)
                values.push_back(value);
            if (values.size() != count)
                throw runtime_error("Template " + line.substr(0, keyend) + " in " + path.string() + " does not hold its declared number of values");
        }
        return result;
    }();
    return data;
}

const vector<int>& templatevalues(const string& key)
{
    const auto found = templatedata().find(key);
    if (found == templatedata().end())
        throw runtime_error("Missing terrain shape template: " + key);
    return found->second;
}

// A square shape stored cell by cell: the span, then span x span values.
template <typename Shape>
void loadgridshapes(Shape shapes[], const string& family, int first, int last)
{
    for (int n = first; n <= last; n++)
    {
        const vector<int>& raw = templatevalues(family + "/" + to_string(n));

        int index = 0;
        int span = raw[index];

        shapes[n].setxsize(span);
        shapes[n].setysize(span);

        shapes[n].clear();

        for (int i = 0; i < span; i++)
        {
            for (int j = 0; j < span; j++)
            {
                index++;

                shapes[n].setpoint(i, j, raw[index]);
            }
        }
    }
}

// A square shape stored as runs: the span, then pairs of run length and value.
void loadrunshapes(boolshapetemplate shapes[], const string& family, int last)
{
    for (int n = 0; n <= last; n++)
    {
        const vector<int>& raw = templatevalues(family + "/" + to_string(n));

        int index = 0;
        int span = raw[index];
        int rawspan = static_cast<int>(raw.size());

        shapes[n].setxsize(span);
        shapes[n].setysize(span);

        shapes[n].clear();

        index = 1;

        int i = 0;
        int j = 0;

        bool alldone = 0;

        do
        {
            int total = raw[index];
            int val = raw[index + 1];

            for (int pets = 0; pets < total; pets++)
            {
                if (i < span && j < span)
                    shapes[n].setpoint(i, j, val);
                else
                    alldone = 1;

                j++;

                if (j >= span)
                {
                    j = 0;
                    i++;
                }
            }

            index = index + 2;

            if (index >= rawspan - 1)
                alldone = 1;

        } while (alldone == 0);
    }
}
}

void loadpeaktemplates(peaktemplate& peaks)
{
    for (int n = 0; n <= 3; n++)
    {
        const vector<int>& raw = templatevalues("peak/" + to_string(n));

        int index = 0;
        int span = raw[index];

        index++;

        int cx = raw[index];

        index++;

        int cy = raw[index];

        peaks.setspan(n, span);
        peaks.setcentrex(n, cx);
        peaks.setcentrey(n, cy);

        for (int i = 0; i <= span; i++)
        {
            for (int j = 0; j <= span; j++)
            {
                index++;

                peaks.setpeakmap(n, i, j, raw[index]);
            }
        }
    }
}

// Mountain range templates: tempno from 1 to MOUNTAINTEMPLATESTOTAL, each in eight directions.
// Directions are: 1 - ENE; 2 - ESE; 3 - SSE; 4 - SSW; 5 - WSW; 6 - WNW; 7 - NNW; 8 - NNE
// Returns the span, or 0 (templates untouched) for a template or direction outside those ranges.

int createmountainrangetemplate(vector<vector<unsigned char>>& dirtemplate, vector<vector<unsigned char>>& heighttemplate, int tempno, int dir)
{
    if (tempno < 1 || tempno > MOUNTAINTEMPLATESTOTAL || dir < 1 || dir > 8)
        return 0;

    const vector<int>& raw = templatevalues("mountain-range/" + to_string(tempno) + "/" + to_string(dir));

    int index = 0;
    int span = raw[index];

    for (int i = 0; i <= span; i++)
    {
        for (int j = 0; j <= span; j++)
        {
            index++;
            dirtemplate[i][j] = static_cast<unsigned char>(raw[index]);

            index++;
            heighttemplate[i][j] = static_cast<unsigned char>(raw[index]);
        }
    }

    return span;
}

int createsmalllaketemplates(boolshapetemplate smalllake[])
{
    loadgridshapes(smalllake, "small-lake", 0, 11);

    return 0;
}

int createislandtemplates(boolshapetemplate island[])
{
    loadgridshapes(island, "island", 0, 11);

    return 0;
}

int createsmudgetemplates(byteshapetemplate smudge[])
{
    loadgridshapes(smudge, "smudge", 1, 5);

    return 0;
}

int createsmallsmudgetemplates(byteshapetemplate smallsmudge[])
{
    loadgridshapes(smallsmudge, "small-smudge", 1, 5);

    return 0;
}

int createlargelaketemplates(boolshapetemplate largelake[])
{
    loadrunshapes(largelake, "large-lake", 9);

    return 0;
}

int createlandshapetemplates(boolshapetemplate landshape[])
{
    loadrunshapes(landshape, "land-shape", 11);

    return 0;
}

int createchainlandtemplates(boolshapetemplate chainland[])
{
    loadrunshapes(chainland, "chain-land", 1);

    return 0;
}
