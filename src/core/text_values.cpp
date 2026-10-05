// Text forms of values: numbers with separators, values read from saved text, direction names.
// Moved unchanged from misc.cpp (Jonathan Hill, 2022).

#include "core/text_values.hpp"

#include <cmath>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include "core/classes.hpp"
#include "core/deterministic_random.hpp"
#include "core/grid_coordinates.hpp"
#include "core/procedural_noise.hpp"

using namespace std;

// This toggles a bool.

void toggle(bool &val)
{
    if (val == 0)
    {
        val = 1;
        return;
    }

    val = 0;
    return;
}

// This converts a number to a string containing commas. Taken from here: https://9to5answer.com/c-format-number-with-commas

string formatnumber(int val)
{
    string s = to_string(val);

    if (val<1000 && val>-1000)
        return s;
    int n = s.length() - 3;
    while (n > 0)
    {
        s.insert(n, ",");
        n -= 3;
    }

    return s;
}

// This converts strings into bools.

bool stob(string const& instring)
{
    if (instring == "1")
        return 1;

    return 0;
}

// This converts strings into shorts. (Note: no error trapping.)

short stos(string const& instring)
{
    int val = stoi(instring);

    short s = (short)val;

    return s;
}

// This converts strings into unsigned shorts. (Note: no error trapping.)

unsigned short stous(string const& instring)
{
    int val = stoi(instring);

    unsigned short s = (unsigned short)val;

    return s;
}

// This converts strings into chars. (Note: no error trapping.)

char stoc(string const& instring)
{
    int val = stoi(instring);

    char s = (char)val;

    return s;
}

// This converts strings into unsigned chars. (Note: no error trapping.)

unsigned char stouc(string const& instring)
{
    int val = stoi(instring);

    unsigned char s = (unsigned char)val;

    return s;
}

// This function translates a direction from an integer into a text string.

string getdirstring(int dir)
{
    string direction;

    switch (dir)
    {
    case 1:
        direction = "north";
        break;

    case 2:
        direction = "northeast";
        break;

    case 3:
        direction = "east";
        break;

    case 4:
        direction = "southeast";
        break;

    case 5:
        direction = "south";
        break;

    case 6:
        direction = "southwest";
        break;

    case 7:
        direction = "west";
        break;

    case 8:
        direction = "northwest";
        break;
    }

    return (direction);
}
