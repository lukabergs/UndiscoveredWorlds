#include "pipeline/generation_progress.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>

#include "pipeline/world_generation_debug.hpp"

using namespace std;

namespace
{
bool timedreportingenabled = false;
bool hastimedreport = false;
string timedreporttext;
chrono::steady_clock::time_point timedreportstart;

void flushtimedreport()
{
    if (hastimedreport == false)
        return;

    const double elapsedms = chrono::duration<double, milli>(chrono::steady_clock::now() - timedreportstart).count();
    const streamsize oldprecision = cout.precision();
    const ios::fmtflags oldflags = cout.flags();

    cout << '\r' << left << setw(48) << timedreporttext << " [" << fixed << setprecision(2) << elapsedms << "ms]" << endl;
    cout.flags(oldflags);
    cout.precision(oldprecision);

    hastimedreport = false;
    timedreporttext.clear();
}
}

void updatereport(string text)
{
    if (timedreportingenabled == false)
    {
        cout << text << endl;
        return;
    }

    flushtimedreport();

    if (text.empty())
    {
        cout << endl;
        return;
    }

    timedreporttext = text;
    timedreportstart = chrono::steady_clock::now();
    hastimedreport = true;

    cout << text << "..." << flush;
}

void begintimedreporting()
{
    timedreportingenabled = true;
    hastimedreport = false;
    timedreporttext.clear();
}

void endtimedreporting()
{
    flushtimedreport();
    timedreportingenabled = false;
}
