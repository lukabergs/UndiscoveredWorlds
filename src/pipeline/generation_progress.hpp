#pragma once

#include <string>

// Console progress for generation stages. With timed reporting on, each report
// line is closed with its elapsed time.
void updatereport(std::string text);
void begintimedreporting();
void endtimedreporting();
