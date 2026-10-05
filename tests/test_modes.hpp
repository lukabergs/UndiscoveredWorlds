#pragma once

// Replay and benchmark modes of the climate tests (tier 4). Each is a disabled
// test that ctest never runs; its arguments come from UW_TEST_ARGUMENTS:
//
//   $env:UW_TEST_ARGUMENTS = "--atmosphere-benchmark 256"
//   world_tests --gtest_also_run_disabled_tests --gtest_filter=ClimateAtmosphere.DISABLED_Modes
//
// A mode returns 0 on success and -1 when the arguments name none of its modes.

#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

inline int runtestmode(int (*mode)(int, char**))
{
    std::vector<std::string> words{ "world_tests" };
    if (const char* text = std::getenv("UW_TEST_ARGUMENTS"))
    {
        std::istringstream stream(text);
        for (std::string word; stream >> word;)
            words.push_back(word);
    }
    std::vector<char*> arguments;
    for (std::string& word : words)
        arguments.push_back(word.data());
    return mode(static_cast<int>(arguments.size()), arguments.data());
}
