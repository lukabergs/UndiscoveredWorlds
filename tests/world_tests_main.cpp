// world_tests runs every TEST by default. A few climate test files also carry
// explicit replay and benchmark modes (tier 4, never run by ctest); each
// returns -1 when the arguments do not name one of its modes.

#include <iostream>

#include <gtest/gtest.h>

int climateatmospheretool(int argc, char** argv);
int climateoceandynamicstool(int argc, char** argv);
int climateoceanheattool(int argc, char** argv);
int climatetransporttool(int argc, char** argv);

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc > 1)
    {
        for (const auto tool : {&climateatmospheretool, &climateoceandynamicstool, &climateoceanheattool, &climatetransporttool})
        {
            const int code = tool(argc, argv);
            if (code >= 0)
                return code;
        }
        // GoogleTest has removed its own flags; anything left is a mistyped mode.
        std::cerr << "Unknown arguments. Modes: --atmosphere-replay FIXTURE [TOLERANCE] | --atmosphere-benchmark WIDTH | "
                     "--momentum-replay FIXTURE CSV EKMAN GYRE | --coastal-replay FIXTURE DAYS [STEP_HOURS] | "
                     "--heat-kernels WIDTH DAYS | --accuracy WIDTH | --ocean-heat-benchmark WIDTH DAYS | --transport-benchmark CSV\n";
        return 2;
    }
    return RUN_ALL_TESTS();
}
