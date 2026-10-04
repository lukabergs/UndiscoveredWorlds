// world_tests runs every TEST by default. A few climate test files also carry
// explicit replay and benchmark modes (tier 4, never run by ctest); each
// returns -1 when the arguments do not name one of its modes.

#include <gtest/gtest.h>

int climateatmospheretool(int argc, char** argv);
int climateoceandynamicstool(int argc, char** argv);
int climateoceanheattool(int argc, char** argv);
int climatetransporttool(int argc, char** argv);

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc > 1)
        for (const auto tool : {&climateatmospheretool, &climateoceandynamicstool, &climateoceanheattool, &climatetransporttool})
        {
            const int code = tool(argc, argv);
            if (code >= 0)
                return code;
        }
    return RUN_ALL_TESTS();
}
