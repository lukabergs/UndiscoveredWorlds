// Tests of landscape/regional_solve.hpp (multiscale plan S3.0: the solve step of the hook
// interface, regional_pass.cpp's solve moved unchanged; S3.2 changes it), suite Landscape:
// TEST_F(Landscape, RegionalSolve<What>) on the synthetic world of regional_test_world.hpp, the
// level's steps up to the channels run through the hook interface (regional_pass.hpp).

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_solve.hpp"
#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_valleys.hpp"

TEST_F(Landscape, RegionalSolveKeepsThePinnedCellsAndTheBeds)
{
    // At refinement.channel_cell_slope 0 (S3.2's term off) the solve leaves every pinned cell
    // (channels, sea, the outer ring) at the working surface, so a land channel cell keeps its
    // bed, writes solved_m as the result and leaves bed_m unchanged.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalParameters parameters;
    const RegionalWindow window = windowa(4);
    const landscape::FineWarp warp = landscape::regionalwarp(window, inputs.seed, parameters);
    const landscape::GlobalIncision incision = landscape::globalincision(inputs, parameters);
    const landscape::RegionalDetails details{};
    const landscape::LevelContext context{inputs, window, parameters, warp, incision, details};
    landscape::LevelState state;
    landscape::noiselevel(context, state);
    landscape::presmoothinglevel(context, state);
    landscape::coastlevel(context, state);
    landscape::incisionlevel(context, state);
    landscape::channelslevel(context, state);
    const std::vector<float> surface = state.surface_m;
    const std::vector<float> beds = state.bed_m;

    landscape::solvelevel(context, state);
    const std::size_t count = surface.size();
    ASSERT_EQ(state.elevation_m.size(), count);
    ASSERT_EQ(state.pinned.size(), count);
    EXPECT_EQ(state.solved_m, state.elevation_m);
    EXPECT_EQ(state.bed_m, beds);
    EXPECT_EQ(state.surface_m, surface);
    std::size_t channelcells = 0, unpinned = 0;
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (state.pinned[cell] != 0)
            EXPECT_EQ(state.elevation_m[cell], surface[cell]) << cell;
        else
            unpinned++;
        if (state.channels.channel[cell] != 0 && state.sea[cell] == 0)
        {
            EXPECT_EQ(state.elevation_m[cell], state.bed_m[cell]) << cell;
            channelcells++;
        }
    }
    EXPECT_GT(channelcells, 0u);
    EXPECT_GT(unpinned, 0u);
    EXPECT_EQ(state.solveddischarge_m3_yr.size(), count);
}
