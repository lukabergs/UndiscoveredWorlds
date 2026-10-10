// Tests of landscape/regional_tributaries.hpp (multiscale plan S3.0: the hook, inert; S3.3 grows
// the tributaries), suite Landscape: TEST_F(Landscape, RegionalTributaries<What>) on the
// synthetic world of regional_test_world.hpp, the level's steps up to the channels run through
// the hook interface (regional_pass.hpp).

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_tributaries.hpp"
#include "simulations/landscape/regional_valleys.hpp"

TEST_F(Landscape, RegionalTributariesOffLeaveTheChannelsUnchanged)
{
    // tributaries.enabled 0 (the default): the step changes no channel field and counts none.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalParameters parameters;
    ASSERT_EQ(parameters.tributaries.enabled, 0);
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
    state.channels = landscape::rasterisechannels(window, inputs.network, parameters.channels, inputs.sealevel,
                                                  state.surface_m, state.sea);
    ASSERT_GT(std::count(state.channels.channel.begin(), state.channels.channel.end(), uint8_t{1}), 0);
    const landscape::FineChannels before = state.channels;
    const std::vector<float> surface = state.surface_m;
    state.tributaries = 7;
    state.tributarycells = 9;

    landscape::tributarieslevel(context, state);
    EXPECT_EQ(state.channels.channel, before.channel);
    EXPECT_EQ(state.channels.direction, before.direction);
    EXPECT_EQ(state.channels.bed_m, before.bed_m);
    EXPECT_EQ(state.channels.discharge_m3s, before.discharge_m3s);
    EXPECT_EQ(state.channels.januaryflow, before.januaryflow);
    EXPECT_EQ(state.channels.julyflow, before.julyflow);
    EXPECT_EQ(state.channels.globalcell, before.globalcell);
    EXPECT_EQ(state.surface_m, surface);
    EXPECT_EQ(state.tributaries, 0);
    EXPECT_EQ(state.tributarycells, 0);
}
