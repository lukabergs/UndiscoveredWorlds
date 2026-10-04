/******************************************************************************
 *  plate-tectonics, a plate tectonics simulation library
 *  Copyright (C) 2012-2013 Lauri Viitanen
 *  Copyright (C) 2014-2015 Federico Tomassetti, Bret Curtis
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, see http://www.gnu.org/licenses/
 *****************************************************************************/

#include "plate.hpp"
#include "material.hpp"
#include "gtest/gtest.h"
#include "noise.hpp"
#include "simplexnoise.hpp"

void initializeHeightmapWithNoise(long seed, float *heightmap, const WorldDimension& wd)
{
    createNoise(heightmap, wd, SimpleRandom(seed), true);
    for (uint32_t i=0; i<wd.getArea(); i++) {
        if (heightmap[i]<0.0f) {
            heightmap[i] *= -1.0f;
        }
    }
}

TEST(SimpleRandom, NextRepeatability)
{
    SimpleRandom sr1(1);
    EXPECT_EQ(81414, sr1.next());
    EXPECT_EQ(1328228615, sr1.next());
    EXPECT_EQ(3215746516, sr1.next());

    SimpleRandom sr999(999);
    EXPECT_EQ(69012276, sr999.next());
    EXPECT_EQ(3490172125, sr999.next());
    EXPECT_EQ(3364058674, sr999.next());
}

TEST(Noise, SimplexRawNoiseRepeatability)
{
    EXPECT_FLOAT_EQ(-0.12851511f, raw_noise_4d(0.3f, 0.78f, 1.677f, 0.99f));
    EXPECT_FLOAT_EQ(-0.83697641f, raw_noise_4d(-0.3f, 0.78f, 1.677f, 0.99f));
    EXPECT_FLOAT_EQ(-0.5346418f, raw_noise_4d(7339.3f, 0.78f, 1.677f, 0.99f));
    EXPECT_FLOAT_EQ(0.089452535f, raw_noise_4d(0.3f, 70.78f, 1.677f, 0.0009f));
    EXPECT_FLOAT_EQ(-0.063593753f, raw_noise_4d(0.3f, 500.78f, 1.677f, 500.99f));
}

TEST(Noise, SimplexNoiseRepeatability)
{
    WorldDimension wd(233, 111);
    float *heightmap = new float[wd.getArea()];
    initializeHeightmapWithNoise(123, heightmap, wd);

    EXPECT_FLOAT_EQ(0.50098729f, heightmap[0]);
    EXPECT_FLOAT_EQ(0.39222634f, heightmap[1000]);
    EXPECT_FLOAT_EQ(0.51659518f, heightmap[2000]);
    EXPECT_FLOAT_EQ(0.5479334f, heightmap[5000]);
    EXPECT_FLOAT_EQ(0.59222502f, heightmap[8000]);
    EXPECT_FLOAT_EQ(0.36362505f, heightmap[11000]);
    EXPECT_FLOAT_EQ(0.57599854f, heightmap[13000]);

    delete[] heightmap;
}

TEST(CreatePlate, SquareDoesNotExplode)
{
    float *heightmap = new float[40000]; // 200 x 200
    initializeHeightmapWithNoise(678, heightmap, WorldDimension(200, 200));
    plate p = plate(123, heightmap, 100, 3, 50, 23, 18, WorldDimension(200, 200));
}

TEST(CreatePlate, NotSquareDoesNotExplode)
{
    float *heightmap = new float[80000]; // 200 x 400
    initializeHeightmapWithNoise(678, heightmap, WorldDimension(200, 400));
    plate p = plate(123, heightmap, 100, 3, 50, 23, 18, WorldDimension(200, 400));
}

TEST(Plate, ClassifiesOceanicAndContinentalTypes)
{
    constexpr uint32_t width = 8;
    constexpr uint32_t height = 8;
    float *oceanicMap = new float[width * height];
    float *continentalMap = new float[width * height];
    std::vector<uint8_t> basalt(width * height,
        tectonics::material::to_index(tectonics::material::Type::Basalt));

    for (uint32_t i = 0; i < width * height; ++i) {
        oceanicMap[i] = 0.18f;
        continentalMap[i] = oceanicMap[i];
    }

    // Material composition determines plate type independently of surface height.
    plate oceanic = plate(11, oceanicMap, width, height, 0, 0, 1,
                          WorldDimension(64, 64), basalt.data());
    plate continental =
        plate(22, continentalMap, width, height, 16, 16, 1, WorldDimension(64, 64));

    EXPECT_TRUE(oceanic.isOceanicPlate());
    EXPECT_TRUE(continental.isContinentalPlate());
    EXPECT_GT(continental.buoyancy(), oceanic.buoyancy());
}

TEST(Plate, OceanicBuoyancyTracksRelativeLightness)
{
    constexpr uint32_t width = 8;
    constexpr uint32_t height = 8;
    float *deepOceanMap = new float[width * height];
    float *lightOceanMap = new float[width * height];
    std::vector<uint8_t> basalt(width * height,
        tectonics::material::to_index(tectonics::material::Type::Basalt));

    for (uint32_t i = 0; i < width * height; ++i) {
        deepOceanMap[i] = 0.08f;
        lightOceanMap[i] = 0.22f;
    }

    plate deepOcean =
        plate(33, deepOceanMap, width, height, 0, 0, 1, WorldDimension(64, 64), basalt.data());
    plate lightOcean =
        plate(44, lightOceanMap, width, height, 20, 20, 1, WorldDimension(64, 64), basalt.data());

    EXPECT_TRUE(deepOcean.isOceanicPlate());
    EXPECT_TRUE(lightOcean.isOceanicPlate());
    EXPECT_GT(lightOcean.buoyancy(), deepOcean.buoyancy());
}

// TODO test also when plate is large as world
