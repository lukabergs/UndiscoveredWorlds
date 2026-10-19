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

#include "simulations/geology/plate_tectonics/tectonic_api.hpp"
#include "gtest/gtest.h"
#include <cstdlib>
#include "simulations/geology/plate_tectonics/simplerandom.hpp"

using namespace std;

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

TEST(Randomness, SequenceDoubles)
{
    long seed = 3;
    SimpleRandom randsource(seed);

    EXPECT_FLOAT_EQ(5.1118433e-05f, randsource.next_float());
    EXPECT_FLOAT_EQ(0.53070194f, randsource.next_float());
    EXPECT_FLOAT_EQ(0.053402752f, randsource.next_float());
}

TEST(Randomness, Maximum)
{
    long seed = 3;
    SimpleRandom randsource(seed);

    EXPECT_FLOAT_EQ(4.2949673e+09f, (float)randsource.maximum());
}
