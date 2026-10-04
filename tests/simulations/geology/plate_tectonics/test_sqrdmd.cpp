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

#include "sqrdmd.hpp"
#include <cstdio>
#include "gtest/gtest.h"

#include <cmath>
#include <vector>

TEST(SquareDiamond, SameSeedGivesTheSameFiniteMap)
{
    constexpr int side = 17;
    std::vector<float> first(side * side, 0.0f), second(first), other(first);
    ASSERT_EQ(sqrdmd(1, first.data(), side, 0.5f), 0);
    ASSERT_EQ(sqrdmd(1, second.data(), side, 0.5f), 0);
    ASSERT_EQ(sqrdmd(2, other.data(), side, 0.5f), 0);
    EXPECT_EQ(first, second);
    EXPECT_NE(first, other);
    for (float value : first)
        EXPECT_TRUE(std::isfinite(value));
}
