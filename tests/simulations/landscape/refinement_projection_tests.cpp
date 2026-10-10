// Tests of landscape/refinement_projection.hpp (docs/architecture/multiscale-plan.md,
// section 3, "Tests"), suite Landscape: Landscape.RefinementProjection*. S0.4 shows that
// a child and a parent window built by hand (refinement_test_windows.hpp) find each other's
// cells, which is all the projection needs from refinement_grid.hpp; S2.2 adds exactness,
// the weights, the border step and the resample.

#include <gtest/gtest.h>

#include "landscape_fixture.hpp"
#include "refinement_test_windows.hpp"
#include "simulations/landscape/refinement_grid.hpp"
#include "simulations/landscape/refinement_projection.hpp"

#include <cstddef>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

// A 16-column world at factor 2: the child at level 2 (64 x 32 cells) crosses the date
// line; the parent at level 1 (32 x 16 cells) keeps every parent cell of the child's
// domain. The parent's field holds each cell's absolute key, so the value a child reads
// names the cell it reached.
TEST_F(Landscape, RefinementProjectionHandBuiltWindowsFindTheirParents)
{
    const landscape::LevelWindow child = handlevelwindow(16, 2, 2, 1, 60, 8, 8, 4);
    const landscape::LevelWindow parent = handlevelwindow(16, 2, 1, 1, 28, 8, 2, 6);
    ASSERT_EQ(child.domaincolumn0, 58);
    ASSERT_EQ(child.domaincolumns, 12);
    ASSERT_EQ(child.domainrow0, 6);
    ASSERT_EQ(child.domainrows, 8);
    ASSERT_EQ(parent.domaincolumn0, 26);
    ASSERT_EQ(parent.domainrows, 10);

    std::vector<float> parentfield(parent.domainsize());
    for (int y = 0; y < parent.domainrows; y++)
        for (int x = 0; x < parent.domaincolumns; x++)
            parentfield[parent.index(x, y)] = static_cast<float>(parent.absolutekey(x, y));

    std::map<std::pair<std::int64_t, std::int64_t>, int> children;
    for (int y = 0; y < child.domainrows; y++)
        for (int x = 0; x < child.domaincolumns; x++)
        {
            const std::int64_t column = child.parentcolumn(x), row = child.parentrow(y);
            // The domain's p-th parent column is the parent of its cells factor p ..
            // factor p + factor - 1.
            EXPECT_EQ(column, landscape::wraplevel(child.domaincolumn0 / child.factor + x / child.factor,
                                                   parent.levelwidth()));
            EXPECT_EQ(row, child.domainrow0 / child.factor + y / child.factor);
            // Every parent cell of the child's domain lies in the parent's kept box.
            ASSERT_TRUE(parent.containscolumn(column));
            ASSERT_TRUE(parent.containsrow(row));
            const int px = parent.localcolumn(column);
            const int py = static_cast<int>(parent.localrow(row));
            EXPECT_GE(px, parent.keptcolumnoffset());
            EXPECT_LT(px, parent.keptcolumnoffset() + parent.columns);
            EXPECT_GE(py, parent.keptrowoffset());
            EXPECT_LT(py, parent.keptrowoffset() + parent.rows);
            EXPECT_EQ(parentfield[parent.index(px, py)], static_cast<float>(landscape::levelkey(column, row, 32)));
            EXPECT_EQ(child.globalcolumn(x), parent.globalcolumn(px));
            EXPECT_EQ(child.globalrow(y), parent.globalrow(py));
            children[{column, row}]++;
        }
    // Six parent columns (29 .. 31 and 0 .. 2) by four rows, each with factor^2 children.
    EXPECT_EQ(children.size(), static_cast<std::size_t>(child.domainparentcolumns() * child.domainparentrows()));
    EXPECT_EQ(children.size(), 24u);
    for (const auto& entry : children)
        EXPECT_EQ(entry.second, 4);
    EXPECT_EQ(children.count({29, 3}), 1u);
    EXPECT_EQ(children.count({2, 6}), 1u);
}
