#pragma once

// Tributaries per level (docs/architecture/multiscale-plan.md, §3 "Kernels" and "Phase S3
// re-planned", step S3.3; the hook written in S3.0). Génevaux et al. 2013 grow a network level
// by level from the one before, with Horton-Strahler numbers as growth priorities (§6 "Learned in
// S3"); here every level keeps its parent's tributaries and grows one more generation from all of
// its channels, so the drainage density rises with the level. Planet-free.
//
// Tributary cells are channel cells (FineChannels) with globalcell -1: the links of the global
// network (rasterisechannels) always carry their upstream global cell, so the two are told apart
// without another field, and the troughs (which read a link's ice catchment) leave tributaries
// alone. "Parent cell": a block of `edge` x `edge` level cells, edge = context.factor in the
// refinement loop (the parent level's cell) and k in the single pass (the global cell). With an
// edge below 3 no length lies in [2, edge) and nothing is grown. The keys are
// RegionalParameters::tributaries (regional_types.hpp, landscape.regional.tributary_*).
//
// 1. Inherited (level >= 2, context.parent set). Every land tributary cell p of the parent inside
//    the parent's kept box whose parent chain (its directions through tributary cells) reaches a
//    land link cell of the parent inside that box is redrawn; the chain gives its depth (tributary
//    cells down to the link) and the link's bed. p's node is one of its children, hashed from p's
//    absolute key, floor(nodemargin x factor) cells or more inside p. Its segment is the
//    8-connected straight line from the node towards its downstream parent cell q's node (q a
//    tributary) or towards the nearest land link cell of the level among the children of p and
//    its 8 neighbours (q a link), ending at the first cell that is the sea, a link or q's node (the
//    receiver), or after the first cell that 8-neighbours a land link cell (that link receives
//    it: the lowest bed, then the smallest key). A segment that leaves the domain or p's 3 x 3
//    parent cells, starts in the sea or finds no link is not drawn. Segments are computed each on
//    its own; a cell on several belongs to the smallest depth, then the smallest parent key. An
//    inherited cell that then drains onto neither a channel nor the sea (q's segment not drawn)
//    ends there with direction 0, and a loop of inherited cells (crossing segments) is cut at its
//    cell of the smallest key (TributaryReport::inheritedends). Discharge and flows are p's,
//    capped at q's (parent) or the receiving link's.
// 2. Grown (every level). Candidates sit on every land channel cell J with a direction (links and
//    inherited cells), one on each side: the tributary leaves J at the junction angle (the mean
//    times 1 + anglespread x a hashed value in [-1, 1)) from the upstream direction (from J's main
//    upstream channel neighbour, the largest discharge, to J's downstream cell; J to its
//    downstream cell for a source), a straight ray from J's centre of target length
//    Lt = max(2, min(edge, Lh / c) x (lengthminfraction + (1 - lengthminfraction) u)) cells (u
//    hashed), where Lh is Hack's law L = C A^h closed over the strip a tributary drains, A = s L:
//    Lh = (C s^h)^(1 / (1 - h)) km, s = max(spacing, edge c) km (c the level's cell size). The
//    ray's cells are kept while their path length from J (cell steps, sqrt(2) per diagonal) is
//    below edge, and stop before a cell outside the domain, in the sea, on a channel, whose bed
//    surface lies below J's bed, or (from the second cell on) 8-neighbouring a channel; the first
//    cell may neighbour only channel cells within one cell of J. A candidate with a path of
//    length >= 2 is valid. So a tributary's length lies in [2 cells, one parent cell), it stays
//    inside J's parent cell and that cell's 8 neighbours, and it never touches another channel
//    except at J. Each parent cell and side is one slot, open with probability
//    min(1, edge c / spacing) (hashed per slot; 1 at every level run while edge c exceeds the
//    spacing: at most one tributary per side and parent cell edge), won by its valid candidate of
//    highest priority: links first, then J's Strahler order among the inherited cells inside J's
//    parent cell and its 8 neighbours (Génevaux's priority, bounded), then a hash. A slot's winner
//    is grown unless the path of a higher-priority winner comes within one cell of its own (one
//    round). spacing is tributaries.spacing_km, or inputs.valleyspacing_km of the junction's
//    global cell where that field is given and positive (S3.6).
// 3. Beds and discharge. No tributary bed stands above the bed surface it crosses
//    (state.bedsurface_m with the incision active, else state.surface_m: at or below the start
//    surface; LevelState::tributarybedsabove counts the exceptions, 0 by this rule). A grown
//    path's bed is b + w (max(r, b) - b): b J's bed, r the running minimum of the bed surface from
//    the tip down to the cell (at least b, since the ray stops before lower ground), w the cell's
//    path length from J over the tip's (1 at the tip): it rises from J to the surface at the tip
//    and never rises downstream. An inherited segment's bed is min(L, r), L the bed of the
//    parent's link its chain reaches and r the running minimum of the bed surface along the
//    segment from its node: deep like the link it joins and never above the ground (a cell where
//    the ground dips makes the next segment's start a step up). (Measured on run 9's window A at
//    k 64, level 3, border step max against 12.7 m without tributaries: beds on the surface's
//    running minimum everywhere 1,123 m, inherited beds capped at the parent's 244 m, inherited
//    beds at the receiver's bed 76 m.) A grown cell's discharge is the runoff of J's global cell
//    (inputs.runoff_mm_yr, 0 without it) over the Hack catchment of the path length above it
//    plus one cell, A = (l / C)^(1 / h) km2, capped at J's discharge; its flows are J's scaled by
//    the discharge ratio (0 where J's discharge is 0). A tributary may cross a global cell's
//    edge, so an outlet of a global cell may carry a tributary's discharge, not the cell's.
//
// Reach. A grown tributary's decision reads at most tributaryreach parent cells around its cells:
// its junction lies within one parent cell, the slot's other candidates within one more, a
// conflicting winner's path within one cell of its own and that winner's slot and rank within two
// more, each rank inside its junction's parent cell and its 8 neighbours. An inherited segment
// reads the parent's kept box (the same in every window where both hold it; its chains stay
// within the tree a level grew, under two global cells) and the level's links and sea inside p's
// 3 x 3 parent cells. The halo rules (regionalhalo, refinementhalo) add tributaryreach global
// cells with the tributaries on, at least as many parent cells at every level.
// Throws std::invalid_argument with tributaries enabled when spacing_km or hackcoefficient is not
// positive and finite, hackexponent not in (0, 1), junctionangle_deg not in (0, 180), anglespread
// not in [0, 1), nodemargin not in [0, 0.5) or lengthminfraction not in [0, 1].

#include <cstdint>
#include <vector>

namespace landscape
{
struct LevelContext;  // regional_pass.hpp
struct LevelState;
struct RegionalWindow;  // regional_grid.hpp
struct FineChannels;    // channel_network.hpp

// The reach of a grown tributary's decision in parent cells (header, "Reach").
constexpr int tributaryreach = 5;

// What one call did (domain counts).
struct TributaryReport
{
    int candidates = 0;        // valid candidates (both sides counted)
    int slotwinners = 0;       // the slots' winners
    int grown = 0;             // grown tributaries (state.tributaries)
    int growncells = 0;
    int inheritedparents = 0;  // the parent's tributary cells read (kept box)
    int deadparents = 0;       // of them, not drawn
    int inheritedcells = 0;
    int inheritedends = 0;     // inherited cells given direction 0 (a dangling end or a loop cut)
    int bedsabove = 0;         // tributary cells whose bed stands above state.surface_m (state.tributarybedsabove)
    int maxorder = 0;          // the largest Strahler order of a tributary cell after the call
};

// Called by channelslevel after the rasterisation and before the beds are pinned into
// state.surface_m. Adds tributary cells to state.channels (channel 1, direction, bed_m,
// discharge_m3s, flows, globalcell -1) and sets state.tributaries (the grown tributaries),
// state.tributarycells (the domain's tributary cells, inherited and grown) and
// state.tributarybedsabove. With tributaries.enabled 0 nothing is changed and the counts are 0.
void tributarieslevel(const LevelContext& context, LevelState& state);

// tributarieslevel with its report.
TributaryReport growtributaries(const LevelContext& context, LevelState& state);

// The Strahler order of every tributary cell (channel, globalcell < 0, not sea) on the window's
// domain grid: 1 at a tributary cell no tributary cell drains into, else the largest order among
// the tributary cells draining into it, plus one where two or more of them have that order; 0 on
// every other cell (the links are the trunks, outside the tributaries' orders).
std::vector<uint8_t> tributaryorders(const RegionalWindow& window, const FineChannels& channels,
                                     const std::vector<uint8_t>& sea);
}
