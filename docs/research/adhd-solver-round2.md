# ADHD round 2: core algorithm beyond SA (2026-10-29)

Problem given to every branch: single-threaded fixed-outline floorplanner, 110 s, n ≤ 200, 6–8% whitespace, interior preplaced blocks, boundary and MIB constraints, T = median of 10 reference seeds with 0 points at 1.05T. Measured state: SA from its own seeds ends at 1.2–1.6× construction WL with per-seed max/median up to 1.28; SA from the construction topology reaches 0.92–0.95×; analytical placement is 1.03–1.07× but overlapping, and two direct legalizers failed (docs/research/legalizer-nogo.md).

Frames: logistics, speedrunner, markets, infinite budget, 10-year-old. Scores are [novelty viability fit], weighted 0.35/0.40/0.25.

## Clusters

**A. Stay legal, improve by local rip-up (LNS)**
- Rip out 1–6 blocks and re-insert each at the best legal position [N6 V8 F9]
- Bankruptcy / re-listing of blocks whose nets are too long [N6 V8 F8] (same idea)
- Ejection chains when the target pocket is short by δ [N7 V6 F8]
- Whitespace "hole" agent: every move is a sliding-puzzle swap with the hole [N8 V3 F5]: trap, one hole cannot move heterogeneous blocks at 6% whitespace

**B. Build legally from the start (constructive)**
- Rocks then water: preplaced, then rigid blocks, then soft blocks reshaped to fit pockets [N7 V5 F8]
- Legal partial placement plus an unplaced queue whose area is the only defect [N7 V6 F8]
- Balloon growing from the analytical spots [N8 V4 F7]
- Strip scheduling (Abacus with preplaced blocks as pre-booked slots) [N6 V5 F6]
- Pixel land-grab on a grid, then snap regions to rectangles [N7 V3 F6]: trap, the rectangle snap is the whole difficulty

**C. Fight seed variance**
- Population of elites + geometric window crossover + per-block median consensus [N7 V6 F8]
- Majority vote on SP relations [N6 V3 F7]: trap, pairwise votes can form cycles and far-pair votes are noise
- Fund-of-funds successive halving [N5 V7 F7]: trap, measured Spearman of early vs final rank is about −0.15
- Fixed internal RNG so every seed returns the same output [N6 V9 F3]: trap, it only hides the reference's variance; students' variance is unchanged and T becomes one run's luck

**D. Partition by the preplaced blocks (rooms)**
- Rooms / catchments / clearing house / area trucking between rooms [N6 V4 F7]: trap, the obstacle-aware bisection experiment showed about 60 of 90 rooms cannot fit their blocks
- Overlap-debt netting by min-cost flow [N7 V4 F7]: trap, same fragmentation failure as the constraint-graph legalizer

**E. Shape and slack tricks**
- Shrink soft blocks to 0.99·A during search [N5 V10 F6]: check whether shape.h already allows it
- Co-quantize soft shapes to shared strip heights to avoid slivers [N7 V6 F7]
- Neighbour pairs trade aspect-ratio slack [N4 V9 F7]
- Filler "air" blocks so packing is exact [N7 V5 F6]

**F. Other**
- Learned move priors from offline training [N8 V2 F6]: trap, no time or data
- Rectangular dual of the analytical adjacency + convex sizing [N8 V3 F7]: trap, preplaced faces and whitespace break the dual
- Read the SP straight off the overlapping analytical placement [N4 V8 F7]: trap, already the current seed path (1.6–1.9× at first legal)
- Pin the boundary ring first, search only the interior [N5 V6 F6]
- Multilevel: cluster tightly connected soft blocks into super-blocks, place small blocks last [N6 V6 F7]

## Shortlist
1. ★ **LNS with an exact insertion oracle** (A). Rip a block out of the SP, compute the legal x/y window of every one of the (n+1)² insertion cells with four ASAP/ALAP quadrant sweeps (about 100 µs), score legal cells by the block's net-weighted median, evaluate the top ~8 with a real pack + warm LP. About 1 ms per move. Every state is legal, so no time is lost regaining legality. Can also run inside SA as an always-legal move.
2. **Legal-first constructive seed** (B). Free space kept as maximal empty rectangles; place rocks, MIB groups, rigid and boundary blocks, then pour soft blocks with widths snapped to the pocket; queued area is the only defect; repair evicts 2–8 soft neighbours in a window and repacks it. Gives a legal start that keeps the analytical topology.
3. **Variance attack, gated by a measurement** (C). First decompose the variance: fork 10 continuations from one run's first-legal state. If they stay within about 1.03, variance comes from early topology and an elite pool + geometric crossover + median-centroid consensus is justified. If they spread about 1.2×, attack the tail instead (LNS, best-of-k tail branches).

## Deepened

### 1. LNS with insertion oracle
- **Risk:** at 6–8% whitespace the neighbourhood may collapse to "put it back where it was". If fewer than about 10% of moves find a feasible, improving, non-original cell, it stagnates like SA. The oracle must include paths through preplaced pins and boundary source/sink terminals, or it will accept illegal cells.
- **First step:** build `InsertionOracle(sp_without_b, outline, pins)` and test it against brute-force pack-and-check on random legal SPs (n = 30–60, with pins and boundary blocks), 100% agreement required. Then run it on finished SA results at n = 200 and log the fraction of blocks with a feasible non-original cell and the median proxy WL gain.
- **Children:** ejection chains to depth 3; oracle-sampled moves inside SA; LNS as a legalizer for the analytical placement; reshape soft blocks along critical paths after each move; path relinking between elites.

### 2. Legal-first constructive seed
- **Risk:** the tail. The last 10–15% of soft blocks meet whatever pockets the greedy left; area-feasible pockets can be shape-infeasible; repairs can cascade. Dead-space-first scoring can pull blocks off their analytical targets.
- **First step:** maximal-empty-rectangle structure with incremental split and dominance pruning, plus a WL-blind corner-occupying greedy. Go/no-go: queued area under about 1% on most hidden-like cases.
- **Children:** beam search over the placement order; pocket reservation by min-cost transport before pouring; window rip-up-and-repour as an LNS operator; extract an SP from the result and hand it to SA + LP polish; multi-front pouring from edges and rocks.

### 3. Population + consensus (gated)
- **Risk:** if variance comes from the stochastic tail rather than early topology choices, a pool cannot average it away; crossover-child repair may cost as much as reaching legality from scratch.
- **First step:** the variance-decomposition experiment above, plus timing the repair of 20 window-crossover children.
- **Children:** path relinking (about 80 ms per path); median centroids as soft anchors for the analytical seeder; freeze relations shared by ≥90% of elites; steer moves to high-disagreement blocks; restart from a pool child on stagnation.

## Provocation
The construction topology gets SA to 0.92–0.95×. The generator builds that topology by recursive dissection with zero whitespace and then shrinks blocks. Could the solver do the same: find a zero-whitespace slicing/pinwheel dissection of the outline around the preplaced blocks that matches the analytical placement, then shrink?

## Results

- **Variance source (item 3), measured 2026-10-29** on dev_tight hidden2/3/4 at ws 0.06–0.07, Phase C binary, 110 s. Ten continuations forked from one run's first-legal state spread as widely as ten independent runs (max/median 1.41–1.53 vs 1.33–1.43). Continuations forked at 50 s spread only 1.01–1.14. So seed spread comes from the late stochastic tail, not from early topology. The population/consensus plan is dropped; the tail is attacked instead (LNS, best-of-k tail branches).
- **Restart portfolio**, dev_tight × 6 portfolios per case, 110 s total each: 1×110 s vs 3×36 s vs 2×55 s best-of. Median WL/construction 1.41 / 1.44 / 1.40; max/median 1.26 / 1.15 / 1.14; p90/median 1.17 / 1.11 / 1.10. Splitting the budget cuts spread but barely moves the median. Scripts and data: scratchpad `port/`, `vardec/`.
- **Legal-first constructive (item 2): no-go.** Maximal-empty-rectangle greedy (`solver/construct.h`, `test_construct.cpp`) leaves 4–6% of block area unplaced at 6–7% whitespace even with the specified repair (an extended 20 s repair still leaves 1.2–1.3%). Complete outputs are 1.96–3.51× construction WL, worse than the current seed path, because packing needs wall-hugging scoring that pulls blocks off their targets. Free space fragments into pockets too small for the large soft tail and slivers just under the soft blocks' min side.
- **LNS insertion oracle (item 1): no-go.** `solver/oracle.h` is exact (100% agreement with brute-force pack-and-check over 7.6M cells at n = 30–60 and 8M cells at n = 200) and takes 130–230 µs per block at n = 200. But on finished hidden4 outputs at 6–8% whitespace only 4–11% of movable blocks have an improving relocation, the median gain is 0.001–0.08% of WL, and the sum of all per-block best gains is 0.03–1.1%. Converged SA outputs are already local optima under single-block relocation; the 1.2–1.6× gap and the seed spread are basin-level, not local.

## Conclusion
All three shortlisted directions, plus the two earlier legalizer attempts, fail with measured blockers that agree: at 6–7% whitespace with interior preplaced blocks, the hard part is the global basin the anneal ends in, and the spread is created in the post-feasible anneal. Phase E therefore tunes the post-feasible schedule at 110 s and adds an internal 2-way restart instead of opening new algorithm fronts.
