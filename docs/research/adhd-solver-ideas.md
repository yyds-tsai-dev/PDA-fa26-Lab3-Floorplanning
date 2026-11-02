# ADHD ideation: reference solver architecture (2026-10-29)

Five isolated divergent branches (hardware engineer, competitor, inversion, remove-the-assumption, biology), 6 ideas each. Scores use novelty/viability/fit, weighted 0.35/0.40/0.25.

## Clusters

**A. Topology search + exact geometry (LP / min-cost flow)**
- Two-clock: a cheap packer runs in the SA loop, and an LP/MCF on frozen constraint graphs runs periodically and as a post-pass [N7 V8 F9]
- Post-legal WL compaction LP with the topology fixed [N5 V9 F10]
- Soft shapes as an inner deterministic solve (slack distribution / Lagrangian sizing) [N7 V6 F9]
- Slicing tree + Stockmeyer shape curves [N6 V6 F6]: trap, because interior preplaced blocks break slicing

**B. Constraints by construction**
- Shape interval per soft block; MIB shares one w variable [N4 V10 F9]
- Boundary as an SP source/sink invariant, with sinks flushed to W−w / H−h [N7 V7 F9]
- Boundary perimeter slots matched to terminal pull [N7 V6 F8]

**C. Throughput / anytime engineering**
- Incremental CSR WL delta [N3 V10 F9]
- Always-legal incumbent from second 0, plus a watchdog [N3 V10 F8]
- Profile-guided budget scheduler (measure moves/s, size the schedule) [N5 V9 F7]
- Embedding-seeded successive-halving portfolio [N7 V8 F8]
- Adaptive move-type selection (clonal) [N6 V8 F7]
- Restarts with an early-abort predictor [N5 V8 F8]

**D. Preplaced as skeleton**
- Partition the outline into channels around preplaced blocks and solve the sub-floorplans separately [N7 V5 F8]: trap, because fragmentation and cross-region migration hurt at 8% whitespace
- Tissue growth: soft blocks inflate around preplaced seeds [N8 V5 F7]

**E. Representation-free**
- Direct x/y local search with overlap repair [N6 V6 F7]: risky at 8% whitespace
- Differential-adhesion cell sorting [N8 V4 F6], quorum-sensing legalization [N7 V4 F5], mosaic-grid assignment [N7 V4 F6]: traps (heterogeneous sizes, no legality guarantee)
- Hebbian net pruning [N6 V5 F4]: trap, because n ≤ 200 with ~700 connections makes WL delta cheap already

## Deepened

### 1. Two-clock LP geometry
- **Freeze.** Build H/V constraint graphs from the packed geometry: an H edge if y-projections overlap, a V edge if x-projections overlap, and the larger-slack direction for diagonal pairs. Transitive reduction with 200-bit reach bitsets leaves about 2–3n edges per graph.
- **Integer handling.** With d = x_a − x_b and δ = w_a − w_b, each connection's cost is 2|cx_a − cx_b| = |2d + δ|. Even δ gives a hinge with weight 2c and breakpoint −δ/2. Odd δ splits into two integer-breakpoint hinges. So the LP runs on integer lower-left coordinates, with no doubled coordinates needed.
- **Solve.** Use the dual min-cost circulation: nodes are the blocks plus a root, each hinge is an arc with capacity 2c, constraint and outline edges are infinite-capacity arcs, and preplaced equalities are two root arcs. Boundary equalities are hinges with capacity M, so they cannot cause infeasibility. Solve with a hand-written network simplex (port of LEMON's block-search pivot, ~300 lines). Expect ~0.5–3 ms per dimension at n = 200.
- **Soft shapes.** Use arc flows as ∂WL/∂w gradients, move each block along its constant-area curve (summed over MIB groups), check with Bellman-Ford, halve the step if infeasible, and do 2–3 rounds.
- **Risk.** At ≤ 8% whitespace the slack is small. The gain may be under 3%, or the packed WL may rank topologies poorly against LP WL.
- **First step.** Prototype the LP in Python on about 50 packings, check with ./verifier, and measure % gain and rank correlation before writing C++.
- **Children.**
  - Warm-started network simplex as the cost function in late SA.
  - Dual-guided H↔V edge flips.
  - Soft-outline LP legalizer that minimizes overflow·M + WL.
  - Joint convex sizing (w and positions) via Lagrangian iteration.

### 2. Embedding-seeded successive-halving portfolio
- **Embed.** Quadratic embedding with terminals and preplaced blocks as anchors, plus IRLS toward L1. SimPL-style rough-legalization anchors (area bisection with preplaced blockages) prevent collapse. All of this takes only a few ms.
- **Convert.** Γ− sorted by x+y and Γ+ by x−y, with diversity from snapshots, jitter, and axis rotation.
- **Why SP over B*-tree.** SP can pin preplaced blocks at their coordinates; contour packing cannot hold an interior (x, y).
- **Race.** K = 32/16/8 runs at n = 50/100/200 under successive halving. Temperature follows global round progress, 60% of the time goes to halving and 40% to the final winner, ranking is lexicographic (gate, boundary, MIB, WL), and a global incumbent is kept throughout.
- **Validation.** A 2×2 ablation (seed × schedule) on the dev set with a paired Wilcoxon test.
- **Risk.** Depth versus diversity at the tightest whitespace. Check the Spearman correlation between round-1 rank and final rank.
- **Children.**
  - Exact MCF WL compaction (same as #1).
  - Legalize the embedding directly (Tetris/Abacus with soft reshaping), then extract the SP from that legal placement, which gives an always-legal incumbent in under 1 s.
  - MIB as a contiguous SP sub-run.
  - Consensus-anchored re-embedding for the dropped runs.
  - Hyperband-style adaptive brackets.

### 3. Constraints by construction + deterministic shaping
- **Soft shapes.** No table is needed: the legal widths form one interval, w ∈ [ceil(sqrt(0.99A/3)), floor(sqrt(3.03A))], with hmin(w) and hmax(w) given in closed form. The interval holds 117 widths at A = 1e4 and 1166 at A = 1e6. MIB members share one w.
- **Boundary in SP.** LEFT means an HCG source and RIGHT an HCG sink; TOP/BOTTOM are the same in the VCG. Sinks are flushed to W−w after packing, which is always overlap-safe. A violation therefore cannot be written down; it turns into outline overflow. The generator construction satisfies this invariant.
- **Shaping cadence.** Per move, try about 8 candidate w for the perturbed block and its critical neighbours. Periodically, run an ALAP/ASAP critical-path slack pass. Keep about 10% random shape moves.
- **Risk.** At hidden4 the invariant plus pinned preplaced blocks shrinks the reachable space. Compare time-to-feasible against a penalty-only baseline.
- **Children.**
  - Invariant-preserving move windows.
  - Final LP WL squeeze.
  - Use hmax to close 1–3 unit edge gaps.
  - Corner-first geometry-consistent initial SP.
  - Lagrangian sizing in closed form, w* = sqrt(A·λv/λh).

## Converged shortlist
1. **SP-based SA** with preplaced blocks pinned, boundary as a source/sink invariant with flush, and MIB sharing one w. The search space only contains layouts that satisfy MIB/boundary.
2. ★ **Exact MCF/network-simplex WL compaction on the frozen topology**, used as a post-pass and on every new best (later possibly in the late-SA loop). It is the non-obvious but viable pick, and both deepen #1 and #2 named it the largest single WL gain.
3. **Embedding seed + legal incumbent + throughput-sized schedule.** Turn on successive halving only if the Spearman check passes.
4. **Incremental WL delta, adaptive move selection, and Lagrangian soft sizing** as engineering multipliers.

Traps: preplaced channel partitioning, slicing/Stockmeyer, representation-free cell sorting / quorum / mosaic, net pruning.
