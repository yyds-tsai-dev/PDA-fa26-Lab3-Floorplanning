# E3: exact LP WL polish on a fixed topology (2026-10-29)

Gating experiment for deepened idea #1 in `adhd-solver-ideas.md`. The question: how much WL does an exact per-axis LP gain over packed positions when the H/V constraint graphs are fixed? Should we write a C++ network simplex?

Prototype: `tools/lp_polish_proto.py` (Python, scipy 1.18.1 HiGHS, networkx 3.7). Reproduce with:

```
python tools/lp_polish_proto.py --variants 8 --mcf 2 --out-dir <dir> --json <dir>/results.json
```

The full run takes 13 s for all 20 dev cases.

## Method

1. **Topologies.** Each of the 20 `cases/dev/*.in` starts from `<case>.construction.out`, using the construction's shapes. For every block pair:
   - positive y-overlap gives an H edge (left to right);
   - positive x-overlap gives a V edge (bottom to top);
   - disjoint on both axes (touching counts as disjoint) is a diagonal pair, resolved in one of two ways. Variant 0 uses the larger-gap axis. Variants 1–7 flip a seeded coin per pair, with seed `crc32("<case>/<v>")`.

   Every edge follows the construction's geometric order, so both graphs are DAGs. The construction satisfies every variant's constraints. The LP and packing run on the transitive reduction (`nx.transitive_reduction`).
2. **Packed placement.** ASAP longest path from 0 in topological order. Preplaced blocks are pinned at their input coordinates. A variant is infeasible if a predecessor ends after a preplaced block's start. RIGHT/TOP blocks are flushed to `W−w` / `H−h`; they are always graph sinks. LEFT/BOTTOM blocks are sources and land at 0. Every packing was checked with `./verifier/verifier`: gate FEASIBLE, full boundary and full MIB, and reported WL = computed WL.
   - **Attrition was 0/160.** Every variant is construction-consistent, so ASAP ≤ construction coordinates and the packing always fits. This gives 160 feasible packings (20 cases × 8).
3. **LP, per axis** (§Formulation): `linprog(method='highs-ds')`, which returns a vertex. All 160 × 2 solutions were integral, with 0 fractional coordinates, so no rounding was needed. Every LP placement passed the verifier with gate FEASIBLE, full boundary, full MIB, and verifier WL = the WL written to the file = (obj_x + obj_y)/2.
4. **MCF cross-check.** The dual min-cost circulation was solved with `nx.network_simplex` on variants 0 and 1 of every case (40 packings, 80 axis solves). The x placement was recovered from shortest-path potentials on the residual graph.
   - The optimum matched the LP exactly on 40/40 packings.
   - The recovered placement had the same WL and passed the verifier on 40/40.
5. **Extra: construction → LP.** The LP depends only on the graphs, not on the starting positions. The construction-as-is LP is therefore the variant-0 LP. It is compared against the construction WL, which equals the `WLTHRESHOLD` placeholder.

## Results

Gain = (WL_packed − WL_LP) / WL_packed. ρ = Spearman(WL_packed, WL_LP) across the 8 variants of a case. Times are scipy wall-clock medians per axis.

| case | n | #var | median gain % | gain range % | ρ | construction→LP gain % | LP time x / y (ms) |
|---|---|---|---|---|---|---|---|
| hidden1_s1001 | 80 | 8 | 4.16 | 3.57–5.47 | -0.33 | 4.08 | 8.3 / 7.3 |
| hidden1_s1002 | 80 | 8 | 2.75 | 2.36–3.84 | -0.17 | 2.99 | 7.4 / 7.2 |
| hidden1_s1003 | 80 | 8 | 3.76 | 3.07–4.54 | -0.60 | 3.21 | 8.0 / 7.9 |
| hidden1_s1004 | 80 | 8 | 3.63 | 3.39–4.96 | -0.21 | 3.60 | 8.4 / 7.5 |
| hidden2_s1001 | 120 | 8 | 3.26 | 2.98–4.75 | -0.38 | 2.83 | 10.9 / 10.7 |
| hidden2_s1002 | 120 | 8 | 3.00 | 2.79–4.25 | -0.69 | 2.57 | 10.7 / 10.5 |
| hidden2_s1003 | 120 | 8 | 2.00 | 1.73–3.90 | -0.52 | 3.23 | 11.5 / 11.0 |
| hidden2_s1004 | 120 | 8 | 2.93 | 2.81–3.74 | -0.07 | 3.25 | 11.5 / 10.7 |
| hidden3_s1001 | 160 | 8 | 2.09 | 2.01–2.89 | +0.29 | 1.95 | 14.1 / 13.4 |
| hidden3_s1002 | 160 | 8 | 2.06 | 1.84–2.66 | -0.07 | 2.38 | 13.6 / 12.6 |
| hidden3_s1003 | 160 | 8 | 2.14 | 2.07–3.30 | +0.02 | 2.52 | 13.2 / 13.3 |
| hidden3_s1004 | 160 | 8 | 2.39 | 2.26–3.23 | -0.17 | 2.65 | 14.8 / 12.9 |
| hidden4_s1001 | 200 | 8 | 1.70 | 1.59–2.10 | -0.12 | 2.12 | 18.4 / 16.6 |
| hidden4_s1002 | 200 | 8 | 1.77 | 1.68–2.85 | -0.33 | 2.13 | 18.4 / 17.7 |
| hidden4_s1003 | 200 | 8 | 1.87 | 1.78–2.49 | -0.86 | 1.90 | 16.8 / 16.5 |
| hidden4_s1004 | 200 | 8 | 2.63 | 2.48–3.95 | -0.21 | 2.56 | 16.1 / 16.2 |
| public_s1001 | 50 | 8 | 4.16 | 3.64–6.52 | -0.43 | 5.57 | 5.0 / 5.1 |
| public_s1002 | 50 | 8 | 3.56 | 2.99–4.50 | +0.07 | 4.06 | 5.2 / 5.1 |
| public_s1003 | 50 | 8 | 3.97 | 3.68–5.13 | +0.48 | 5.11 | 5.0 / 5.3 |
| public_s1004 | 50 | 8 | 3.60 | 3.13–5.05 | +0.05 | 4.28 | 4.9 / 5.3 |

### Overall

- **Packed → LP gain**, over all 160 packings: median **2.90 %**, p10 1.81 %, p90 4.22 %, min 1.59 %, max 6.52 %, mean 2.98 %.
  - The deterministic (larger-gap) variant alone has a median gain of 3.93 %.
  - The gain shrinks with size and tightness:
    - public about 3.6–4.2 %;
    - hidden1 about 3.6 %;
    - hidden2 about 3.0 %;
    - hidden3 about 2.1 %;
    - hidden4 about 1.8 %.
- **Construction → LP gain** (20 cases): median **2.91 %**, p10 2.10 %, p90 4.36 %, range 1.90–5.57 %. The LP beats the construction WL, i.e. the current `WLTHRESHOLD` placeholder, on every case.
- **Spearman ρ**: median **−0.19**, max +0.48. ρ < 0.8 on 20/20 cases. Across the variants of one case, packed WL carries essentially no information about LP WL. The pooled, per-case-normalised ρ is −0.02.
  - In 20/20 cases, variant 0 (larger-gap rule) had the best LP WL. Its packed WL ranked 8th of 8 (the worst) in 14/20 cases, with ranks `[8,8,6,8,8,8,8,6,8,6,8,8,5,8,8,8,8,5,3,8]`.
  - ρ is negative on 15/20 cases. The sign is systematic, not noise: packed WL tends to be anti-aligned with LP WL.
  - Hypothesis, not verified: the gap rule leaves fewer binding diagonal constraints. ASAP then compacts harder toward the origin, which hurts packed WL because terminals sit on all four sides. The same freedom is what the LP exploits.
  - **Caveat:** all variants share one construction and differ only in diagonal-pair choices. The spread is small: WL_LP varies by 0.4–1.6 % within a case and WL_packed by 0.2–1.6 %. This is a ranking test between near-identical topologies, which is the regime that late SA and SA move acceptance live in. It is not a test across unrelated topologies.
- **Problem sizes** (x axis, variant 0; y is similar):

  | n | LP vars (n + hinges) | LP rows | hinges | H edges full → TR | V edges full → TR |
  |---|---|---|---|---|---|
  | 50 | 223–237 | 474–510 | 173–187 | 550–657 → 114–136 | 568–675 → 117–126 |
  | 80 | 369–379 | 787–831 | 289–299 | 1552–1865 → 209–239 | 1295–1608 → 190–259 |
  | 120 | 522–534 | 1128–1235 | 402–414 | 2996–3997 → 324–413 | 3143–4144 → 343–368 |
  | 160 | 671–716 | 1540–1645 | 511–556 | 5882–7212 → 501–575 | 5508–6838 → 473–524 |
  | 200 | 861–884 | 1958–2096 | 661–684 | 8470–10862 → 626–745 | 9038–11430 → 623–681 |

  After transitive reduction, each graph has about 2.3–3.7 n edges, against 25–55 n before.
- **Solve times** (scipy HiGHS dual simplex, including model setup inside HiGHS): median 11.0 ms per axis, p90 16.9 ms, max 21.0 ms (n = 200). networkx `network_simplex`, in pure Python, takes 7–50 ms per axis. A C++ network simplex should be about 100× faster than networkx, roughly 0.1–0.5 ms per axis at n = 200.

## Formulation that matched (LP = MCF on 40/40)

Per axis, work in integer lower-left coordinates `x_i` with a root `r` (`x_r = 0`), in units of 2·WL.

**Hinges.** Take a connection of weight `c` between `a` and `b`, and let `d = x_a − x_b`.
- Block–block: `δ = w_a − w_b`.
- Block–terminal: `b := r` and `δ = w_a − 2·t_x`.

Then `2·c·|cx_a − cx_b| = c·|2d + δ|`, which becomes:
- **δ even:** one hinge, weight `2c`, offset `k = δ/2`: `2c·|d + k|`.
- **δ odd:** two hinges, weight `c` each, offsets `k = (δ−1)/2` and `k = (δ+1)/2`. This is exact for integer `d`, because `|u| + |u+1| = |2u+1|`.

All offsets are integers. Doubled coordinates are not needed. With doubled coordinates an odd `X` would give a half-integer `x`.

**Primal LP** (2·WL_axis = objective):

```
min  Σ_h wt_h · t_h
s.t. t_h ≥  (x_a − x_b + k_h),  t_h ≥ −(x_a − x_b + k_h)       per hinge
     x_j − x_i ≥ w_i                                           per TR edge i→j
     0 ≤ x_i ≤ L − w_i                                         (L = W or H)
     x_i = fixed_i       preplaced; LEFT/BOTTOM → 0; RIGHT/TOP → L − w_i
```

The hinge rows include `t_h`, so the LP matrix itself is not totally unimodular: for one hinge, columns `{x_a, t_h}` give the determinant `−2`. Integrality comes from duality instead:
- The LP is the dual of the min-cost circulation below, which has integer costs and capacities.
- Optimal potentials from shortest paths on its residual graph are therefore integral, and they are an optimal `x`.
- A C++ network simplex computes exactly these potentials, so integrality is free there.

All 320 HiGHS vertices were integral. That is empirical confirmation, not the proof.

**Dual min-cost circulation.** `nx.MultiDiGraph`, all node demands 0. Arc cost is the negated "profit"; "no cap" means infinite capacity.

| primal term | arc(s) | cost | capacity |
|---|---|---|---|
| hinge `wt·\|x_a − x_b + k\|` | `a → b` | `−k` | `wt` |
| | `b → a` | `+k` | `wt` |
| constraint edge `x_j − x_i ≥ w_i` | `i → j` | `−w_i` | ∞ |
| lower bound `x_i ≥ lo_i` | `r → i` | `−lo_i` | ∞ |
| upper bound `x_i ≤ hi_i` | `i → r` | `+hi_i` | ∞ |

- Lower and upper bounds come from the table above: `lo = 0` and `hi = L − w_i`, or `lo = hi = fixed_i` for fixed blocks.
- A terminal hinge uses `b = r`.

Optimum: `−flowCost = LP objective` (2·WL_axis). The arc pair `a→b` / `b→a` with capacity `wt` each is equivalent to one bidirectional flow `g ∈ [−wt, wt]` with profit `k·g`.

**Recovering x.** Build the residual graph:
- an arc `u→v` with cost `w` appears when `f < cap`;
- a reverse arc `v→u` with cost `−w` appears when `f > 0`.

Take Bellman-Ford distances `π` from `r` and set `x_i = −(π_i − π_r)`. This is an integral optimal primal, which the verifier confirmed on 40/40. An infeasible topology shows up as `NetworkXUnbounded`, a positive-profit cycle of infinite arcs. It did not occur here.

For C++, this is exactly LEMON-style network simplex input: n + 1 nodes, about `2·#hinges + #TR edges + 2n` arcs. At n = 200 that is about 1350 + 700 + 400 ≈ 2.5k arcs per axis.

## Caveats

- The baselines are generator-construction layouts, and packing them barely changes WL (packed ≈ construction ± 1 %). An SA incumbent that already optimises WL through packing may leave less on the table, but the effect can go either way:
  - SA's packed positions are ASAP-biased exactly as here;
  - the construction is a hand-anchored layout, not WL-optimised.
- Shapes were frozen. The soft-shape gradient step, the next idea in deepened #1, was not tested.
- Spearman was measured across small perturbations of one topology, not across SA-scale diversity.

## Recommendation: **GO** for a C++ network simplex

- **Gate criterion:** median gain ≥ 3 % **or** ρ < 0.8.
  - Median gain is 2.90 %, just under 3 %. The construction→LP gain is 2.91 %.
  - ρ < 0.8 holds on 20/20 cases (median −0.19). The criterion is met through the ρ branch.
- **The correlation result matters more than the gain.** Packed WL does not rank neighbouring topologies the way the exact geometry does. It tends to rank them in the opposite order: the LP-best topology had the worst packed WL in 14/20 cases. A post-pass LP alone is therefore not enough to recover the WL lost to the ranking error. The LP should become the SA cost evaluator late in the run (the two-clock design), at least for acceptance among the top candidates. That requires a fast solver: about 22 ms per evaluation with scipy versus an estimated sub-ms in C++.
- **Expected value as a post-pass:**
  - about 2–4 % WL at no feasibility risk, since boundary, MIB and preplaced are all preserved;
  - smallest on hidden4 (about 1.8 %). The tightest case is where T will be set, so this still matters at the margin.
- **The port is de-risked.**
  - The arc list above matched HiGHS exactly on 40/40 instances.
  - Integrality holds by construction (the odd/even hinge split).
  - The instances are small: at most 201 nodes and about 2.5k arcs per axis.
