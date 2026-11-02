// fp.cpp -- TA reference solver for Lab3 fixed-outline floorplanning (Phase G).
//
//   fp <input> <output> [--time SEC=24] [--seed N=1] [--mode weak|strong]
//      [--init-from PLACEMENT.out] [--pack-only] [--set key=value]...
//
// --mode weak is the frozen Phase A solver (the calibration yardstick,
// docs/design.md §9): plain SA from a random SP, Solver::run(). Nothing below
// marked "strong" runs in weak mode, and weak mode's RNG stream is unchanged.
//
// Phase A (both modes): sequence-pair simulated annealing with constraints by construction
// (docs/research/adhd-solver-ideas.md, deepened #3).
//   * Representation: SP (Γ+, Γ-) over all blocks, including preplaced ones.
//     Convention in extract.h: a<b in both => a left-of b; a after b in Γ+ and
//     before b in Γ- => a below b.
//   * Packing: longest path per axis, computed in O(n log n) with a prefix-max
//     Fenwick tree over Γ- positions (FAST-SP style). A preplaced block p is
//     pinned at (x0, y0); its successors start at >= x0 + w; if its
//     predecessors end past x0 the excess is the preplaced violation.
//   * Soft shapes: one index into the block's legal-width list (shape.h),
//     h = hmin(w). An MIB group shares one index (identical A => identical
//     list => identical (w, h)).
//   * Boundary invariant: LEFT blocks are HCG sources, RIGHT HCG sinks, BOTTOM
//     VCG sources, TOP VCG sinks (corners both). Sinks are flushed after
//     packing (x := W - w / y := H - h) when that moves them right/up. Moves
//     that break the invariant, or that give two preplaced blocks a relation
//     contradicting their fixed geometry, are rejected.
//   * Overflow (documented choice): per-block sum
//       ovf = Σ_b max(0, x_b + w_b - W) + max(0, y_b + h_b - H),
//     which still has a gradient when the bounding box is set by one block.
//   * SA cost: wl2 / wl_norm + λ (ovf + viol) / (W + H), λ adaptive (raised
//     while the walk is mostly infeasible, lowered slowly while mostly
//     feasible); exponential temperature schedule over the time budget.
//   * A state is gate-feasible iff ovf == 0 and viol == 0; a candidate
//     incumbent is re-checked with a mirror of verifier gate checks 3-8
//     before it is accepted.
//
// Phase C (--mode strong, the default), Solver::run_strong():
//   * Forbidden relations: a block may never be left-of / right-of / below /
//     above a preplaced block p when p's gap to that outline edge is smaller
//     than the block's narrowest legal width / height (build_forbidden). No
//     move may add such a relation (the count may only fall; SP moves are
//     rejected otherwise). This removes the sliver plateaus of FloorSet cases.
//   * Analytical seeds (GlobalPlacer): quadratic placement + IRLS toward L1 +
//     SimPL-style spreading by recursive area bisection with preplaced area
//     excluded; the SP is extracted (extract.h) from a virtual placement of
//     grown rectangles around the bisection-leaf centres, ties broken by
//     x - y / x + y, so the boundary invariant holds by construction. Soft
//     shapes follow the leaf aspect ratio, then a greedy shape descent. Eight
//     seeds, ranked by packed penalty + WL; the best starts SA at 0.05 x the
//     probed temperature.
//   * Plateaus: s_plateau_k * n moves without a better penalty in the current
//     attempt. Near-feasible (best penalty <= 0.005 (W + H)): a soft-shape
//     descent and a short reheat from the attempt's best state (until 3 in a
//     row bring no better penalty). Otherwise restart from the next
//     seed (lambda reset, T = 0.3 x the re-probed T), then from random SPs
//     once the pool is used up.
//   * Post-feasible phase as Phase A.
//
// Phase D (strong only): exact WL polish (solver/mcf, min-cost-flow dual of
// the fixed-topology WL LP), Solver::on_incumbent().
//   * The packed incumbent (lowest packed WL; it drives SA bookkeeping) is
//     kept apart from the output (best_rects), which may be a polished
//     placement. The SA state is never changed by a polish.
//   * Each new packed incumbent is polished, at most once per s_pol_ms ms
//     and while polishing has used <= s_pol_budget of the elapsed time,
//     on its SP constraint graphs (sp_graphs: transitively reduced left-of /
//     below relations) and on its geometric graphs (mcf::polish); a skipped
//     one is polished later. The final incumbent is always polished, the
//     output is re-polished on its own geometric graphs until no gain, then
//     a soft-shape descent with exact LP re-solves on that topology
//     (shape_lp, s_pol_shape) and a last re-polish.
//   * A polished placement replaces the output only if its WL is lower, it
//     passes the gate mirror and its boundary/MIB counts are >= those of the
//     placement it came from (shapes are unchanged by a polish; boundary
//     edges already touched stay pinned).
//   * s_lpcost (off by default; it did not pay on dev): over the last
//     s_lp_frac of the post-feasible phase, the WL term of a feasible
//     candidate is its LP WL on the SP graphs (warm-started AxisSolver per
//     axis); an accepted LP placement is offered as output.
//
// Phase E (strong only; knobs in Config, all tuned at --time 110):
//   * Post-feasible schedule: Tend = s_tend * T0 (3e-3; Phase D used tend =
//     1e-4, which spent the late budget frozen), T = T2 (Tend/T2)^(g^s_cool_p); s_p2_cyc equal-time
//     segments, segment k >= 1 restarting at s_p2_rh * T2 from the best packed
//     incumbent (s_p2_best); s_stall_k: after s_stall_k * n moves without a
//     new packed incumbent, re-anneal from the incumbent at
//     min(T2, s_stall_t * T) over the remaining time.
//   * Restart portfolio (s_rs, run_portfolio()), only when run 1 is feasible
//     at the split point (else run 1 continues alone to the deadline):
//     s_rs=1: run 1 anneals to s_rs_frac * deadline and is finalized, run 2
//     (another seed) runs to the deadline, the better output wins; s_rs=2:
//     runs 1 and 2 each run s_rs_frac * deadline on schedules planned to
//     (1 - s_rs_frac) * deadline, then the better one is resumed. A run
//     writes the output file only while it is the better one. RESULT ttf is
//     run 1's (the first feasible time of the process), moves the sum.
//   * Final polish rounds (s_sp_repol): while a round gains, re-polish the
//     output on the SP graphs of SPs extracted from it, then soft-shape
//     descent and geometric re-polishes again.
//   * s_wl0=1 (calibration sanity run, docs/design.md 9.5): WL weight 0 in the
//     SA cost, random start (s_seeds=0), no polish, and the run stops at the
//     first feasible state (the output is that state).
//
// Phase F (strong only; knobs in Config, tuned at --time 110 on dev_tight):
//   * GLS legalizer (solver/gls.h, coordinate-space overlap removal with
//     guided-local-search weights), Solver::gls_try(): legalize a placement,
//     re-extract the SP (keys X2 - Y2 / X2 + Y2 on doubled centroids), accept
//     it only if it keeps the invariant, adds no forbidden relation and packs
//     feasible (+ gate mirror); consider() then sets ttf / FIRSTFEAS.
//       - GLS-first (s_gls=1): right after the seeds, from the best seed's
//         packing, budget min(s_gls_t=30 s, s_gls_frac=0.25 x remaining);
//       - plateau GLS: on a pre-feasible plateau restart, from the attempt's
//         best state, budget min(10 s, 25% remaining), <= s_gls_plat=3 calls;
//       - after a hand-off the post-feasible anneal starts at s_gls_t2 x T0
//         (0.02; hotter starts, 0.2 / 0.5, were worse).
//   * Homotopy (s_shrink=0.92), run_homotopy(); phases are fractions of the
//     budget Rem left when it starts:
//       1. the strong flow (seeds, GLS-first, post-feasible SA) on the relaxed
//          instance relax_problem(s_shrink) (fixed dims floor(true sqrt s),
//          soft shape lists scaled by sqrt s with the same indices, the TRUE
//          forbidden table) until s_h1=0.35 x Rem; a muted Solver (no_gate,
//          never writes): each relaxed incumbent is packed at true dims and
//          offered to the true solver only if that packing is feasible;
//       2. s_h_mode=0: predictor-corrector ramp s -> 1 (ported from
//          test_homotopy.cpp) for s_h2=0.15 x Rem, then at a stall (s >=
//          0.95) or at the end GLS from the packing at s = 1; s_h_mode=1
//          (default): no ramp, GLS from the relaxed packed incumbent with true
//          shapes; GLS budget min(s_h_gls=10 s, 0.15 x Rem);
//       3. post-feasible anneal (start s_gls_t2 x T0) from the best packed
//          true incumbent, then the usual final polish. If nothing legal came
//          out, the normal true-instance flow runs on the remaining time.
//   * Hybrid (s_hyb=1, default): the step-1 flow (seeds, GLS-first, SA) runs
//     until s_hyb_at=0 s after its first feasible state; if >= s_hyb_min=0.5
//     x the deadline remains then, the rest of the budget goes to the
//     homotopy. The homotopy only offers candidates through consider()
//     (feasible + gate mirror + lower WL) and the switch happens only once
//     feasible, so feasibility equals the step-1 flow's (same prefix).
//     s_hyb=0 with s_shrink > 0 is the standalone homotopy (whole budget;
//     ttf 40-50 s on dev_tight; its Phase F "feasibility losses" were
//     crashes, an empty SP packed after a failed absorber, fixed in Phase G);
//     s_hyb=0 s_shrink=0 is the step-1 flow.
//
// Phase G (strong only; A/Bs at --time 110 on dev_tight (16 cases) + dev_cal (8
// cases, ws 3.0-6.5%) x 5 seeds, same-window mixed queues against the Phase F
// binary; accept rule: feasible >= control and 0 verifier disagreements on both
// sets, then pooled median WL/construction, then the tail (p90/median, max/median,
// share of runs above 1.05 x their case median; the score gives 0 WL points above
// 1.05 T). Noise floor between identical binaries at 5 seeds: 1-2.5% median.
//   * Fixes: the standalone homotopy (s_hyb=0, s_shrink>0) crashed (SIGSEGV) when
//     a GLS absorber failed and an empty SP was packed; S now always holds a valid
//     state and gls_try() restores K only for a non-empty S. Phase F's "feasibility
//     losses" of the standalone path were these crashes (corrected numbers in
//     docs/research/adhd-solver-round3.md).
//   * Output refresh: the true output file is kept fresh while the relaxed
//     homotopy run owns the clock (Solver::maybe_write through tru) and flushed
//     at every homotopy phase boundary, so a SIGKILL at any time leaves the best
//     feasible placement found so far (tested at 30 / 60 / 90 s).
//   * Defaults are unchanged. s_gls_frac 0.35 and s_rng_split 1 looked better in
//     the first 5-seed sweep (dev_tight 1.052 vs 1.063, dev_cal 40/40 vs 37/40) but
//     did not replicate: pooled over 40 seeds per case against Phase F they are
//     dev_tight 640/640 vs 638/640, median 1.056 vs 1.055, tail share 0.172 vs
//     0.196; dev_cal 300/320 vs 305/320, 1.125 vs 1.120, 0.267 vs 0.252 (at
//     --time 110 s_gls_frac only moves the GLS-first budget from 27 s to the 30 s
//     cap, and s_rng_split cannot affect feasibility). Phase G behaves like
//     Phase F by construction (same defaults, same RNG streams; only per-move
//     statistics bookkeeping was added, 369 k vs 383 k moves/s); its value is the
//     fixes, the measurements below and the knobs / instrumentation they left
//     behind. Battery of this binary vs Phase F (--time 110): dev_tight x 10
//     160/160 vs 159/160, median 1.054 vs 1.037; dev_cal x 10 79/80 vs 76/80,
//     1.094 vs 1.128; dev x 5 100/100 both; benchmarks x 5 at 60 s 25/25 both.
//   * Measured and rejected (knobs kept, default off; numbers are the change of
//     the median WL/construction against the bench's base arm, dev_tight /
//     dev_cal: Phase F for the acceptance-rule and schedule arms, C for the
//     racing and directed-move arms; C and F differ by <= 1.2% in those benches):
//       - acceptance rules s_acc=1 LAHC (s_lahc_len 2000 / 20000: +9% / +4% median
//         on dev_tight, +8% / +16% on dev_cal), s_acc=2 Modified Lam (a0 0.06 /
//         0.15: +7% / +51%); s_quench 0.15 / 0.2 (greedy tail: +2% tight, +3% /
//         +7% cal); s_tend 1e-4 / 1e-3 / 6e-3 / 1e-2 (+6 / +1.6 / +1 / +4% tight,
//         +8 / +2 / +3 / +7% cal). The post-feasible anneal runs from T2 = 0.02 T0
//         (after the GLS hand-off) to Tend = 0.003 T0, i.e. its uphill acceptance
//         only falls from 5.3% to 4.8% (STATS phase2 up_by_tenth); both a colder
//         and a hotter end are worse, so the schedule is not the lever.
//       - directed phase-2 moves s_p_med (net-median insertion) / s_p_mib (MIB
//         member swap) at 0.05 / 0.1: within noise on dev_tight (-1.5% .. +0%),
//         +1.8 .. +3.8% on dev_cal with a larger tail.
//       - phase-1 race s_race_k / s_race_t / s_race_keep (k relaxed runs on their
//         own RNG streams and clocks, paused after s_race_t s, the best resumed):
//         k3@8s / k3@16s / k5@8s / k5@8s keep 2 / k5@16s: no arm within +0.5%
//         median with a smaller tail on both sets (k5@8s has the smallest dev_cal
//         tail measured, p90/median 1.087 vs 1.150, at +2.4% median and a worse
//         dev_tight tail). k5@16s ran before the race budget was clamped to the
//         time left (run_homotopy) and overran the deadline on 12 of 120 runs
//         (wall up to 136.8 s after a 45 s first-feasible time); its numbers
//         (+7% / +11%) are tainted.
//         Diagnostic (s_trace=1, 6 tail-heavy cases x 10 seeds): Spearman rho of the relaxed incumbent WL vs the final WL is 0.32 / 0.39 /
//         0.48 / 0.61 at 10 / 20 / 30 s / end of phase 1, of the true incumbent
//         0.01 / 0.03 / 0.54 / 0.94 at 10 / 20 / 50 / 100 s: the step-1 state does
//         not steer the result (the relaxed run starts from its own seeds), and
//         half of the seed-to-seed rank is decided during the true anneal.
//       - s_h_gls_feas=3 (absorber cap when a true incumbent exists): worse.
//
// stderr ends with one line:
//   RESULT feasible=<0|1> ttf=<s|-1> wl=<x.x> boundary=<sat>/<p> mib=<sat>/<q> moves=<N>
#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "extract.h"
#include "lab3_io.h"
#include "mcf/wl_lp.h"
#include "shape.h"

namespace {

using Clock = std::chrono::steady_clock;

// ------------------------------------------------------------------ config

struct Config {
    double time = 24.0;
    uint64_t seed = 1;
    bool strong = true;          // false: frozen Phase A (weak) path
    std::string init_from;       // seed SP from a placement (E1 / Phase C hook)
    bool pack_only = false;      // with init_from: pack + flush, no SA
    // Tuning knobs (--set key=value).
    double p_shape = 0.15;       // share of shape moves
    double p_local = 0.50;       // SP partner chosen near in the sequence
    double p_target = 0.05;      // targeted repair move while infeasible
    double cycle_k = 20000;      // pre-feasible cycle length = cycle_k * n moves
    double reheat = 0.30;        // pre-feasible cycle restart temperature / T0
    double t2_floor = 0.05;      // post-feasible anneal starts at >= t2_floor * T0
    double t2_cap = 1.0;         // ... and at <= t2_cap * T0
    double tend = 1e-4;          // final temperature / T0
    double acc0 = 0.30;          // initial uphill acceptance probability
    double lam_max = 1e3;        // cap on the adaptive penalty weight
    // Strong-mode knobs (Phase C; ignored by --mode weak).
    double s_seeds = 8;          // analytical seeds generated at start (0: random start)
    double s_gp_iters = 25;      // spreading iterations per analytical seed
    double s_gp_anchor = 0.05;   // anchor weight growth per spreading iteration
    double s_grow = 1;           // seed SP from grown rectangles (1) or points (0)
    double s_geo_prio = 1;       // seed SP tie-breaks from geometry (1) or the slicing tree (0)
    double s_shape_desc = 1;     // greedy shape descent on each seed
    double s_first = 0;          // pool index of the first start (rank experiments)
    double s_tstart = 0.05;      // seeded start temperature / probed T
    double s_cycle_k = 20000;    // pre-feasible cycle length = s_cycle_k * n moves
    double s_p_target = 0.05;    // targeted repair share while infeasible
    double s_p_ins = 0;          // share of valid-window insert moves
    double s_forbid = 1;         // never add a geometrically impossible block-preplaced relation
    double s_plateau_k = 1500;   // plateau = s_plateau_k * n moves without a better penalty
    double s_plat_restart = 1;   // plateau: restart from the next start (1) or reheat (0)
    double s_restart_t = 0.3;    // restart start T / probed T
    double s_kick_cycle_k = 1500;  // reheat cycle length = s_kick_cycle_k * n moves
    double s_near = 0.005;       // plateau with best penalty <= s_near * (W + H): reheat, not restart
    double s_polish = 1;         // exact LP WL polish of incumbents (solver/mcf) for the output
    double s_pol_ms = 250;       // min ms between two polishes during the run (final polish always)
    double s_pol_budget = 0.02;  // ... and only while run-time polishing used <= this share of the elapsed time
    double s_pol_geo = 1;        // also polish on the geometric constraint graphs (1) or SP graphs only (0)
    double s_pol_shape = 1;      // final: soft-shape descent on the polished output's topology (exact LP)
    double s_lpcost = 0;         // late SA: cost of a feasible candidate = LP WL of its SP topology
    double s_lp_frac = 0.2;      // ... over the last s_lp_frac of the post-feasible phase
    // Phase E (strong): post-feasible schedule, restart portfolio, sanity run.
    double s_tend = 3e-3;        // post-feasible final temperature / T0 (strong; weak uses tend)
    double s_cool_p = 1;         // post-feasible T = T2 (Tend/T2)^(g^s_cool_p), g = segment time fraction
    double s_p2_cyc = 1;         // post-feasible anneal split into this many equal-time segments
    double s_p2_rh = 0.3;        // ... segment k >= 1 starts at s_p2_rh * T2
    double s_p2_best = 1;        // ... and (and a stall reheat) from the best packed incumbent (1) or the current state (0)
    double s_stall_k = 0;        // post-feasible stall: s_stall_k * n moves without a new packed incumbent (0: off) ...
    double s_stall_t = 3;        // ... re-anneals from min(T2, s_stall_t * T) to Tend over the remaining time
    double s_rs = 0;             // restart portfolio: 0 off, 1 two half-budget runs (best-of), 2 race then continue the better
    double s_rs_frac = 0.5;      // split point / budget (mode 1: end of run 1; mode 2: race length of each run)
    double s_wl0 = 0;            // sanity: WL weight 0 (random start, no polish, stop at the first feasible state)
    double s_sp_repol = 1;       // final: re-polish the output on the SP graphs of SPs extracted from it
    // Phase F (strong): GLS legalizer (gls.h) and the homotopy flow.
    double s_gls = 1;            // GLS-first from the best seed + GLS on pre-feasible plateaus (0: off)
    double s_gls_t = 30;         // GLS-first budget: min(s_gls_t, s_gls_frac * remaining time) seconds
    double s_gls_frac = 0.25;
    double s_gls_plat = 3;       // max GLS calls on pre-feasible plateaus (budget min(10 s, 25% of the remaining))
    double s_gls_t2 = 0.02;      // post-feasible start temperature / T0 right after a GLS / homotopy hand-off
    double s_shrink = 0.92;      // homotopy: relaxed dimension scale s0 (movable areas x s0; 0: off)
    double s_h1 = 0.35;          // ... phase 1 (strong flow on the relaxed instance) ends at s_h1 * deadline
    double s_h2 = 0.15;          // ... phase 2 (ramp s0 -> 1) gets s_h2 * deadline
    double s_h_mode = 1;         // ... 0: ramp, then GLS absorber at a stall; 1: no ramp, GLS from the relaxed placement
    double s_h_gls = 10;         // ... absorber GLS budget: min(s_h_gls s, 0.15 x the homotopy's budget)
    double s_hyb = 1;            // hybrid: step-1 flow until s_hyb_at s after the first feasible state, then the homotopy (s_shrink)
    double s_hyb_at = 0;         // ... seconds after the first feasible state
    double s_hyb_min = 0.5;      // ... only if >= s_hyb_min x the deadline remains then (else the step-1 flow continues)
    // Phase G (strong). s_acc, s_p_med, s_p_mib act in the true-instance phase-2 anneal only;
    // s_quench and s_trace also in the relaxed homotopy run; s_race_* in run_homotopy. All
    // measured worse than or equal to the defaults (header, Phase G) and kept off.
    double s_rng_split = 0;      // homotopy: relaxed run on its own RNG seed (seed ^ 0xA5A5A5A5)
    double s_h_gls_feas = -1;    // homotopy absorber budget cap when a true incumbent exists (-1: s_h_gls rule only)
    double s_acc = 0;            // 0 Metropolis on the time schedule, 1 Late Acceptance Hill Climbing, 2 Modified Lam
    double s_lahc_len = 2000;    // LAHC history length (moves)
    double s_lam_a0 = 1.0;       // Modified Lam target acceptance at the phase start ...
    double s_lam_mid = 0.44;     // ... plateau value (f in [0.15, 0.65]); decays to s_lam_mid / 440 at the end
    double s_quench = 0;         // last s_quench of phase 2 at T = 0 (improving and neutral moves only)
    double s_p_med = 0;          // directed: share of net-median insertions in phase 2 (block next to the block nearest its net median)
    double s_p_mib = 0;          // directed: share of MIB-member swaps (same shape, both sequences)
    double s_trace = 0;          // diagnostic: TRACE line (incumbent WL) every s_trace seconds (0: off)
    double s_race_k = 0;         // homotopy phase 1 race: this many relaxed runs (own RNG streams), 0/1: off
    double s_race_t = 8;         // ... each paused after s_race_t s of its own run
    double s_race_keep = 1;      // ... the best s_race_keep resume (one after the other) on the rest of phase 1
    double s_race_min = 0.1;     // ... phase 1 lasts max(s_h1, (k s_race_t) / Rem + s_race_min) x Rem
    bool set(const std::string &kv) {
        size_t e = kv.find('=');
        if (e == std::string::npos) return false;
        std::string k = kv.substr(0, e);
        double v = std::atof(kv.c_str() + e + 1);
        if (k == "p_shape") p_shape = v;
        else if (k == "p_local") p_local = v;
        else if (k == "p_target") p_target = v;
        else if (k == "cycle_k") cycle_k = v;
        else if (k == "reheat") reheat = v;
        else if (k == "t2_floor") t2_floor = v;
        else if (k == "t2_cap") t2_cap = v;
        else if (k == "tend") tend = v;
        else if (k == "acc0") acc0 = v;
        else if (k == "lam_max") lam_max = v;
        else if (k == "s_seeds") s_seeds = v;
        else if (k == "s_gp_iters") s_gp_iters = v;
        else if (k == "s_gp_anchor") s_gp_anchor = v;
        else if (k == "s_grow") s_grow = v;
        else if (k == "s_geo_prio") s_geo_prio = v;
        else if (k == "s_shape_desc") s_shape_desc = v;
        else if (k == "s_first") s_first = v;
        else if (k == "s_tstart") s_tstart = v;
        else if (k == "s_cycle_k") s_cycle_k = v;
        else if (k == "s_p_target") s_p_target = v;
        else if (k == "s_p_ins") s_p_ins = v;
        else if (k == "s_forbid") s_forbid = v;
        else if (k == "s_plateau_k") s_plateau_k = v;
        else if (k == "s_plat_restart") s_plat_restart = v;
        else if (k == "s_restart_t") s_restart_t = v;
        else if (k == "s_kick_cycle_k") s_kick_cycle_k = v;
        else if (k == "s_near") s_near = v;
        else if (k == "s_polish") s_polish = v;
        else if (k == "s_pol_ms") s_pol_ms = v;
        else if (k == "s_pol_budget") s_pol_budget = v;
        else if (k == "s_pol_geo") s_pol_geo = v;
        else if (k == "s_pol_shape") s_pol_shape = v;
        else if (k == "s_lpcost") s_lpcost = v;
        else if (k == "s_lp_frac") s_lp_frac = v;
        else if (k == "s_tend") s_tend = v;
        else if (k == "s_cool_p") s_cool_p = v;
        else if (k == "s_p2_cyc") s_p2_cyc = v;
        else if (k == "s_p2_rh") s_p2_rh = v;
        else if (k == "s_p2_best") s_p2_best = v;
        else if (k == "s_stall_k") s_stall_k = v;
        else if (k == "s_stall_t") s_stall_t = v;
        else if (k == "s_rs") s_rs = v;
        else if (k == "s_rs_frac") s_rs_frac = v;
        else if (k == "s_wl0") s_wl0 = v;
        else if (k == "s_sp_repol") s_sp_repol = v;
        else if (k == "s_gls") s_gls = v;
        else if (k == "s_gls_t") s_gls_t = v;
        else if (k == "s_gls_frac") s_gls_frac = v;
        else if (k == "s_gls_plat") s_gls_plat = v;
        else if (k == "s_gls_t2") s_gls_t2 = v;
        else if (k == "s_shrink") s_shrink = v;
        else if (k == "s_h1") s_h1 = v;
        else if (k == "s_h2") s_h2 = v;
        else if (k == "s_h_mode") s_h_mode = v;
        else if (k == "s_h_gls") s_h_gls = v;
        else if (k == "s_hyb") s_hyb = v;
        else if (k == "s_hyb_at") s_hyb_at = v;
        else if (k == "s_hyb_min") s_hyb_min = v;
        else if (k == "s_acc") s_acc = v;
        else if (k == "s_quench") s_quench = v;
        else if (k == "s_p_med") s_p_med = v;
        else if (k == "s_p_mib") s_p_mib = v;
        else if (k == "s_trace") s_trace = v;
        else if (k == "s_race_k") s_race_k = v;
        else if (k == "s_race_t") s_race_t = v;
        else if (k == "s_race_keep") s_race_keep = v;
        else if (k == "s_race_min") s_race_min = v;
        else if (k == "s_rng_split") s_rng_split = v;
        else if (k == "s_h_gls_feas") s_h_gls_feas = v;
        else if (k == "s_lahc_len") s_lahc_len = v;
        else if (k == "s_lam_a0") s_lam_a0 = v;
        else if (k == "s_lam_mid") s_lam_mid = v;
        else return false;
        return true;
    }
};

// ------------------------------------------------------------------ rng

struct Rng {
    uint64_t s[4];
    static uint64_t splitmix(uint64_t &x) {
        uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    explicit Rng(uint64_t seed) {
        uint64_t x = seed;
        for (auto &v : s) v = splitmix(x);
    }
    static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    uint64_t next() {
        uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t; s[3] = rotl(s[3], 45);
        return r;
    }
    int below(int n) { return (int)(((next() >> 32) * (uint64_t)n) >> 32); }
    double uni() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
};

// ------------------------------------------------------------------ problem

enum { K_SOFT = 0, K_FIXED = 1, K_PRE = 2 };

struct BBConn { int a, b, w; };
struct BTConn { int a, tx2, ty2, w; };

struct Problem {
    lab3::Instance inst;
    int n = 0;
    int W = 0, H = 0;
    std::vector<int> kind, fw, fh, px, py;   // fixed dims / preplaced position
    std::vector<int> group;                  // shape group per block, -1 if not soft
    std::vector<std::vector<int>> gmem;      // members per shape group
    std::vector<fp::ShapeRange> gshape;      // legal widths per shape group
    std::vector<int> bmask;                  // boundary edges per block
    std::vector<int> pre, right, top;        // preplaced / RIGHT / TOP block lists
    std::vector<BBConn> bb;
    std::vector<BTConn> bt;
    // Strong mode only (build_forbidden): block-preplaced relations that no
    // legal placement can have, e.g. "a left-of p" when p.x < a's narrowest
    // width. fmask[a * |pre| + preidx[p]] bit r forbids rel(a, p) == r.
    std::vector<int> preidx;
    std::vector<uint8_t> fmask;
    std::vector<std::array<int, 3>> forb;  // (a, p, mask), mask != 0
};

// SP relation of a to b: 0 a left-of b, 1 b left-of a, 2 a below b, 3 b below a.
inline int sp_rel(const std::vector<int> &pp, const std::vector<int> &pm, int a, int b) {
    bool bp = pp[a] < pp[b], bm = pm[a] < pm[b];
    return bp ? (bm ? 0 : 3) : (bm ? 2 : 1);
}

void build_forbidden(Problem &P) {
    const int n = P.n, np = (int)P.pre.size();
    P.preidx.assign(n, -1);
    for (int i = 0; i < np; ++i) P.preidx[P.pre[i]] = i;
    P.fmask.assign((size_t)n * std::max(1, np), 0);
    for (int a = 0; a < n; ++a) {
        if (P.kind[a] == K_PRE) continue;
        int wmin, hmin;
        if (P.kind[a] == K_SOFT) {
            const fp::ShapeRange &g = P.gshape[P.group[a]];
            wmin = g.widths.front();
            hmin = *std::min_element(g.heights.begin(), g.heights.end());
        } else {
            wmin = P.fw[a];
            hmin = P.fh[a];
        }
        for (int i = 0; i < np; ++i) {
            int p = P.pre[i], m = 0;
            if (wmin > P.px[p]) m |= 1;
            if (wmin > P.W - P.px[p] - P.fw[p]) m |= 2;
            if (hmin > P.py[p]) m |= 4;
            if (hmin > P.H - P.py[p] - P.fh[p]) m |= 8;
            P.fmask[(size_t)a * np + i] = (uint8_t)m;
            if (m) P.forb.push_back({a, p, m});
        }
    }
}

int64_t forb_count(const Problem &P, const struct State &S);

bool build_problem(const std::string &path, Problem &P, std::string &err) {
    if (!lab3::parse_input_file(path, P.inst, err)) return false;
    const lab3::Instance &I = P.inst;
    P.n = (int)I.blocks.size();
    if (I.W > (1 << 28) || I.H > (1 << 28)) { err = "outline too large for 32-bit coordinates"; return false; }
    P.W = (int)I.W;
    P.H = (int)I.H;
    const int n = P.n;
    P.kind.assign(n, 0); P.fw.assign(n, 0); P.fh.assign(n, 0); P.px.assign(n, 0); P.py.assign(n, 0);
    P.group.assign(n, -1); P.bmask.assign(n, 0);
    std::vector<int> mib_of(n, -1);
    for (size_t g = 0; g < I.mibs.size(); ++g) for (int b : I.mibs[g]) mib_of[b] = (int)g;
    std::vector<int> mib_group(I.mibs.size(), -1);
    for (int i = 0; i < n; ++i) {
        const lab3::Block &b = I.blocks[i];
        if (b.kind == lab3::Kind::Soft) {
            P.kind[i] = K_SOFT;
            int g;
            if (mib_of[i] >= 0 && mib_group[mib_of[i]] >= 0) {
                g = mib_group[mib_of[i]];
            } else {
                g = (int)P.gmem.size();
                P.gmem.emplace_back();
                P.gshape.push_back(fp::make_shape_range(b.area, I.R));
                if (P.gshape.back().widths.empty()) { err = "soft block " + b.name + " has no legal shape"; return false; }
                if (mib_of[i] >= 0) mib_group[mib_of[i]] = g;
            }
            P.group[i] = g;
            P.gmem[g].push_back(i);
        } else {
            P.kind[i] = (b.kind == lab3::Kind::Fixed) ? K_FIXED : K_PRE;
            P.fw[i] = (int)b.w;
            P.fh[i] = (int)b.h;
            if (P.kind[i] == K_PRE) {
                P.px[i] = (int)b.x;
                P.py[i] = (int)b.y;
                P.pre.push_back(i);
            }
        }
    }
    for (const lab3::Boundary &bd : I.boundaries) P.bmask[bd.block] = bd.edges;
    for (int i = 0; i < n; ++i) {
        if (P.bmask[i] & lab3::EDGE_RIGHT) P.right.push_back(i);
        if (P.bmask[i] & lab3::EDGE_TOP) P.top.push_back(i);
    }
    for (const lab3::Conn &c : I.conns) {
        if (c.a >= 0 && c.b >= 0) {
            P.bb.push_back({c.a, c.b, (int)c.weight});
        } else {
            int blk = c.a >= 0 ? c.a : c.b, t = c.a >= 0 ? c.b : c.a;
            const lab3::Terminal &tm = I.terminals[-t - 1];
            P.bt.push_back({blk, (int)(2 * tm.x), (int)(2 * tm.y), (int)c.weight});
        }
    }
    return true;
}

// ------------------------------------------------------------------ state

struct State {
    std::vector<int> gp, gm;   // Γ+, Γ- (block ids)
    std::vector<int> pp, pm;   // positions of each block in Γ+, Γ-
    std::vector<int> gidx;     // per shape group: index into its widths list
    std::vector<int> w, h;     // current block dimensions

    void swap_in(State &o) { gp.swap(o.gp); gm.swap(o.gm); pp.swap(o.pp); pm.swap(o.pm); gidx.swap(o.gidx); w.swap(o.w); h.swap(o.h); }
    void sync_pos() {
        for (int i = 0; i < (int)gp.size(); ++i) pp[gp[i]] = i;
        for (int i = 0; i < (int)gm.size(); ++i) pm[gm[i]] = i;
    }
};

int64_t forb_count(const Problem &P, const State &S) {
    int64_t c = 0;
    for (const auto &f : P.forb) c += (f[2] >> sp_rel(S.pp, S.pm, f[0], f[1])) & 1;
    return c;
}

// Transitively reduced constraint graphs of an SP (strong mode, LP polish):
// h holds (i, j) for i left-of j, v holds (i, j) for i below j, only for
// pairs with no block between them on that axis. For i, the blocks after i
// in Γ+ are scanned in order; a right neighbour j (after i in Γ-) is direct
// iff no earlier-scanned right neighbour precedes it in Γ- (and mirrored for
// the blocks below i). O(n^2), no bitsets.
void sp_graphs(const std::vector<int> &gp, const std::vector<int> &pm, mcf::ConstraintGraphs &g) {
    const int n = (int)gp.size();
    g.h.clear();
    g.v.clear();
    thread_local std::vector<int> q;   // Γ- position by Γ+ position
    q.resize(n);
    for (int c = 0; c < n; ++c) q[c] = pm[gp[c]];
    for (int a = 0; a < n; ++a) {
        const int i = gp[a], mi = q[a];
        int minr = n, maxb = -1;
        for (int c = a + 1; c < n && (minr > mi + 1 || maxb < mi - 1); ++c) {
            const int mj = q[c];
            if (mj > mi ? mj < minr : mj > maxb) {
                if (mj > mi) { g.h.emplace_back(i, gp[c]); minr = mj; }
                else { g.v.emplace_back(gp[c], i); maxb = mj; }
            }
        }
    }
}

void apply_group(const Problem &P, State &S, int g) {
    int w = P.gshape[g].widths[S.gidx[g]], h = P.gshape[g].heights[S.gidx[g]];
    for (int b : P.gmem[g]) { S.w[b] = w; S.h[b] = h; }
}

void init_dims(const Problem &P, State &S) {
    S.w.assign(P.n, 0);
    S.h.assign(P.n, 0);
    for (int i = 0; i < P.n; ++i) if (P.kind[i] != K_SOFT) { S.w[i] = P.fw[i]; S.h[i] = P.fh[i]; }
    for (int g = 0; g < (int)P.gmem.size(); ++g) apply_group(P, S, g);
}

// Boundary source/sink invariant + preplaced-pair geometric consistency.
int64_t g_rej_boundary = 0, g_rej_pre = 0;  // rejection causes (stats only)
bool invariant_ok(const Problem &P, const State &S) {
    const int n = P.n;
    int mn = INT_MAX, mx = -1;
    for (int i = 0; i < n; ++i) {
        int b = S.gp[i], m = S.pm[b], k = P.bmask[b];
        if (k) {
            if ((k & lab3::EDGE_LEFT) && mn < m) { ++g_rej_boundary; return false; }  // someone left-of b
            if ((k & lab3::EDGE_TOP) && mx > m) { ++g_rej_boundary; return false; }   // someone above b
        }
        if (m < mn) mn = m;
        if (m > mx) mx = m;
    }
    mn = INT_MAX; mx = -1;
    for (int i = n - 1; i >= 0; --i) {
        int b = S.gp[i], m = S.pm[b], k = P.bmask[b];
        if (k) {
            if ((k & lab3::EDGE_RIGHT) && mx > m) { ++g_rej_boundary; return false; }   // someone right-of b
            if ((k & lab3::EDGE_BOTTOM) && mn < m) { ++g_rej_boundary; return false; }  // someone below b
        }
        if (m < mn) mn = m;
        if (m > mx) mx = m;
    }
    const int np = (int)P.pre.size();
    for (int i = 0; i < np; ++i) {
        for (int j = i + 1; j < np; ++j) {
            int a = P.pre[i], b = P.pre[j];
            if (S.pp[a] > S.pp[b]) std::swap(a, b);           // a before b in Γ+
            if (S.pm[a] < S.pm[b]) {                          // a left-of b
                if (P.px[a] + P.fw[a] > P.px[b]) { ++g_rej_pre; return false; }
            } else {                                          // b below a
                if (P.py[b] + P.fh[b] > P.py[a]) { ++g_rej_pre; return false; }
            }
        }
    }
    return true;
}

// ------------------------------------------------------------------ packer

struct Packer {
    int n = 0;
    std::vector<int> bit;
    std::vector<int> x, y;
    int64_t ovfx = 0, ovfy = 0, viol = 0, wl2 = 0;

    void init(int nn) { n = nn; bit.assign(n + 1, 0); x.assign(n, 0); y.assign(n, 0); }
    int query(int k) const {  // max over Fenwick indices 1..k
        int r = 0;
        for (; k > 0; k -= k & -k) if (bit[k] > r) r = bit[k];
        return r;
    }
    void update(int k, int v) {
        for (; k <= n; k += k & -k) if (bit[k] < v) bit[k] = v;
    }

    void pack(const Problem &P, const State &S) {
        viol = 0;
        std::fill(bit.begin(), bit.end(), 0);
        for (int i = 0; i < n; ++i) {  // x: predecessors = before in Γ+ and Γ-
            int b = S.gp[i], k = S.pm[b], q = query(k), xb;
            if (P.kind[b] == K_PRE) {
                xb = P.px[b];
                if (q > xb) viol += q - xb;
            } else {
                xb = q;
            }
            x[b] = xb;
            update(k + 1, xb + S.w[b]);
        }
        std::fill(bit.begin(), bit.end(), 0);
        for (int i = n - 1; i >= 0; --i) {  // y: below = after in Γ+, before in Γ-
            int b = S.gp[i], k = S.pm[b], q = query(k), yb;
            if (P.kind[b] == K_PRE) {
                yb = P.py[b];
                if (q > yb) viol += q - yb;
            } else {
                yb = q;
            }
            y[b] = yb;
            update(k + 1, yb + S.h[b]);
        }
        // Flush sinks (no successor on that axis, so moving them out is safe).
        for (int b : P.right) if (x[b] + S.w[b] <= P.W) x[b] = P.W - S.w[b];
        for (int b : P.top) if (y[b] + S.h[b] <= P.H) y[b] = P.H - S.h[b];
        ovfx = ovfy = 0;
        for (int b = 0; b < n; ++b) {
            int ex = x[b] + S.w[b] - P.W, ey = y[b] + S.h[b] - P.H;
            if (ex > 0) ovfx += ex;
            if (ey > 0) ovfy += ey;
        }
        // 2*WL with doubled centroids, equal to lab3::wirelength2.
        int64_t s = 0;
        for (const BBConn &c : P.bb) {
            int dx = (2 * x[c.a] + S.w[c.a]) - (2 * x[c.b] + S.w[c.b]);
            int dy = (2 * y[c.a] + S.h[c.a]) - (2 * y[c.b] + S.h[c.b]);
            s += (int64_t)c.w * ((dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy));
        }
        for (const BTConn &c : P.bt) {
            int dx = (2 * x[c.a] + S.w[c.a]) - c.tx2;
            int dy = (2 * y[c.a] + S.h[c.a]) - c.ty2;
            s += (int64_t)c.w * ((dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy));
        }
        wl2 = s;
    }
    bool feasible() const { return ovfx == 0 && ovfy == 0 && viol == 0; }
};

// ------------------------------------------------------------------ gate mirror

struct GateInfo {
    bool ok = false;
    std::string why;
    int bsat = 0, msat = 0;
};

// Mirror of verifier gate checks 3-8 (design.md §6) plus boundary/MIB counts.
GateInfo gate_check(const Problem &P, const std::vector<lab3::Rect> &r) {
    const lab3::Instance &I = P.inst;
    GateInfo g;
    const int n = P.n;
    for (int i = 0; i < n; ++i) {
        const lab3::Block &b = I.blocks[i];
        const lab3::Rect &q = r[i];
        if (q.w <= 0 || q.h <= 0) { g.why = "size " + b.name; return g; }
        if (b.kind == lab3::Kind::Fixed && (q.w != b.w || q.h != b.h)) { g.why = "fixed " + b.name; return g; }
        if (b.kind == lab3::Kind::Preplaced && (q.x != b.x || q.y != b.y || q.w != b.w || q.h != b.h)) {
            g.why = "preplaced " + b.name; return g;
        }
        if (b.kind == lab3::Kind::Soft) {
            __int128 d = (__int128)q.w * q.h - b.area;
            if (d < 0) d = -d;
            if (100 * d > (__int128)b.area) { g.why = "area " + b.name; return g; }
            if ((__int128)I.R * q.h < q.w || (__int128)I.R * q.w < q.h) { g.why = "aspect " + b.name; return g; }
        }
        if (q.x < 0 || q.y < 0 || q.x + q.w > I.W || q.y + q.h > I.H) { g.why = "outline " + b.name; return g; }
    }
    for (int i = 0; i < n; ++i) {
        const lab3::Rect &a = r[i];
        for (int j = i + 1; j < n; ++j) {
            const lab3::Rect &b = r[j];
            if (std::max(a.x, b.x) < std::min(a.x + a.w, b.x + b.w) &&
                std::max(a.y, b.y) < std::min(a.y + a.h, b.y + b.h)) {
                g.why = "overlap " + I.blocks[i].name + " " + I.blocks[j].name;
                return g;
            }
        }
    }
    g.ok = true;
    for (const lab3::Boundary &bd : I.boundaries) {
        const lab3::Rect &q = r[bd.block];
        bool ok = true;
        if ((bd.edges & lab3::EDGE_LEFT) && q.x != 0) ok = false;
        if ((bd.edges & lab3::EDGE_RIGHT) && q.x + q.w != I.W) ok = false;
        if ((bd.edges & lab3::EDGE_BOTTOM) && q.y != 0) ok = false;
        if ((bd.edges & lab3::EDGE_TOP) && q.y + q.h != I.H) ok = false;
        g.bsat += ok;
    }
    for (const auto &grp : I.mibs) {
        bool ok = true;
        for (int b : grp) if (r[b].w != r[grp[0]].w || r[b].h != r[grp[0]].h) ok = false;
        g.msat += ok;
    }
    return g;
}

std::vector<lab3::Rect> to_rects(const Packer &K, const State &S) {
    std::vector<lab3::Rect> r(K.n);
    for (int i = 0; i < K.n; ++i) r[i] = lab3::Rect{K.x[i], K.y[i], S.w[i], S.h[i]};
    return r;
}

bool write_output(const Problem &P, const std::vector<lab3::Rect> &r, const std::string &path) {
    std::string tmp = path + ".tmp";
    {
        std::ofstream f(tmp);
        if (!f) return false;
        f << "WIRELENGTH " << lab3::format_half(lab3::wirelength2(P.inst, r)) << "\n";
        f << "BLOCKS " << P.n << "\n";
        for (int i = 0; i < P.n; ++i)
            f << P.inst.blocks[i].name << ' ' << r[i].x << ' ' << r[i].y << ' ' << r[i].w << ' ' << r[i].h << "\n";
        f.close();
        if (!f) { std::remove(tmp.c_str()); return false; }
    }
    return std::rename(tmp.c_str(), path.c_str()) == 0;
}

// ------------------------------------------------------------------ initial SP

// SP from a placement (all blocks). Soft (w,h) from the placement set the
// group width index (nearest legal width; MIB groups take the first member).
bool sp_from_rects(const Problem &P, State &S, const std::vector<fp::XRect> &r, Rng &rng, std::string &err) {
    const int n = P.n;
    std::vector<int64_t> kp(n), km(n);
    for (int i = 0; i < n; ++i) { kp[i] = (int64_t)rng.below(1 << 30); km[i] = (int64_t)rng.below(1 << 30); }
    if (!fp::extract_sp(r, kp, km, S.gp, S.gm, err)) return false;
    S.pp.assign(n, 0);
    S.pm.assign(n, 0);
    S.sync_pos();
    return true;
}

// Random initial SP that satisfies the invariant: extract it from a virtual
// non-overlapping placement on a doubled grid. Preplaced blocks keep their
// (doubled) rectangles; movable non-boundary blocks are 1x1 cells at random
// odd coordinates inside the outline and outside every preplaced block;
// boundary blocks sit in 1-wide lanes just outside their outline edge(s).
bool random_initial_sp(const Problem &P, State &S, Rng &rng, std::string &err) {
    const int n = P.n;
    const int64_t W2 = 2LL * P.W, H2 = 2LL * P.H;
    std::vector<fp::XRect> r(n);
    std::vector<std::pair<int64_t, int64_t>> used;
    auto inside_pre = [&](int64_t cx, int64_t cy) {
        for (int p : P.pre) {
            int64_t x0 = 2LL * P.px[p], y0 = 2LL * P.py[p];
            if (cx >= x0 && cx < x0 + 2LL * P.fw[p] && cy >= y0 && cy < y0 + 2LL * P.fh[p]) return true;
        }
        return false;
    };
    std::vector<char> usedx, usedy;  // keep cells distinct per axis (avoids ties)
    usedx.assign(P.W, 0);
    usedy.assign(P.H, 0);
    auto pick = [&](std::vector<char> &u, int lim) {
        for (int t = 0; t < 64; ++t) {
            int v = rng.below(lim);
            if (!u[v]) { u[v] = 1; return v; }
        }
        return rng.below(lim);
    };
    for (int i = 0; i < n; ++i) {
        if (P.kind[i] == K_PRE) {
            r[i] = {2LL * P.px[i], 2LL * P.py[i], 2LL * P.fw[i], 2LL * P.fh[i]};
            continue;
        }
        int k = P.bmask[i];
        int64_t cx, cy;
        if (k == 0) {
            int tries = 0;
            do {
                cx = 2LL * rng.below(P.W) + 1;
                cy = 2LL * rng.below(P.H) + 1;
            } while (inside_pre(cx, cy) && ++tries < 10000);
        } else {
            cx = 2LL * pick(usedx, P.W) + 1;
            cy = 2LL * pick(usedy, P.H) + 1;
            if (k & lab3::EDGE_LEFT) cx = -4;
            if (k & lab3::EDGE_RIGHT) cx = W2 + 3;
            if (k & lab3::EDGE_BOTTOM) cy = -4;
            if (k & lab3::EDGE_TOP) cy = H2 + 3;
        }
        r[i] = {cx, cy, 1, 1};
    }
    (void)used;
    return sp_from_rects(P, S, r, rng, err);
}

bool read_placement(const Problem &P, const std::string &path, std::vector<lab3::Rect> &r, std::string &err) {
    std::ifstream f(path);
    if (!f) { err = "cannot open " + path; return false; }
    std::string tok, wl;
    int64_t cnt;
    if (!(f >> tok >> wl) || tok != "WIRELENGTH") { err = "bad WIRELENGTH line"; return false; }
    if (!(f >> tok >> cnt) || tok != "BLOCKS" || cnt != P.n) { err = "bad BLOCKS line"; return false; }
    r.assign(P.n, lab3::Rect());
    std::vector<char> seen(P.n, 0);
    for (int i = 0; i < P.n; ++i) {
        std::string name;
        lab3::Rect q;
        if (!(f >> name >> q.x >> q.y >> q.w >> q.h)) { err = "bad block line"; return false; }
        auto it = P.inst.block_index.find(name);
        if (it == P.inst.block_index.end() || seen[it->second]) { err = "bad block name " + name; return false; }
        seen[it->second] = 1;
        r[it->second] = q;
    }
    return true;
}

bool initial_from_placement(const Problem &P, State &S, const std::string &path, Rng &rng, std::string &err) {
    std::vector<lab3::Rect> r;
    if (!read_placement(P, path, r, err)) return false;
    std::vector<fp::XRect> xr(P.n);
    for (int i = 0; i < P.n; ++i) xr[i] = {r[i].x, r[i].y, r[i].w, r[i].h};
    // Width index per shape group: the given width if legal, else the nearest.
    S.gidx.assign(P.gmem.size(), 0);
    for (int g = 0; g < (int)P.gmem.size(); ++g) {
        const auto &ws = P.gshape[g].widths;
        int w = (int)r[P.gmem[g][0]].w;
        int best = 0;
        for (int i = 0; i < (int)ws.size(); ++i)
            if (std::abs(ws[i] - w) < std::abs(ws[best] - w)) best = i;
        S.gidx[g] = best;
    }
    init_dims(P, S);
    return sp_from_rects(P, S, xr, rng, err);
}

// ------------------------------------------------------------------ analytical seed (strong)
//
// Quadratic placement of block centres (terminals and preplaced centres fixed;
// a boundary block's coordinate across its edge fixed at the edge), L1 via
// IRLS net weights 1/max(|d|, delta), spread SimPL-style: each iteration
// rough-legalizes by recursive area bisection of the outline (cut placed so
// the free area, preplaced area excluded, splits in proportion to the block
// area on each side) and pulls blocks toward their bisection targets with a
// growing anchor weight. The last bisection is a slicing floorplan of the
// movable blocks; its tree gives tie-break orders for extract_sp, its leaf
// regions give soft shapes, and its leaf centres are the points of the
// virtual placement (preplaced blocks as real rectangles, boundary blocks in
// lanes outside their edge), so the SP satisfies the invariant by construction.
struct GlobalPlacer {
    const Problem &P;
    Rng &rng;
    int n;
    std::vector<double> area, x, y, tx, ty, leafw, leafh;
    std::vector<char> fixx, fixy;
    std::vector<int> seqp, seqm;   // slicing-tree orders (movable blocks only)
    std::vector<double> sumc;      // per block: sum of connection weights
    double delta;
    bool rand_axis = false;

    GlobalPlacer(const Problem &p, Rng &r) : P(p), rng(r), n(p.n) {
        area.assign(n, 0); x.assign(n, 0); y.assign(n, 0); tx.assign(n, 0); ty.assign(n, 0);
        leafw.assign(n, 0); leafh.assign(n, 0); fixx.assign(n, 0); fixy.assign(n, 0); sumc.assign(n, 0);
        delta = std::max(1.0, 0.005 * (P.W + P.H));
        for (int i = 0; i < n; ++i) {
            double w, h;
            if (P.kind[i] == K_SOFT) { area[i] = (double)P.gshape[P.group[i]].A; w = h = std::sqrt(area[i]); }
            else { w = P.fw[i]; h = P.fh[i]; area[i] = w * h; }
            if (P.kind[i] == K_PRE) {
                x[i] = tx[i] = P.px[i] + w / 2; y[i] = ty[i] = P.py[i] + h / 2; fixx[i] = fixy[i] = 1;
                continue;
            }
            x[i] = P.W / 2.0 + (rng.uni() - 0.5) * 0.2 * P.W;
            y[i] = P.H / 2.0 + (rng.uni() - 0.5) * 0.2 * P.H;
            int k = P.bmask[i];
            if (k & lab3::EDGE_LEFT) { x[i] = w / 2; fixx[i] = 1; }
            if (k & lab3::EDGE_RIGHT) { x[i] = P.W - w / 2; fixx[i] = 1; }
            if (k & lab3::EDGE_BOTTOM) { y[i] = h / 2; fixy[i] = 1; }
            if (k & lab3::EDGE_TOP) { y[i] = P.H - h / 2; fixy[i] = 1; }
        }
        for (const BBConn &c : P.bb) { sumc[c.a] += c.w; sumc[c.b] += c.w; }
        for (const BTConn &c : P.bt) sumc[c.a] += c.w;
    }

    // One weighted quadratic solve on one axis (CG, Jacobi preconditioner).
    void solve_axis(bool ax, const std::vector<double> *tgt, double anc) {
        std::vector<double> &p = ax ? y : x;
        const std::vector<char> &fx = ax ? fixy : fixx;
        std::vector<double> d(n, 0), b(n, 0);
        std::vector<std::vector<std::pair<int, double>>> off(n);
        for (const BBConn &c : P.bb) {
            double wgt = c.w / std::max(std::fabs(p[c.a] - p[c.b]), delta);
            bool ma = !fx[c.a], mb = !fx[c.b];
            if (ma) { d[c.a] += wgt; if (mb) off[c.a].push_back({c.b, wgt}); else b[c.a] += wgt * p[c.b]; }
            if (mb) { d[c.b] += wgt; if (ma) off[c.b].push_back({c.a, wgt}); else b[c.b] += wgt * p[c.a]; }
        }
        for (const BTConn &c : P.bt) {
            if (fx[c.a]) continue;
            double t = (ax ? c.ty2 : c.tx2) / 2.0;
            double wgt = c.w / std::max(std::fabs(p[c.a] - t), delta);
            d[c.a] += wgt; b[c.a] += wgt * t;
        }
        double mid = (ax ? P.H : P.W) / 2.0;
        for (int i = 0; i < n; ++i) {
            if (fx[i]) continue;
            double eps = 1e-3 * (sumc[i] + 1) / (ax ? P.H : P.W);
            d[i] += eps; b[i] += eps * mid;
            if (tgt) {
                double aw = anc * (sumc[i] + 1) / std::max(std::fabs(p[i] - (*tgt)[i]), delta);
                d[i] += aw; b[i] += aw * (*tgt)[i];
            }
        }
        // CG on movable entries.
        std::vector<double> r(n, 0), z(n, 0), q(n, 0), Ap(n, 0);
        auto mul = [&](const std::vector<double> &v, std::vector<double> &o) {
            for (int i = 0; i < n; ++i) {
                if (fx[i]) { o[i] = 0; continue; }
                double s = d[i] * v[i];
                for (auto &e : off[i]) s -= e.second * v[e.first];
                o[i] = s;
            }
        };
        std::vector<double> v(n);
        for (int i = 0; i < n; ++i) v[i] = fx[i] ? 0 : p[i];
        mul(v, Ap);
        double rz = 0, bn = 0;
        for (int i = 0; i < n; ++i) {
            if (fx[i]) continue;
            r[i] = b[i] - Ap[i]; z[i] = r[i] / d[i]; q[i] = z[i]; rz += r[i] * z[i]; bn += b[i] * b[i];
        }
        for (int it = 0; it < 200 && rz > 1e-12 * bn; ++it) {
            mul(q, Ap);
            double qAq = 0;
            for (int i = 0; i < n; ++i) if (!fx[i]) qAq += q[i] * Ap[i];
            if (qAq <= 0) break;
            double al = rz / qAq, rz2 = 0;
            for (int i = 0; i < n; ++i) {
                if (fx[i]) continue;
                v[i] += al * q[i]; r[i] -= al * Ap[i]; z[i] = r[i] / d[i]; rz2 += r[i] * z[i];
            }
            double be = rz2 / rz;
            rz = rz2;
            for (int i = 0; i < n; ++i) if (!fx[i]) q[i] = z[i] + be * q[i];
        }
        double lim = ax ? P.H : P.W;
        for (int i = 0; i < n; ++i) if (!fx[i]) p[i] = std::min(lim, std::max(0.0, v[i]));
    }

    double free_area(double x0, double x1, double y0, double y1) const {
        double a = (x1 - x0) * (y1 - y0);
        for (int p : P.pre) {
            double ox = std::min(x1, (double)P.px[p] + P.fw[p]) - std::max(x0, (double)P.px[p]);
            double oy = std::min(y1, (double)P.py[p] + P.fh[p]) - std::max(y0, (double)P.py[p]);
            if (ox > 0 && oy > 0) a -= ox * oy;
        }
        return std::max(a, 0.0);
    }

    // Recursive area bisection; fills tx/ty/leafw/leafh and appends the
    // slicing-tree orders of the blocks in ids to sp/sm.
    void bisect(std::vector<int> &ids, double x0, double x1, double y0, double y1, std::vector<int> &sp,
                std::vector<int> &sm) {
        if (ids.empty()) return;
        if (ids.size() == 1) {
            int b = ids[0];
            tx[b] = (x0 + x1) / 2; ty[b] = (y0 + y1) / 2; leafw[b] = x1 - x0; leafh[b] = y1 - y0;
            sp.push_back(b); sm.push_back(b);
            return;
        }
        double wr = x1 - x0, hr = y1 - y0;
        bool vert = wr >= hr;  // vertical cut line (split along x)
        if (rand_axis && wr < 1.5 * hr && hr < 1.5 * wr) vert = rng.uni() < 0.5;
        const std::vector<double> &c = vert ? x : y;
        std::sort(ids.begin(), ids.end(), [&](int a, int b) { return c[a] < c[b]; });
        double tot = 0;
        for (int b : ids) tot += area[b];
        double acc = 0;
        size_t k = 0;
        while (k + 1 < ids.size() && acc + area[ids[k]] / 2 < tot / 2) acc += area[ids[k++]];
        if (k == 0) { acc = area[ids[0]]; k = 1; }
        double frac = acc / tot, ftot = free_area(x0, x1, y0, y1);
        double lo = vert ? x0 : y0, hi = vert ? x1 : y1;
        for (int it = 0; it < 40; ++it) {
            double m = (lo + hi) / 2;
            double f = vert ? free_area(x0, m, y0, y1) : free_area(x0, x1, y0, m);
            if (f < frac * ftot) lo = m; else hi = m;
        }
        double cut = (lo + hi) / 2;
        std::vector<int> A(ids.begin(), ids.begin() + k), B(ids.begin() + k, ids.end());
        std::vector<int> pa, ma, pb, mb;
        if (vert) {
            bisect(A, x0, cut, y0, y1, pa, ma);
            bisect(B, cut, x1, y0, y1, pb, mb);
            sp.insert(sp.end(), pa.begin(), pa.end()); sp.insert(sp.end(), pb.begin(), pb.end());
            sm.insert(sm.end(), ma.begin(), ma.end()); sm.insert(sm.end(), mb.begin(), mb.end());
        } else {  // A below B: Γ+ = B then A, Γ- = A then B
            bisect(A, x0, x1, y0, cut, pa, ma);
            bisect(B, x0, x1, cut, y1, pb, mb);
            sp.insert(sp.end(), pb.begin(), pb.end()); sp.insert(sp.end(), pa.begin(), pa.end());
            sm.insert(sm.end(), ma.begin(), ma.end()); sm.insert(sm.end(), mb.begin(), mb.end());
        }
    }

    void spread() {
        std::vector<int> ids;
        for (int i = 0; i < n; ++i) if (P.kind[i] != K_PRE) ids.push_back(i);
        seqp.clear(); seqm.clear();
        bisect(ids, 0, P.W, 0, P.H, seqp, seqm);
    }

    void run(int iters, double anc0) {
        for (int r = 0; r < 4; ++r) { solve_axis(false, nullptr, 0); solve_axis(true, nullptr, 0); }
        for (int k = 1; k <= iters; ++k) {
            spread();
            double anc = anc0 * k;
            solve_axis(false, &tx, anc);
            solve_axis(true, &ty, anc);
        }
        spread();
    }
};

// SP + soft shapes from a GlobalPlacer result (see above).
bool sp_from_gp(const Problem &P, State &S, const GlobalPlacer &G, bool grow, bool geo_prio, std::string &err) {
    const int n = P.n;
    // Soft shapes: the leaf region's aspect ratio (MIB groups: mean over members).
    S.gidx.assign(P.gmem.size(), 0);
    for (int g = 0; g < (int)P.gmem.size(); ++g) {
        double lr = 0;
        for (int b : P.gmem[g]) lr += std::log(std::max(1e-9, G.leafw[b]) / std::max(1e-9, G.leafh[b]));
        lr /= P.gmem[g].size();
        const auto &ws = P.gshape[g].widths;
        int target = (int)std::lround(std::sqrt((double)P.gshape[g].A * std::exp(lr)));
        int i = (int)(std::lower_bound(ws.begin(), ws.end(), target) - ws.begin());
        S.gidx[g] = std::min(i, (int)ws.size() - 1);
    }
    init_dims(P, S);
    // Virtual placement on the doubled grid.
    std::vector<fp::XRect> r(n);
    std::vector<int64_t> kp(n, 0), km(n, 0);
    for (int i = 0; i < (int)G.seqp.size(); ++i) kp[G.seqp[i]] = i + 1;
    for (int i = 0; i < (int)G.seqm.size(); ++i) km[G.seqm[i]] = i + 1;
    if (geo_prio) {  // tie-breaks from target geometry: Γ+ by x - y, Γ- by x + y
        for (int i = 0; i < n; ++i) {
            kp[i] = (int64_t)std::llround(G.tx[i] - G.ty[i]);
            km[i] = (int64_t)std::llround(G.tx[i] + G.ty[i]);
        }
    }
    std::vector<char> usedx(2 * P.W + 8, 0), usedy(2 * P.H + 8, 0);
    auto inside_pre = [&](int64_t cx, int64_t cy) {
        for (int p : P.pre) {
            int64_t x0 = 2LL * P.px[p], y0 = 2LL * P.py[p];
            if (cx >= x0 && cx < x0 + 2LL * P.fw[p] && cy >= y0 && cy < y0 + 2LL * P.fh[p]) return p;
        }
        return -1;
    };
    // Free odd coordinate near v in [lo, hi] (odd), not used on that axis.
    auto take = [&](std::vector<char> &u, int64_t v, int64_t lo, int64_t hi) {
        v = std::max(lo, std::min(hi, v | 1));
        for (int64_t d = 0; d <= hi - lo + 2; d += 2) {
            for (int s = 0; s < 2; ++s) {
                int64_t c = s ? v - d : v + d;
                if (c >= lo && c <= hi && !u[c + 4]) { u[c + 4] = 1; return c; }
            }
        }
        return v;
    };
    for (int p : P.pre) r[p] = {2LL * P.px[p], 2LL * P.py[p], 2LL * P.fw[p], 2LL * P.fh[p]};
    for (int i = 0; i < n; ++i) {
        if (P.kind[i] == K_PRE) continue;
        int k = P.bmask[i];
        int64_t cx = (int64_t)std::floor(2 * G.tx[i]), cy = (int64_t)std::floor(2 * G.ty[i]);
        if (!k) {
            int p = inside_pre(cx | 1, cy | 1);
            if (p >= 0) {  // push out across the nearest preplaced edge
                int64_t x0 = 2LL * P.px[p], y0 = 2LL * P.py[p], x1 = x0 + 2LL * P.fw[p], y1 = y0 + 2LL * P.fh[p];
                int64_t dl = cx - x0, dr = x1 - cx, db = cy - y0, dt = y1 - cy;
                int64_t m = std::min(std::min(dl, dr), std::min(db, dt));
                if (m == dl) cx = x0 - 1; else if (m == dr) cx = x1 + 1; else if (m == db) cy = y0 - 1; else cy = y1 + 1;
            }
            int64_t ox = take(usedx, cx, -1, 2LL * P.W + 1), oy = take(usedy, cy, -1, 2LL * P.H + 1);
            for (int t = 0; inside_pre(ox, oy) >= 0; ++t) {  // walk out of preplaced blocks
                if (t & 1) { usedx[ox + 4] = 0; ox = take(usedx, ox + (t % 4 == 1 ? 2 : -2) * (t / 2 + 1), -1, 2LL * P.W + 1); }
                else { usedy[oy + 4] = 0; oy = take(usedy, oy + (t % 4 == 0 ? 2 : -2) * (t / 2 + 1), -1, 2LL * P.H + 1); }
                if (t > 4 * (P.W + P.H)) { err = "no free cell"; return false; }
            }
            cx = ox; cy = oy;
        } else {
            cx = take(usedx, cx, 1, 2LL * P.W - 1);
            cy = take(usedy, cy, 1, 2LL * P.H - 1);
            if (k & lab3::EDGE_LEFT) cx = -4;
            if (k & lab3::EDGE_RIGHT) cx = 2LL * P.W + 3;
            if (k & lab3::EDGE_BOTTOM) cy = -4;
            if (k & lab3::EDGE_TOP) cy = 2LL * P.H + 3;
        }
        r[i] = {cx, cy, 1, 1};
    }
    if (grow) {
        // Grow each point toward its real (doubled) size around its centre,
        // shrinking overlapping pairs until the set is disjoint, so aligned
        // pairs get the relation their real extents imply.
        std::vector<double> ex(n, 0), ey(n, 0);
        std::vector<int64_t> c0(n), c1(n);
        for (int i = 0; i < n; ++i) {
            if (P.kind[i] == K_PRE) continue;
            c0[i] = r[i].x; c1[i] = r[i].y;
            int k = P.bmask[i];
            ex[i] = (k & (lab3::EDGE_LEFT | lab3::EDGE_RIGHT)) ? 0 : S.w[i];
            ey[i] = (k & (lab3::EDGE_BOTTOM | lab3::EDGE_TOP)) ? 0 : S.h[i];
            if (!(k & (lab3::EDGE_LEFT | lab3::EDGE_RIGHT))) ex[i] = std::min(ex[i], (double)std::min<int64_t>(c0[i] + 1, 2LL * P.W + 1 - c0[i]));
            if (!(k & (lab3::EDGE_BOTTOM | lab3::EDGE_TOP))) ey[i] = std::min(ey[i], (double)std::min<int64_t>(c1[i] + 1, 2LL * P.H + 1 - c1[i]));
        }
        auto mk = [&](int i) {
            int64_t a = (int64_t)ex[i], b = (int64_t)ey[i];
            return fp::XRect{c0[i] - a, c1[i] - b, 2 * a + 1, 2 * b + 1};
        };
        for (int i = 0; i < n; ++i) if (P.kind[i] != K_PRE) r[i] = mk(i);
        auto ovl = [&](const fp::XRect &A, const fp::XRect &B) {
            return A.x < B.x + B.w && B.x < A.x + A.w && A.y < B.y + B.h && B.y < A.y + A.h;
        };
        for (int round = 0; round < 200; ++round) {
            bool any = false;
            for (int a = 0; a < n; ++a) {
                for (int b = a + 1; b < n; ++b) {
                    if (!ovl(r[a], r[b])) continue;
                    any = true;
                    if (P.kind[a] != K_PRE) { ex[a] *= 0.85; ey[a] *= 0.85; if (ex[a] < 1) ex[a] = 0; if (ey[a] < 1) ey[a] = 0; r[a] = mk(a); }
                    if (P.kind[b] != K_PRE) { ex[b] *= 0.85; ey[b] *= 0.85; if (ex[b] < 1) ex[b] = 0; if (ey[b] < 1) ey[b] = 0; r[b] = mk(b); }
                }
            }
            if (!any) break;
        }
    }
    if (!fp::extract_sp(r, kp, km, S.gp, S.gm, err)) return false;
    S.pp.assign(n, 0);
    S.pm.assign(n, 0);
    S.sync_pos();
    return true;
}

double now_s(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

// ------------------------------------------------------------------ SA

enum MoveType { M_SWAP_P, M_SWAP_M, M_SWAP_B, M_INS_P, M_INS_M, M_SHAPE, M_NTYPES };

struct Undo {
    int type = 0, a = 0, b = 0, g = -1, old_idx = 0;
};

void seq_swap(std::vector<int> &seq, std::vector<int> &pos, int i, int j) {
    std::swap(seq[i], seq[j]);
    pos[seq[i]] = i;
    pos[seq[j]] = j;
}
void seq_move(std::vector<int> &seq, std::vector<int> &pos, int i, int j) {  // element at i goes to j
    if (i < j) {
        std::rotate(seq.begin() + i, seq.begin() + i + 1, seq.begin() + j + 1);
        for (int k = i; k <= j; ++k) pos[seq[k]] = k;
    } else if (j < i) {
        std::rotate(seq.begin() + j, seq.begin() + i, seq.begin() + i + 1);
        for (int k = j; k <= i; ++k) pos[seq[k]] = k;
    }
}

}  // namespace

// GLS legalizer (gls.h): it is written against gp.h's copies of these types
// in namespace fp; alias fp.cpp's own definitions instead of including gp.h.
#define FP_GP_H
namespace fp {
using ::Problem; using ::Rng; using ::BBConn; using ::BTConn;
using ::K_SOFT; using ::K_FIXED; using ::K_PRE;
}  // namespace fp
// The legalizer's types refer to fp.cpp's TU-local (anonymous-namespace) types; fine in this one TU.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsubobject-linkage"
#include "gls.h"
#pragma GCC diagnostic pop

namespace {

struct Solver {
    const Problem &P;
    Config cfg;
    Rng rng;
    Clock::time_point t0;
    double deadline;
    State S;
    Packer K;
    int ngroups;

    // incumbent
    bool have_feasible = false;
    double ttf = -1;
    int64_t best_wl2 = INT64_MAX;
    std::vector<lab3::Rect> best_rects;
    int64_t best_pen = INT64_MAX;          // fallback while infeasible
    std::vector<lab3::Rect> fallback_rects;
    double last_write = -1e9;
    bool dirty = false;
    int64_t cycles = 0, moves = 0, rejected = 0, accepted = 0, gate_checks = 0, gate_disagree = 0;

    Solver(const Problem &p, const Config &c, Clock::time_point start)
        : P(p), cfg(c), rng(c.seed), t0(start) {
        double margin = std::max(0.15, 0.02 * cfg.time);
        deadline = cfg.time - margin;
        fin_margin = margin;
        ngroups = (int)P.gmem.size();
        K.init(P.n);
    }

    // Phase E (strong). fin_margin: time the final polish may use after the
    // deadline. rival: another run of a restart portfolio; this run writes
    // the output file only while it has no feasible solution or a better one.
    // pause_at: run_strong() returns at this time if feasible (pause_hard:
    // also when infeasible), saving its post-feasible schedule in rs so a
    // later run_strong() call resumes it; if infeasible it runs on alone, to
    // ext_deadline when that is later than the deadline.
    double fin_margin = 0;
    const Solver *rival = nullptr;
    double pause_at = 1e18, ext_deadline = 0;   // ext_deadline: deadline when infeasible at pause_at
    bool pause_hard = false, paused = false, final_done = false;
    Clock::time_point pause_wall;
    struct Resume { bool valid = false; double T0 = 0, T2 = 0, t2 = 0, lambda = 1, wl_norm = 1; int seg = 0; } rs;
    bool gls_first_done = false;
    double hyb_at = -1;                    // hybrid: stop run_strong this long after the first feasible state (-1: off)
    bool hyb_go = false;                   // ... and it did
    State inc_S;                           // SP state of the packed incumbent
    int64_t inc_moves = 0;                 // move count when it was found
    bool no_gate = false;                  // relaxed instance (homotopy): no gate mirror, never writes
    Solver *tru = nullptr;                 // relaxed run: the true-instance solver (output, ttf)
    int64_t tru_base = 0;                  // ... its move count when the relaxed run started
    // Relaxed run (homotopy phase 1): offer the current SP and shape indices,
    // at TRUE dimensions, to the true solver; it becomes an output candidate
    // only if that packing is feasible (and passes the gate mirror there).
    void tru_offer(const std::string &out) {
        State T = S;
        init_dims(tru->P, T);
        tru->K.pack(tru->P, T);
        if (!tru->K.feasible()) return;
        // tru->S is left holding this state on purpose: the true solver is not
        // running while the relaxed run is, and phase 3 reloads tru->inc_S.
        tru->S.swap_in(T);
        tru->moves = tru_base + moves;   // process move count (FIRSTFEAS / inc_moves)
        tru->consider(out);
    }
    bool may_write() const {
        if (no_gate) return false;
        return !rival || !rival->have_feasible || (have_feasible && best_wl2 < rival->best_wl2);
    }

    // Strong mode: the packed incumbent (lowest packed WL so far, its rects
    // and SP) is tracked apart from the output (best_rects / best_wl2), which
    // may be an LP-polished placement of lower WL. In weak mode both coincide.
    int64_t inc_wl2 = INT64_MAX;
    std::vector<lab3::Rect> inc_rects;
    std::vector<int> inc_gp, inc_pm;
    bool pol_pending = false;              // packed incumbent not polished yet
    double last_pol = -1e9;
    int64_t pol_runs = 0, pol_wins = 0, pol_reject = 0, pol_wsp = 0, pol_wgeo = 0, pol_wre = 0, pol_wsprep = 0, pol_rounds = 0;
    double pol_ms = 0;

    // Output candidate `out` from an LP polish of `in`: kept if it has a lower
    // WL than the output, passes the gate mirror and satisfies at least the
    // boundary/MIB counts of `in`.
    bool offer(const std::vector<lab3::Rect> &in, std::vector<lab3::Rect> &out, int64_t wl2) {
        if (wl2 >= best_wl2) return false;
        GateInfo g1 = gate_check(P, out), g0 = gate_check(P, in);
        if (!g1.ok || !g0.ok || g1.bsat < g0.bsat || g1.msat < g0.msat) { ++pol_reject; return false; }
        best_rects.swap(out);
        best_wl2 = wl2;
        dirty = true;
        return true;
    }

    // LP polish of `in` on graphs g (or, g == nullptr, on its geometric graphs).
    bool polish_one(const std::vector<lab3::Rect> &in, const mcf::ConstraintGraphs *g) {
        mcf::PolishInfo info;
        std::vector<lab3::Rect> r = g ? mcf::polish_topology(P.inst, in, *g, &info) : mcf::polish(P.inst, in, &info);
        ++pol_runs;
        return info.ok && offer(in, r, info.wl2_out);
    }

    // Polish the packed incumbent on its SP graphs (+ geometric graphs).
    void polish_inc() {
        Clock::time_point c0 = Clock::now();
        mcf::ConstraintGraphs g;
        sp_graphs(inc_gp, inc_pm, g);
        bool w = polish_one(inc_rects, &g);
        pol_wsp += w;
        if (cfg.s_pol_geo > 0 && polish_one(inc_rects, nullptr)) { w = true; ++pol_wgeo; }
        pol_wins += w;
        pol_pending = false;
        pol_ms += std::chrono::duration<double, std::milli>(Clock::now() - c0).count();
        last_pol = now_s(t0);
    }

    // Exact LP WL of the current S on its SP graphs (s_lpcost), with one
    // warm-started AxisSolver per axis. K must hold the packing of S (gives
    // the boundary edges to keep). Returns 2*WL, or -1 if a solve fails.
    mcf::AxisSolver lp_ax[2];
    mcf::AxisProblem lp_pr[2];
    mcf::ConstraintGraphs lp_g;
    std::vector<lab3::Rect> lp_r;
    std::vector<int64_t> lp_x[2];
    bool lp_init = false;
    double lp_t = 0;
    int64_t lp_m0 = 0;
    int64_t lp_evals = 0, lp_fail = 0, lp_wins = 0;
    double lp_ms = 0;

    int64_t lp_eval() {
        Clock::time_point c0 = Clock::now();
        const int n = P.n;
        if (!lp_init) {
            for (int a = 0; a < 2; ++a) lp_pr[a].nets = mcf::instance_nets(P.inst, a);
            lp_r.resize(n);
            lp_init = true;
        }
        for (int i = 0; i < n; ++i) lp_r[i] = lab3::Rect{K.x[i], K.y[i], S.w[i], S.h[i]};
        sp_graphs(S.gp, S.pm, lp_g);
        lp_pr[0].edges.swap(lp_g.h);
        lp_pr[1].edges.swap(lp_g.v);
        int64_t tot = 0;
        ++lp_evals;
        for (int a = 0; a < 2; ++a) {
            mcf::AxisProblem &pr = lp_pr[a];
            pr.size.resize(n);
            for (int i = 0; i < n; ++i) pr.size[i] = a == 0 ? S.w[i] : S.h[i];
            mcf::instance_bounds(P.inst, lp_r, a, pr.lo, pr.hi);
            lp_ax[a].set_problem(pr);
            mcf::AxisResult res = lp_ax[a].solve(true);
            if (!res.ok) { ++lp_fail; tot = -1; break; }
            tot += res.obj2;
            lp_x[a].swap(res.x);
        }
        lp_ms += std::chrono::duration<double, std::milli>(Clock::now() - c0).count();
        return tot;
    }

    // Offer the placement of the last successful lp_eval() (2*WL wl2) as output.
    void lp_offer(int64_t wl2) {
        if (wl2 >= best_wl2) return;
        std::vector<lab3::Rect> in = to_rects(K, S), r(P.n);
        for (int i = 0; i < P.n; ++i) r[i] = lab3::Rect{lp_x[0][i], lp_x[1][i], S.w[i], S.h[i]};
        if (lab3::wirelength2(P.inst, r) != wl2) { ++pol_reject; return; }
        lp_wins += offer(in, r, wl2);
    }

    // Soft-shape descent on the output's own topology (final polish): per
    // shape group try nearby widths, re-solve both axes' LPs (warm), keep a
    // change that lowers the exact WL. Boundary blocks that touch their
    // edges stay pinned to them under the new sizes. Stops at time tlim.
    int64_t shp_tries = 0, shp_wins = 0;
    void shape_lp(double tlim) {
        const int n = P.n;
        const std::vector<lab3::Rect> in = best_rects;
        mcf::ConstraintGraphs g;
        if (ngroups == 0 || !mcf::build_constraint_graphs(in, g)) return;
        mcf::AxisProblem pr[2];
        pr[0].edges = mcf::transitive_reduction(n, g.h);
        pr[1].edges = mcf::transitive_reduction(n, g.v);
        GateInfo g0 = gate_check(P, in);
        if (!g0.ok) return;
        std::vector<char> bok(n, 1);   // boundary block whose edges are all touched
        for (const lab3::Boundary &bd : P.inst.boundaries) {
            const lab3::Rect &q = in[bd.block];
            bok[bd.block] = !((bd.edges & lab3::EDGE_LEFT) && q.x != 0) && !((bd.edges & lab3::EDGE_RIGHT) && q.x + q.w != P.W) &&
                            !((bd.edges & lab3::EDGE_BOTTOM) && q.y != 0) && !((bd.edges & lab3::EDGE_TOP) && q.y + q.h != P.H);
        }
        std::vector<int> gi(ngroups, -1);
        for (int k = 0; k < ngroups; ++k) {
            const std::vector<int> &ws = P.gshape[k].widths, &hs = P.gshape[k].heights;
            const lab3::Rect &q = in[P.gmem[k][0]];
            for (int t = 0; t < (int)ws.size(); ++t) if (ws[t] == q.w && hs[t] == q.h) gi[k] = t;
            for (int b : P.gmem[k]) if (in[b].w != q.w || in[b].h != q.h) gi[k] = -1;
        }
        std::vector<int64_t> wv(n), hv(n), x[2];
        for (int i = 0; i < n; ++i) { wv[i] = in[i].w; hv[i] = in[i].h; }
        mcf::AxisSolver ax[2];
        for (int a = 0; a < 2; ++a) pr[a].nets = mcf::instance_nets(P.inst, a);
        auto eval = [&](bool keep) -> int64_t {
            int64_t tot = 0;
            for (int a = 0; a < 2; ++a) {
                mcf::AxisProblem &q = pr[a];
                const std::vector<int64_t> &sz = a == 0 ? wv : hv;
                const int64_t L = a == 0 ? P.W : P.H;
                q.size = sz;
                q.lo.assign(n, 0);
                q.hi.resize(n);
                for (int i = 0; i < n; ++i) {
                    q.hi[i] = L - sz[i];
                    if (P.kind[i] == K_PRE) q.lo[i] = q.hi[i] = a == 0 ? P.px[i] : P.py[i];
                    else if (bok[i]) {
                        int k = P.bmask[i];
                        if (k & (a == 0 ? lab3::EDGE_LEFT : lab3::EDGE_BOTTOM)) q.hi[i] = 0;
                        if (k & (a == 0 ? lab3::EDGE_RIGHT : lab3::EDGE_TOP)) q.lo[i] = L - sz[i];
                    }
                    if (q.lo[i] > q.hi[i]) return -1;
                }
                ax[a].set_problem(q);
                mcf::AxisResult res = ax[a].solve(true);
                if (!res.ok) return -1;
                tot += res.obj2;
                if (keep) x[a].swap(res.x);
            }
            return tot;
        };
        int64_t cur = eval(false);
        if (cur < 0) return;
        bool imp = true;
        for (int sweep = 0; sweep < 4 && imp; ++sweep) {
            imp = false;
            for (int k = 0; k < ngroups && now_s(t0) < tlim; ++k) {
                if (gi[k] < 0) continue;
                int cnt = (int)P.gshape[k].widths.size(), o = gi[k], bt = o;
                for (int d : {-2, -1, 1, 2}) {
                    int t = o + d;
                    if (t < 0 || t >= cnt) continue;
                    for (int b : P.gmem[k]) { wv[b] = P.gshape[k].widths[t]; hv[b] = P.gshape[k].heights[t]; }
                    ++shp_tries;
                    int64_t v = eval(false);
                    if (v >= 0 && v < cur) { cur = v; bt = t; }
                }
                for (int b : P.gmem[k]) { wv[b] = P.gshape[k].widths[bt]; hv[b] = P.gshape[k].heights[bt]; }
                if (bt != o) { gi[k] = bt; imp = true; ++shp_wins; }
            }
        }
        if (eval(true) != cur) return;
        std::vector<lab3::Rect> r(n);
        for (int i = 0; i < n; ++i) r[i] = lab3::Rect{x[0][i], x[1][i], wv[i], hv[i]};
        if (lab3::wirelength2(P.inst, r) == cur) offer(in, r, cur);
    }

    // MCF / LP WL-polish hook (strong mode only). Called with final=false on
    // every new gate-checked packed incumbent (inc_*, S/K hold it), and once
    // with final=true before the final write. Polished placements replace
    // best_rects only when better (offer()); S/K are never touched. During the
    // run a polish happens at most every s_pol_ms ms; a skipped incumbent is
    // polished later from maybe_write() or at the end.
    void on_incumbent(bool final) {
        if (cfg.s_polish <= 0) return;
        if (final) {
            Clock::time_point c0 = Clock::now();
            if (pol_pending) polish_inc();
            // Re-polish the output on its own geometric graphs until no gain.
            // The final steps use at most ~60% of the time margin.
            const double margin = fin_margin, tstop = deadline + 0.6 * margin;
            for (int it = 0; it < 10 && now_s(t0) < tstop; ++it) {
                std::vector<lab3::Rect> cur = best_rects;
                if (!polish_one(cur, nullptr)) break;
                ++pol_wre;
            }
            if (cfg.s_pol_shape > 0) {
                shape_lp(deadline + 0.35 * margin);
                for (int it = 0; it < 10 && now_s(t0) < tstop; ++it) {
                    std::vector<lab3::Rect> cur = best_rects;
                    if (!polish_one(cur, nullptr)) break;
                    ++pol_wre;
                }
            }
            // Further rounds while they gain: re-polish on the SP graphs of
            // SPs extracted from the output (tie-breaks by corner and centre
            // diagonals, then random ones), soft-shape descent again,
            // geometric re-polishes.
            Rng trng(cfg.seed ^ 0x5DEECE66DULL);
            for (int round = 0; cfg.s_sp_repol > 0 && round < 5 && now_s(t0) < tstop; ++round) {
                const int64_t wl_round = best_wl2;
                for (int var = 0; var < 5 && now_s(t0) < tstop; ++var) {
                    const std::vector<lab3::Rect> cur = best_rects;
                    std::vector<fp::XRect> xr(P.n);
                    std::vector<int64_t> kp(P.n), km(P.n);
                    for (int i = 0; i < P.n; ++i) {
                        xr[i] = {cur[i].x, cur[i].y, cur[i].w, cur[i].h};
                        int64_t cx = var == 0 ? cur[i].x : 2 * cur[i].x + cur[i].w;
                        int64_t cy = var == 0 ? cur[i].y : 2 * cur[i].y + cur[i].h;
                        kp[i] = var < 2 ? cx - cy : (int64_t)trng.below(1 << 30);
                        km[i] = var < 2 ? cx + cy : (int64_t)trng.below(1 << 30);
                    }
                    std::vector<int> gp, gm, pm(P.n);
                    std::string err;
                    if (!fp::extract_sp(xr, kp, km, gp, gm, err)) continue;
                    for (int k = 0; k < P.n; ++k) pm[gm[k]] = k;
                    mcf::ConstraintGraphs g;
                    sp_graphs(gp, pm, g);
                    if (polish_one(cur, &g)) ++pol_wsprep;
                }
                if (cfg.s_pol_shape > 0) shape_lp(tstop);
                for (int it = 0; it < 10 && now_s(t0) < tstop; ++it) {
                    std::vector<lab3::Rect> cur = best_rects;
                    if (!polish_one(cur, nullptr)) break;
                    ++pol_wre;
                }
                if (best_wl2 == wl_round) break;
                ++pol_rounds;
            }
            std::fprintf(stderr, "STATS polish runs=%lld wins=%lld rejected=%lld ms=%.1f wins_sp=%lld wins_geo=%lld wins_repolish=%lld shape_tries=%lld shape_wins=%lld wins_sprepol=%lld gain_rounds=%lld packed_wl2=%lld out_wl2=%lld final_ms=%.1f\n",
                         (long long)pol_runs, (long long)pol_wins, (long long)pol_reject, pol_ms,
                         (long long)pol_wsp, (long long)pol_wgeo, (long long)pol_wre, (long long)shp_tries, (long long)shp_wins, (long long)pol_wsprep, (long long)pol_rounds,
                         (long long)inc_wl2, (long long)best_wl2,
                         std::chrono::duration<double, std::milli>(Clock::now() - c0).count());
            return;
        }
        inc_gp = S.gp;
        inc_pm = S.pm;
        pol_pending = true;
        if (pol_due(now_s(t0))) polish_inc();
    }

    bool pol_due(double t) const {
        return (t - last_pol) * 1000.0 >= cfg.s_pol_ms && pol_ms <= cfg.s_pol_budget * 1000.0 * t;
    }

    void consider(const std::string &out) {
        if (K.feasible()) {
            if (K.wl2 < inc_wl2) {
                std::vector<lab3::Rect> r = to_rects(K, S);
                GateInfo gi;
                if (no_gate) gi.ok = true;
                else gi = gate_check(P, r);
                ++gate_checks;
                if (!gi.ok) {
                    ++gate_disagree;
                    std::fprintf(stderr, "WARN internal gate rejected packed state: %s\n", gi.why.c_str());
                    return;
                }
                if (cfg.strong && cfg.s_wl0 > 0 && have_feasible) return;  // sanity run: keep the first feasible state
                inc_wl2 = K.wl2;
                if (tru) tru_offer(out);
                if (cfg.strong) {
                    inc_S = S;
                    inc_moves = moves;
                    if (K.wl2 < best_wl2) { best_wl2 = K.wl2; best_rects = r; }
                    inc_rects.swap(r);
                    on_incumbent(false);
                } else {
                    best_wl2 = K.wl2;
                    best_rects.swap(r);
                }
                if (!have_feasible) {
                    have_feasible = true;
                    ttf = now_s(t0);
                    std::fprintf(stderr, no_gate ? "RELAXFEAS moves=%lld wl2=%lld\n" : "FIRSTFEAS moves=%lld wl2=%lld\n", (long long)moves, (long long)best_wl2);
                    if (may_write()) write_output(P, best_rects, out);
                    last_write = ttf;
                    dirty = false;
                    return;
                }
                dirty = true;
            }
        } else if (!have_feasible) {
            int64_t pen = K.ovfx + K.ovfy + K.viol;
            if (pen < best_pen) {
                best_pen = pen;
                fallback_rects = to_rects(K, S);
                dirty = true;
            }
        }
    }

    // Phase F: legalize placement r0 (soft shapes: width index g0 per shape
    // group) with the GLS legalizer (gls.h; it clamps the blocks into the
    // outline). On success S becomes the SP extracted from the legal placement
    // (keys X2 - Y2 / X2 + Y2 on doubled centroids), which must keep the
    // invariant, add no forbidden relation and pack feasible (gate mirror
    // too); consider() then records it (ttf / FIRSTFEAS) and the next
    // post-feasible start uses s_gls_t2 (gls_handoff). On failure S is unchanged.
    bool gls_handoff = false;
    int64_t gls_calls = 0, gls_wins = 0, gls_plat_calls = 0;
    double gls_secs = 0;
    bool gls_try(const std::vector<lab3::Rect> &r0, const std::vector<int> &g0, double budget, const std::string &out,
                 const char *tag) {
        if (budget < 0.2) return false;
        std::vector<lab3::Rect> r = r0;
        std::vector<int> g = g0;
        Rng grng(cfg.seed * 0x9E3779B97F4A7C15ULL + 12345 + 7919ULL * (uint64_t)gls_calls);
        fp::GlsStats st;
        const double ts = now_s(t0);
        ++gls_calls;
        bool legal = fp::gls_legalize(P, r, g, budget, grng, st), ok = false;
        int d_ext = -1, d_inv = -1, d_feas = -1;
        int64_t d_forb = -1;
        if (legal) {
            State T;
            T.gidx = g;
            init_dims(P, T);
            std::vector<fp::XRect> xr(P.n);
            std::vector<int64_t> kp(P.n), km(P.n);
            for (int i = 0; i < P.n; ++i) {
                xr[i] = {r[i].x, r[i].y, r[i].w, r[i].h};
                int64_t X2 = 2 * r[i].x + r[i].w, Y2 = 2 * r[i].y + r[i].h;
                kp[i] = X2 - Y2;
                km[i] = X2 + Y2;
            }
            std::string err;
            d_ext = fp::extract_sp(xr, kp, km, T.gp, T.gm, err) ? 1 : 0;
            if (d_ext) {
                T.pp.assign(P.n, 0);
                T.pm.assign(P.n, 0);
                T.sync_pos();
                d_inv = invariant_ok(P, T) ? 1 : 0;
                d_forb = forb_count(P, T);
                if (d_inv && d_forb == 0) {
                    K.pack(P, T);
                    d_feas = K.feasible() ? 1 : 0;
                    if (d_feas && (no_gate || gate_check(P, to_rects(K, T)).ok)) { S = T; ok = true; }
                }
            }
        }
        gls_secs += now_s(t0) - ts;
        std::fprintf(stderr, "STATS gls %s legal=%d ok=%d ext=%d inv=%d forb=%lld feas=%d t=%.2f budget=%.2f steps=%lld at=%.2f\n", tag,
                     legal ? 1 : 0, ok ? 1 : 0, d_ext, d_inv, (long long)d_forb, d_feas, now_s(t0) - ts, budget, (long long)st.steps,
                     now_s(t0));
        if (!ok) {   // K may hold the rejected packing; restore it for a valid S (none in standalone homotopy)
            if (!S.gp.empty()) K.pack(P, S);
            return false;
        }
        ++gls_wins;
        gls_handoff = true;
        consider(out);
        return true;
    }

    // Current (accepted) state's coordinates and penalty, for targeted moves.
    std::vector<int> cur_x, cur_y;
    int64_t cur_ox = 0, cur_oy = 0, cur_pen = 0;
    bool strong_inv = false;   // strong: also never add a forbidden relation
    int64_t cur_forb = 0;
    std::vector<std::pair<int, int>> culprits;

    bool shape_move(Undo &u, int g, int dir) {  // dir: -1 narrower, +1 wider, 0 any
        int cnt = (int)P.gshape[g].widths.size();
        int cur = S.gidx[g], nw;
        if (cnt <= 1) return false;
        if (dir < 0) { if (cur == 0) return false; nw = rng.below(cur); }
        else if (dir > 0) { if (cur == cnt - 1) return false; nw = cur + 1 + rng.below(cnt - 1 - cur); }
        else if (rng.uni() < 0.5) {
            int step = std::max(1, cnt / 20);
            nw = cur + rng.below(2 * step + 1) - step;
            nw = std::max(0, std::min(cnt - 1, nw));
            if (nw == cur) nw = (cur + 1 < cnt) ? cur + 1 : cur - 1;
        } else {
            nw = rng.below(cnt);
            if (nw == cur) return false;
        }
        u.type = M_SHAPE;
        u.g = g;
        u.old_idx = cur;
        S.gidx[g] = nw;
        apply_group(P, S, g);
        return true;
    }

    // SP move of the given type on blocks a, b; b < 0 picks a partner.
    bool sp_move(Undo &u, int type, int a, int b) {
        const int n = P.n;
        if (b < 0) {
            if (rng.uni() < cfg.p_local) {
                const std::vector<int> &seq = (type == M_SWAP_M || type == M_INS_M) ? S.gm : S.gp;
                const std::vector<int> &pos = (type == M_SWAP_M || type == M_INS_M) ? S.pm : S.pp;
                int d = 1 + rng.below(std::max(2, n / 16));
                int j = pos[a] + (rng.uni() < 0.5 ? -d : d);
                if (j < 0 || j >= n) j = pos[a] + (j < 0 ? d : -d);
                if (j < 0 || j >= n) return false;
                b = seq[j];
            } else {
                b = rng.below(n - 1);
                if (b >= a) ++b;
            }
        }
        u.type = type;
        switch (type) {
        case M_SWAP_P: { int i = S.pp[a], j = S.pp[b]; u.a = i; u.b = j; seq_swap(S.gp, S.pp, i, j); break; }
        case M_SWAP_M: { int i = S.pm[a], j = S.pm[b]; u.a = i; u.b = j; seq_swap(S.gm, S.pm, i, j); break; }
        case M_SWAP_B: {
            u.a = a; u.b = b;
            int i = S.pp[a], j = S.pp[b]; seq_swap(S.gp, S.pp, i, j);
            i = S.pm[a]; j = S.pm[b]; seq_swap(S.gm, S.pm, i, j);
            break;
        }
        case M_INS_P: { int i = S.pp[a], j = S.pp[b]; u.a = i; u.b = j; seq_move(S.gp, S.pp, i, j); break; }
        case M_INS_M: { int i = S.pm[a], j = S.pm[b]; u.a = i; u.b = j; seq_move(S.gm, S.pm, i, j); break; }
        }
        if (!invariant_ok(P, S) || (strong_inv && forb_count(P, S) > cur_forb)) { undo(u); return false; }
        return true;
    }

    // Move block a in one sequence to just after (after=true) or just before p.
    bool place_next_to(Undo &u, bool in_p, int a, int p, bool after) {
        std::vector<int> &seq = in_p ? S.gp : S.gm;
        std::vector<int> &pos = in_p ? S.pp : S.pm;
        int i = pos[a], q = pos[p], j;
        if (after) j = (i < q) ? q : q + 1;
        else j = (i < q) ? q - 1 : q;
        if (j == i) return false;
        u.type = in_p ? M_INS_P : M_INS_M;
        u.a = i;
        u.b = j;
        seq_move(seq, pos, i, j);
        if (!invariant_ok(P, S) || (strong_inv && forb_count(P, S) > cur_forb)) { undo(u); return false; }
        return true;
    }

    int random_sp_type() {
        double q = rng.uni();
        return q < 0.25 ? M_SWAP_P : q < 0.45 ? M_SWAP_M : q < 0.70 ? M_SWAP_B : q < 0.85 ? M_INS_P : M_INS_M;
    }

    // Repair move on a block that causes the current infeasibility: an
    // outline-overflowing block, or a predecessor of a preplaced block p that
    // ends past p's pinned coordinate (moved to the other relation with p).
    bool propose_target(Undo &u) {
        culprits.clear();
        for (int b = 0; b < P.n; ++b)
            if (cur_x[b] + S.w[b] > P.W || cur_y[b] + S.h[b] > P.H) culprits.push_back({b, -1});
        for (int p : P.pre) {
            int x0 = P.px[p], y0 = P.py[p], qp = S.pp[p], qm = S.pm[p];
            for (int a = 0; a < P.n; ++a) {
                if (S.pm[a] >= qm || a == p) continue;
                if (S.pp[a] < qp) { if (cur_x[a] + S.w[a] > x0) culprits.push_back({a, p}); }
                else if (cur_y[a] + S.h[a] > y0) culprits.push_back({a, p});
            }
        }
        if (culprits.empty()) return false;
        auto c = culprits[rng.below((int)culprits.size())];
        int a = c.first, p = c.second;
        if (p < 0) {
            int g = P.group[a];
            if (g >= 0 && rng.uni() < 0.5) {
                bool ox = cur_x[a] + S.w[a] > P.W, oy = cur_y[a] + S.h[a] > P.H;
                int dir = (ox && !oy) ? -1 : (oy && !ox) ? 1 : 0;
                return shape_move(u, g, dir);
            }
            return sp_move(u, random_sp_type(), a, -1);
        }
        bool left = S.pp[a] < S.pp[p];  // a left-of p (else a below p)
        if (rng.uni() < 0.5) return place_next_to(u, false, a, p, true);  // Γ-: a after p
        return place_next_to(u, true, a, p, left);  // Γ+: a after p (left case) / before p (below case)
    }

    // Applies a random move; returns false if it is rejected (state unchanged).
    bool propose(Undo &u) {
        if (cur_pen > 0 && rng.uni() < cfg.p_target) return propose_target(u);
        if (ngroups > 0 && rng.uni() < cfg.p_shape) {
            int g = rng.below(ngroups);
            int dir = 0;
            if ((cur_ox > 0 || cur_oy > 0) && rng.uni() < 0.4) dir = cur_ox >= cur_oy ? -1 : 1;  // fit
            return shape_move(u, g, dir);
        }
        return sp_move(u, random_sp_type(), rng.below(P.n), -1);
    }

    // ------------------------------------------------ strong-mode moves
    // Insert block a into Γ+ (in_p) or Γ- at a position drawn from the window
    // of positions that keeps the boundary invariant and the preplaced-pair
    // geometry (the relation of a with every other block x depends only on
    // whether a lands before or after x, so the valid final indices form an
    // interval [L, U]). Never rejected for the invariant.
    bool window(bool in_p, int a, int &L, int &U) {
        const int n = P.n;
        const std::vector<int> &pos = in_p ? S.pp : S.pm;
        const std::vector<int> &opos = in_p ? S.pm : S.pp;
        const int i = pos[a], ka = P.bmask[a];
        const bool apre = P.kind[a] == K_PRE;
        const int np = (int)P.pre.size();
        L = 0;
        U = n - 1;
        for (int x = 0; x < n; ++x) {
            if (x == a) continue;
            int kx = P.bmask[x];
            bool xpre = apre && P.kind[x] == K_PRE;
            int fm = 0;  // forbidden rel(a, x) bits
            if (strong_inv && np) {
                if (!apre && P.kind[x] == K_PRE) fm = P.fmask[(size_t)a * np + P.preidx[x]];
                else if (apre && P.kind[x] != K_PRE) {
                    int g = P.fmask[(size_t)x * np + P.preidx[a]];  // rel(x, a): swap 0<->1, 2<->3
                    fm = ((g & 1) << 1) | ((g & 2) >> 1) | ((g & 4) << 1) | ((g & 8) >> 1);
                }
            }
            if (!kx && !ka && !xpre && !fm) continue;
            bool o = opos[a] < opos[x];  // a before x in the other sequence
            // Relation if a lands before x / after x in this sequence:
            // 0 = a left-of x, 1 = x left-of a, 2 = a below x, 3 = x below a.
            int rb, ra;
            if (in_p) { rb = o ? 0 : 3; ra = o ? 2 : 1; }
            else { rb = o ? 0 : 2; ra = o ? 3 : 1; }
            auto bad = [&](int rel) {
                if ((fm >> rel) & 1) return true;
                switch (rel) {
                case 0: if ((kx & lab3::EDGE_LEFT) || (ka & lab3::EDGE_RIGHT)) return true;
                        return xpre && P.px[a] + P.fw[a] > P.px[x];
                case 1: if ((ka & lab3::EDGE_LEFT) || (kx & lab3::EDGE_RIGHT)) return true;
                        return xpre && P.px[x] + P.fw[x] > P.px[a];
                case 2: if ((kx & lab3::EDGE_BOTTOM) || (ka & lab3::EDGE_TOP)) return true;
                        return xpre && P.py[a] + P.fh[a] > P.py[x];
                default: if ((ka & lab3::EDGE_BOTTOM) || (kx & lab3::EDGE_TOP)) return true;
                         return xpre && P.py[x] + P.fh[x] > P.py[a];
                }
            };
            int r = pos[x] - (pos[x] > i);  // rank of x once a is removed
            if (bad(rb)) L = std::max(L, r + 1);  // a must land after x
            if (bad(ra)) U = std::min(U, r);      // a must land before x
        }
        return L < U;
    }

    bool ins_window(Undo &u, bool in_p, int a) {
        const int n = P.n;
        std::vector<int> &seq = in_p ? S.gp : S.gm;
        std::vector<int> &pos = in_p ? S.pp : S.pm;
        const int i = pos[a];
        int L, U;
        if (!window(in_p, a, L, U)) return false;
        int j;
        if (rng.uni() < cfg.p_local) {
            int d = 1 + rng.below(std::max(2, n / 16));
            j = i + (rng.uni() < 0.5 ? -d : d);
            j = std::max(L, std::min(U, j));
        } else {
            j = L + rng.below(U - L + 1);
        }
        if (j == i) return false;
        u.type = in_p ? M_INS_P : M_INS_M;
        u.a = i;
        u.b = j;
        seq_move(seq, pos, i, j);
        return true;
    }

    // Steepest descent on the penalty over soft shapes only: for every shape
    // group try nearby and spread-out widths, apply the best, repeat.
    void shape_polish(const std::string &out) {
        K.pack(P, S);
        int64_t cur = K.ovfx + K.ovfy + K.viol;
        for (int it = 0; it < 50 && cur > 0 && now_s(t0) < deadline; ++it) {
            int64_t bp = cur;
            int bg = -1, bi = 0;
            for (int g = 0; g < ngroups; ++g) {
                int cnt = (int)P.gshape[g].widths.size(), o = S.gidx[g];
                for (int k = 0; k < 24; ++k) {
                    int nw = k < 16 ? o + (k < 8 ? -(k + 1) : k - 7) : (int)((k - 16) * (cnt - 1) / 7.0);
                    if (nw < 0 || nw >= cnt || nw == o) continue;
                    S.gidx[g] = nw;
                    apply_group(P, S, g);
                    K.pack(P, S);
                    int64_t pe = K.ovfx + K.ovfy + K.viol;
                    if (pe < bp) { bp = pe; bg = g; bi = nw; }
                }
                S.gidx[g] = o;
                apply_group(P, S, g);
            }
            if (bg < 0) break;
            S.gidx[bg] = bi;
            apply_group(P, S, bg);
            cur = bp;
        }
        K.pack(P, S);
        consider(out);
    }

    // Phase G directed moves (phase 2 only, knob-gated). Net-median insertion:
    // block a moves next to (one sequence) the block whose centre is nearest
    // the weighted median of a's connection partners (current placement).
    bool in_phase2 = false;
    std::vector<std::vector<std::pair<int, int>>> nbr;   // per block: (partner, weight), partner < 0: terminal ~index
    std::vector<std::pair<int, int>> mx_, my_;
    bool med_move(Undo &u) {
        const int n = P.n;
        if (nbr.empty()) {
            nbr.assign(n, {});
            for (const BBConn &c : P.bb) { nbr[c.a].push_back({c.b, c.w}); nbr[c.b].push_back({c.a, c.w}); }
            for (int e = 0; e < (int)P.bt.size(); ++e) nbr[P.bt[e].a].push_back({-1 - e, P.bt[e].w});
        }
        int a = rng.below(n);
        if (P.kind[a] == K_PRE || nbr[a].empty()) return false;
        mx_.clear(); my_.clear();
        int64_t tw = 0;
        for (auto &e : nbr[a]) {
            int x2, y2;
            if (e.first >= 0) { x2 = 2 * cur_x[e.first] + S.w[e.first]; y2 = 2 * cur_y[e.first] + S.h[e.first]; }
            else { x2 = P.bt[-1 - e.first].tx2; y2 = P.bt[-1 - e.first].ty2; }
            mx_.push_back({x2, e.second}); my_.push_back({y2, e.second});
            tw += e.second;
        }
        auto wmed = [&](std::vector<std::pair<int, int>> &v) {
            std::sort(v.begin(), v.end());
            int64_t acc = 0;
            for (auto &q : v) { acc += q.second; if (2 * acc >= tw) return q.first; }
            return v.back().first;
        };
        const int64_t mxv = wmed(mx_), myv = wmed(my_);
        int c = -1;
        int64_t bd = INT64_MAX;
        for (int b = 0; b < n; ++b) {
            if (b == a) continue;
            int64_t d = std::llabs(2LL * cur_x[b] + S.w[b] - mxv) + std::llabs(2LL * cur_y[b] + S.h[b] - myv);
            if (d < bd) { bd = d; c = b; }
        }
        if (c < 0) return false;
        return place_next_to(u, rng.uni() < 0.5, a, c, rng.uni() < 0.5);
    }
    bool mib_move(Undo &u) {
        const auto &M = P.inst.mibs;
        if (M.empty()) return false;
        const auto &g = M[rng.below((int)M.size())];
        if (g.size() < 2) return false;
        int i = rng.below((int)g.size()), j = rng.below((int)g.size() - 1);
        if (j >= i) ++j;
        return sp_move(u, M_SWAP_B, g[i], g[j]);
    }

    bool propose_strong(Undo &u) {
        if (in_phase2 && cfg.s_p_med > 0 && rng.uni() < cfg.s_p_med) return med_move(u);
        if (in_phase2 && cfg.s_p_mib > 0 && rng.uni() < cfg.s_p_mib) return mib_move(u);
        if (cur_pen > 0 && rng.uni() < cfg.p_target) return propose_target(u);
        if (ngroups > 0 && rng.uni() < cfg.p_shape) {
            int g = rng.below(ngroups);
            int dir = 0;
            if ((cur_ox > 0 || cur_oy > 0) && rng.uni() < 0.4) dir = cur_ox >= cur_oy ? -1 : 1;
            return shape_move(u, g, dir);
        }
        double q = rng.uni();
        int a = rng.below(P.n);
        if (q < cfg.s_p_ins) return ins_window(u, q < cfg.s_p_ins / 2, a);
        return sp_move(u, random_sp_type(), a, -1);
    }

    void undo(const Undo &u) {
        switch (u.type) {
        case M_SHAPE: S.gidx[u.g] = u.old_idx; apply_group(P, S, u.g); break;
        case M_SWAP_P: seq_swap(S.gp, S.pp, u.a, u.b); break;
        case M_SWAP_M: seq_swap(S.gm, S.pm, u.a, u.b); break;
        case M_SWAP_B: {
            int i = S.pp[u.a], j = S.pp[u.b]; seq_swap(S.gp, S.pp, i, j);
            i = S.pm[u.a]; j = S.pm[u.b]; seq_swap(S.gm, S.pm, i, j);
            break;
        }
        case M_INS_P: seq_move(S.gp, S.pp, u.b, u.a); break;
        case M_INS_M: seq_move(S.gm, S.pm, u.b, u.a); break;
        }
    }

    void maybe_write(const std::string &out, double t) {
        if (tru) tru->maybe_write(out, t);   // relaxed run: keep the true output file fresh
        if (pol_pending && pol_due(t)) polish_inc();
        if (dirty && t - last_write >= 2.0) {
            if (may_write()) write_output(P, have_feasible ? best_rects : fallback_rects, out);
            last_write = t;
            dirty = false;
        }
    }

    // Write the current output now (homotopy switch and phase boundaries), so
    // a kill at any instant after FIRSTFEAS leaves a complete feasible file.
    void flush_output(const std::string &out) {
        if (pol_pending) polish_inc();
        if (have_feasible && may_write()) write_output(P, best_rects, out);
        last_write = now_s(t0);
        dirty = false;
    }

    void run(const std::string &out) {
        K.pack(P, S);
        consider(out);
        const double wl_norm = std::max<double>(1.0, (double)K.wl2);
        const double pen_norm = (double)P.W + P.H;
        double lambda = 1.0;
        auto cost_of = [&](int64_t wl2, int64_t pen) { return wl2 / wl_norm + lambda * pen / pen_norm; };
        int64_t cur_wl2 = K.wl2;
        cur_pen = K.ovfx + K.ovfy + K.viol;
        cur_ox = K.ovfx;
        cur_oy = K.ovfy;
        cur_x = K.x;
        cur_y = K.y;
        double cur_cost = cost_of(cur_wl2, cur_pen);

        // Initial temperature from the mean uphill delta of random moves.
        double sum_up = 0;
        int n_up = 0;
        for (int t = 0; t < 400; ++t) {
            Undo u;
            if (!propose(u)) continue;
            K.pack(P, S);
            double d = cost_of(K.wl2, K.ovfx + K.ovfy + K.viol) - cur_cost;
            if (d > 0) { sum_up += d; ++n_up; }
            undo(u);
        }
        K.pack(P, S);  // restore Packer arrays for the current state
        const double T0 = n_up ? (sum_up / n_up) / std::log(1.0 / cfg.acc0) : 1e-3;
        const double Tend = T0 * cfg.tend;
        double T = T0;
        // Schedule. Until the first gate-feasible state: annealing cycles of a
        // fixed move count (so time-to-feasible does not depend on --time),
        // each reheating to reheat*T0 from the current state. After that: one
        // exponential anneal from max(T, t2_floor*T0) to Tend over the remaining
        // time.
        const int64_t cycle_len = std::max<int64_t>(200000, (int64_t)(cfg.cycle_k * P.n));
        int64_t cycle_pos = 0;
        double Tc = T0;
        bool phase2 = false;
        double T2 = 0, t2 = 0;

        int64_t win_moves = 0, win_feas = 0;
        const int64_t WIN = 1000;
        double t = now_s(t0);
        const double t_start = t;
        for (;;) {
            if ((moves & 63) == 0) {
                t = now_s(t0);
                if (t >= deadline) break;
                if (!phase2 && have_feasible) {
                    phase2 = true;
                    // A feasible start (e.g. --init-from) anneals from t2_floor*T0.
                    T2 = moves == 0 ? cfg.t2_floor * T0 : std::max(T, cfg.t2_floor * T0);
                    T2 = std::min(T2, cfg.t2_cap * T0);
                    t2 = t;
                }
                if (phase2) {
                    double f = (t - t2) / std::max(1e-3, deadline - t2);
                    T = T2 * std::pow(Tend / T2, f);
                } else {
                    if (cycle_pos >= cycle_len) { cycle_pos = 0; Tc = cfg.reheat * T0; ++cycles; }
                    T = Tc * std::pow(1e-3, (double)cycle_pos / cycle_len);
                }
                maybe_write(out, t);
            }
            ++cycle_pos;
            ++moves;
            Undo u;
            if (!propose(u)) { ++rejected; continue; }
            K.pack(P, S);
            int64_t pen = K.ovfx + K.ovfy + K.viol;
            double c = cost_of(K.wl2, pen);
            double d = c - cur_cost;
            if (d <= 0 || rng.uni() < std::exp(-d / T)) {
                cur_cost = c; cur_wl2 = K.wl2; cur_pen = pen; cur_ox = K.ovfx; cur_oy = K.ovfy;
                if (pen > 0) { cur_x = K.x; cur_y = K.y; }
                ++accepted;
                consider(out);
            } else {
                undo(u);
            }
            ++win_moves;
            if (cur_pen == 0) ++win_feas;
            if (win_moves >= WIN) {
                double fr = (double)win_feas / win_moves;
                double old = lambda;
                if (fr < 0.5) lambda = std::min(cfg.lam_max, lambda * 1.15);
                else if (fr > 0.9) lambda = std::max(0.05, lambda / 1.02);
                if (lambda != old) cur_cost = cost_of(cur_wl2, cur_pen);
                win_moves = win_feas = 0;
            }
        }
        std::fprintf(stderr, "STATS rej_boundary=%lld rej_preplaced_pair=%lld\n", (long long)g_rej_boundary, (long long)g_rej_pre);
        std::fprintf(stderr, "STATS cycles=%lld T0=%.3g lambda=%.3g accepted=%lld rejected=%lld gate_checks=%lld gate_disagree=%lld time=%.2f moves_per_s=%.0f\n",
                     (long long)cycles, T0, lambda, (long long)accepted, (long long)rejected, (long long)gate_checks,
                     (long long)gate_disagree, now_s(t0), moves / std::max(1e-9, now_s(t0) - t_start));
    }

    // ---------------------------------------------------------------- strong
    // Restart pool (analytical seeds, best first); random SPs once exhausted.
    std::vector<State> pool;
    size_t pool_next = 0;
    int64_t restarts = 0;

    bool next_start(State &s) {
        if (pool_next < pool.size()) { s = pool[pool_next++]; return true; }
        std::string err;
        State r = s;
        if (!random_initial_sp(P, r, rng, err)) return false;
        s = r;
        return true;
    }

    // Like run() with these changes (docs in the file header):
    //   * moves from propose_strong(): never add a forbidden block-preplaced
    //     relation (strong_inv), optional valid-window inserts;
    //   * an analytical seed starts at s_tstart * the probed temperature;
    //   * a penalty plateau (s_plateau_k * n moves without a better penalty
    //     in the current attempt) restarts from the next start in the pool
    //     (s_plat_restart; lambda reset, T = s_restart_t * re-probed T), or,
    //     with s_plat_restart=0, reheats to a T re-probed under the current
    //     lambda for a short s_kick_cycle_k * n cycle.
    void run_strong(const std::string &out, bool seeded) {
        cfg.p_target = cfg.s_p_target;
        strong_inv = cfg.s_forbid > 0;
        K.pack(P, S);
        consider(out);
        if (seeded && !rs.valid && !have_feasible && !gls_first_done && cfg.s_gls > 0 && cfg.s_wl0 <= 0) {
            // Phase F: GLS-first from the start state's packed placement.
            gls_first_done = true;
            gls_try(to_rects(K, S), S.gidx, std::min(cfg.s_gls_t, cfg.s_gls_frac * (deadline - now_s(t0))), out, "first");
        }
        in_phase2 = false;
        const bool resume = rs.valid;
        rs.valid = false;
        paused = false;
        const double wl_norm = resume ? rs.wl_norm : std::max<double>(1.0, (double)K.wl2);
        const double wl_w = cfg.s_wl0 > 0 ? 0.0 : 1.0;
        const double pen_norm = (double)P.W + P.H;
        double lambda = resume ? rs.lambda : 1.0;
        auto cost_of = [&](int64_t wl2, int64_t pen) { return wl_w * wl2 / wl_norm + lambda * pen / pen_norm; };
        int64_t cur_wl2 = 0;
        double cur_cost = 0;
        auto load_cur = [&]() {
            K.pack(P, S);
            cur_forb = forb_count(P, S);
            cur_wl2 = K.wl2;
            cur_pen = K.ovfx + K.ovfy + K.viol;
            cur_ox = K.ovfx;
            cur_oy = K.ovfy;
            cur_x = K.x;
            cur_y = K.y;
            cur_cost = cost_of(cur_wl2, cur_pen);
        };
        // Mean uphill cost delta of 400 random moves -> T with acceptance acc0.
        auto probe = [&]() {
            double sum_up = 0;
            int n_up = 0;
            for (int t = 0; t < 400; ++t) {
                Undo u;
                if (!propose_strong(u)) continue;
                K.pack(P, S);
                double d = cost_of(K.wl2, K.ovfx + K.ovfy + K.viol) - cur_cost;
                if (d > 0) { sum_up += d; ++n_up; }
                undo(u);
            }
            K.pack(P, S);
            return n_up ? (sum_up / n_up) / std::log(1.0 / cfg.acc0) : 1e-3;
        };
        load_cur();
        const double T0 = resume ? rs.T0 : probe();
        const double Tend = T0 * cfg.s_tend;
        double T = T0;
        const int64_t cycle_len = std::max<int64_t>(200000, (int64_t)(cfg.s_cycle_k * P.n));
        const int64_t kick_len = std::max<int64_t>(20000, (int64_t)(cfg.s_kick_cycle_k * P.n));
        const int64_t plateau = std::max<int64_t>(20000, (int64_t)(cfg.s_plateau_k * P.n));
        int64_t cycle_pos = 0, cur_cycle = cycle_len, last_imp = 0;
        double Tc = seeded ? cfg.s_tstart * T0 : T0;
        bool phase2 = resume, lpmode = false;
        double T2 = resume ? rs.T2 : 0, t2 = resume ? rs.t2 : 0;
        // Post-feasible segments (s_p2_cyc) and stall re-anneals (s_stall_k).
        const int nseg = std::max(1, (int)cfg.s_p2_cyc);
        int seg = resume ? rs.seg : 0;
        double T2s = T2;   // start temperature of the current segment
        // Phase G acceptance rules (phase 2 only; not in the relaxed homotopy run).
        const int acc_rule = no_gate ? 0 : (int)cfg.s_acc;
        std::vector<double> lahc;
        size_t lahc_i = 0;
        double lam_T = 0, lam_rate = 0, lam_f = 0;
        int64_t acc_n = 0, acc_acc = 0, up_n = 0, up_acc = 0;
        double up_f[10] = {0}, up_fn[10] = {0};   // uphill acceptance per tenth of phase 2 (stats)
        const int64_t stall = cfg.s_stall_k > 0 ? std::max<int64_t>(20000, (int64_t)(cfg.s_stall_k * P.n)) : INT64_MAX;
        int64_t stalls = 0;
        auto to_incumbent = [&]() {
            if (cfg.s_p2_best > 0 && !inc_S.gp.empty()) { S = inc_S; load_cur(); }
        };
        int64_t run_best = cur_pen;  // best penalty of the current attempt
        State snap = S;              // ... and its state
        int kicks = 0;
        int64_t kick_best = INT64_MAX;

        int64_t win_moves = 0, win_feas = 0;
        const int64_t WIN = 1000;
        double t = now_s(t0);
        const double t_start = t;
        double trace_next = 0;
        for (;;) {
            if ((moves & 63) == 0) {
                t = now_s(t0);
                if (cfg.s_trace > 0 && t >= trace_next) {   // diagnostic: incumbent WL over time
                    trace_next = (std::floor(t / cfg.s_trace) + 1) * cfg.s_trace;
                    std::fprintf(stderr, "TRACE t=%.2f rel=%d inc_wl2=%lld best_wl2=%lld\n", t, no_gate ? 1 : 0,
                                 inc_wl2 == INT64_MAX ? -1LL : (long long)inc_wl2, best_wl2 == INT64_MAX ? -1LL : (long long)best_wl2);
                }
                if (t >= pause_at) {
                    if (have_feasible || pause_hard) {
                        rs.valid = phase2;
                        rs.T0 = T0; rs.T2 = T2; rs.t2 = t2; rs.lambda = lambda; rs.wl_norm = wl_norm; rs.seg = seg;
                        paused = true;
                        pause_wall = Clock::now();
                        break;
                    }
                    pause_at = 1e18;
                    deadline = std::max(deadline, ext_deadline);
                }
                if (t >= deadline) break;
                if (hyb_at >= 0 && have_feasible && t >= ttf + hyb_at) {
                    // Phase F hybrid: hand the remaining budget to the homotopy (run_homotopy).
                    hyb_at = -1;
                    if (deadline - t >= cfg.s_hyb_min * deadline) { hyb_go = true; break; }
                }
                if (cfg.s_wl0 > 0 && have_feasible) break;
                if (!phase2 && have_feasible) {
                    phase2 = true;
                    T2 = moves == 0 ? cfg.t2_floor * T0 : std::max(T, cfg.t2_floor * T0);
                    T2 = std::min(T2, cfg.t2_cap * T0);
                    if (gls_handoff) T2 = cfg.s_gls_t2 * T0;   // Phase F: start temperature after a hand-off
                    gls_handoff = false;
                    t2 = t;
                    if (acc_rule == 1) lahc.assign(std::max<size_t>(1, (size_t)cfg.s_lahc_len), cur_cost);
                    in_phase2 = !no_gate && (cfg.s_p_med > 0 || cfg.s_p_mib > 0);
                    if (in_phase2) { K.pack(P, S); cur_x = K.x; cur_y = K.y; }
                    lam_T = T2;
                    lam_rate = 0.5 * cfg.s_lam_a0;
                }
                if (phase2) {
                    double f = (t - t2) / std::max(1e-3, deadline - t2);
                    int k = std::min(nseg - 1, (int)(f * nseg));
                    if (k > seg) {   // new segment: reheat (from the best packed incumbent)
                        seg = k;
                        to_incumbent();
                    }
                    if (nseg == 1 && moves - std::max(inc_moves, last_imp) > stall) {
                        // Stall: re-anneal from the incumbent over the remaining time.
                        to_incumbent();
                        T2 = std::min(T2, cfg.s_stall_t * T);
                        t2 = t;
                        f = 0;
                        last_imp = moves;
                        ++stalls;
                    }
                    T2s = seg == 0 ? T2 : cfg.s_p2_rh * T2;
                    double g = nseg == 1 ? f : std::min(1.0, f * nseg - seg);
                    T = T2s * std::pow(Tend / T2s, cfg.s_cool_p == 1 ? g : std::pow(g, cfg.s_cool_p));
                    lam_f = f;
                    if (acc_rule == 2) T = lam_T;
                    if (cfg.s_quench > 0 && f >= 1.0 - cfg.s_quench) T = 1e-300;   // quench: greedy finish
                    if (!lpmode && cfg.s_lpcost > 0 && f >= 1.0 - cfg.s_lp_frac) {
                        // Switch the WL term to the LP WL; re-cost the current state.
                        lpmode = true;
                        lp_t = t;
                        lp_m0 = moves;
                        K.pack(P, S);
                        if (cur_pen == 0) {
                            int64_t l = lp_eval();
                            if (l >= 0) { cur_wl2 = l; lp_offer(l); }
                        }
                        cur_cost = cost_of(cur_wl2, cur_pen);
                    }
                } else {
                    bool plat = moves - last_imp > plateau && run_best > 0;
                    if (cycle_pos >= cur_cycle || plat) {
                        cycle_pos = 0;
                        cur_cycle = cycle_len;
                        last_imp = moves;
                        ++cycles;
                        if (plat && run_best < kick_best) kicks = 0;  // count only non-improving reheats
                        bool near = plat && cfg.s_near > 0 && run_best <= cfg.s_near * pen_norm && kicks < 3;
                        if (near) {  // near-feasible: shape descent, then a short reheat from the attempt's best state
                            S = snap;
                            shape_polish(out);
                            load_cur();
                            if (cur_pen < run_best) { run_best = cur_pen; snap = S; }
                            ++kicks;
                            kick_best = run_best;
                            Tc = std::min(cfg.reheat * probe(), 30.0 * lambda / pen_norm);
                            cur_cycle = kick_len;
                        } else if (plat && cfg.s_plat_restart > 0) {
                            kicks = 0;
                            kick_best = INT64_MAX;
                            bool legal = false;
                            if (cfg.s_gls > 0 && cfg.s_wl0 <= 0 && gls_plat_calls < cfg.s_gls_plat) {
                                // Phase F: legalize the attempt's best state with GLS before restarting.
                                const double bud = std::min(10.0, 0.25 * (deadline - t));
                                if (bud >= 1.0) {
                                    ++gls_plat_calls;
                                    S = snap;
                                    K.pack(P, S);
                                    legal = gls_try(to_rects(K, S), S.gidx, bud, out, "plateau");
                                }
                            }
                            if (legal) {
                                load_cur();
                                run_best = cur_pen;
                                snap = S;
                            } else {
                                State s = S;
                                if (next_start(s)) S = s;
                                ++restarts;
                                lambda = 1.0;
                                load_cur();
                                run_best = cur_pen;
                                snap = S;
                                Tc = cfg.s_restart_t * probe();
                            }
                        } else if (plat) {
                            Tc = std::min(cfg.reheat * probe(), 30.0 * lambda / pen_norm);  // <= 30 penalty units
                            cur_cycle = kick_len;
                        } else {
                            Tc = cfg.reheat * T0;
                        }
                    }
                    T = Tc * std::pow(1e-3, (double)cycle_pos / cur_cycle);
                }
                maybe_write(out, t);
            }
            ++cycle_pos;
            ++moves;
            Undo u;
            if (!propose_strong(u)) { ++rejected; continue; }
            K.pack(P, S);
            int64_t pen = K.ovfx + K.ovfy + K.viol, wl = K.wl2;
            bool lpok = false;
            if (lpmode && pen == 0) {
                int64_t l = lp_eval();
                if (l >= 0) { wl = l; lpok = true; }
            }
            double c = cost_of(wl, pen);
            double d = c - cur_cost;
            bool acc;
            if (phase2 && acc_rule == 1) {   // LAHC: better than the current state or than the state L moves ago
                double &h = lahc[lahc_i];
                acc = d <= 0 || c <= h;
                if (++lahc_i == lahc.size()) lahc_i = 0;
                h = acc ? c : cur_cost;   // history: the current cost after this move
            } else {
                acc = d <= 0 || rng.uni() < std::exp(-d / T);
            }
            if (phase2 && acc_rule == 2) {   // Modified Lam: steer T toward the target acceptance rate
                const double f = lam_f, mid = cfg.s_lam_mid;
                const double tgt = f < 0.15 ? mid + (cfg.s_lam_a0 - mid) * std::pow(560.0, -f / 0.15)
                                 : f < 0.65 ? mid : mid * std::pow(440.0, -(f - 0.65) / 0.35);
                if (d > 0) {   // uphill moves only (neutral moves would pin the rate above any low target)
                    lam_rate = 0.998 * lam_rate + 0.002 * (acc ? 1.0 : 0.0);
                    lam_T *= lam_rate > tgt ? 0.999 : 1.0 / 0.999;
                    lam_T = std::max(1e-3 * Tend, std::min(lam_T, T0));
                }
                T = lam_T;
            }
            if (phase2) {
                ++acc_n; acc_acc += acc;
                if (d > 0) { ++up_n; up_acc += acc; up_f[std::min(9, (int)(lam_f * 10))] += acc; up_fn[std::min(9, (int)(lam_f * 10))] += 1; }
            }
            if (acc) {
                if (lpok) lp_offer(wl);
                cur_cost = c; cur_wl2 = wl; cur_pen = pen; cur_ox = K.ovfx; cur_oy = K.ovfy;
                if (strong_inv) cur_forb = forb_count(P, S);
                if (pen > 0 || in_phase2) { cur_x = K.x; cur_y = K.y; }
                ++accepted;
                if (!phase2 && pen < run_best) {
                    run_best = pen;
                    last_imp = moves;
                    if ((cfg.s_near > 0 && pen <= cfg.s_near * pen_norm) || cfg.s_gls > 0) snap = S;
                }
                consider(out);
            } else {
                undo(u);
            }
            ++win_moves;
            if (cur_pen == 0) ++win_feas;
            if (win_moves >= WIN) {
                double fr = (double)win_feas / win_moves;
                double old = lambda;
                if (fr < 0.5) lambda = std::min(cfg.lam_max, lambda * 1.15);
                else if (fr > 0.9) lambda = std::max(0.05, lambda / 1.02);
                if (lambda != old) cur_cost = cost_of(cur_wl2, cur_pen);
                win_moves = win_feas = 0;
            }
        }
        std::fprintf(stderr, "STATS rej_boundary=%lld rej_preplaced_pair=%lld\n", (long long)g_rej_boundary, (long long)g_rej_pre);
        if (lpmode)
            std::fprintf(stderr, "STATS lpcost evals=%lld fail=%lld wins=%lld ms=%.0f phase_moves=%lld phase_s=%.2f\n",
                         (long long)lp_evals, (long long)lp_fail, (long long)lp_wins, lp_ms,
                         (long long)(moves - lp_m0), now_s(t0) - lp_t);
        if (stalls) std::fprintf(stderr, "STATS stalls=%lld\n", (long long)stalls);
        if (acc_n) {
            std::fprintf(stderr, "STATS phase2 rule=%d acc_rate=%.4f up_rate=%.4f lam_T=%.3g up_by_tenth=", acc_rule, (double)acc_acc / acc_n,
                         (double)up_acc / std::max<int64_t>(1, up_n), lam_T);
            for (int k = 0; k < 10; ++k) std::fprintf(stderr, "%.4f%s", up_f[k] / std::max(1.0, up_fn[k]), k < 9 ? "," : "\n");
        }
        if (gls_calls) std::fprintf(stderr, "STATS gls calls=%lld wins=%lld secs=%.2f\n", (long long)gls_calls, (long long)gls_wins, gls_secs);
        std::fprintf(stderr, "STATS cycles=%lld restarts=%lld T0=%.3g lambda=%.3g accepted=%lld rejected=%lld gate_checks=%lld gate_disagree=%lld time=%.2f moves_per_s=%.0f\n",
                     (long long)cycles, (long long)restarts, T0, lambda, (long long)accepted, (long long)rejected,
                     (long long)gate_checks, (long long)gate_disagree, now_s(t0), moves / std::max(1e-9, now_s(t0) - t_start));
    }
};

// Greedy soft-shape descent on a fixed SP: per shape group try a few width
// steps and keep any that lowers the normalized overflow + violation.
void shape_descent(const Problem &P, State &t, Packer &K2) {
    auto sc = [&]() {
        K2.pack(P, t);
        return (double)K2.ovfx / P.W + (double)K2.ovfy / P.H + (double)K2.viol / (P.W + P.H);
    };
    double cur = sc();
    for (int sweep = 0; sweep < 20 && cur > 0; ++sweep) {
        bool imp = false;
        for (int g = 0; g < (int)t.gidx.size(); ++g) {
            for (int d : {-3, -1, 1, 3}) {
                int o = t.gidx[g], nw = o + d;
                if (nw < 0 || nw >= (int)P.gshape[g].widths.size()) continue;
                t.gidx[g] = nw;
                apply_group(P, t, g);
                double c2 = sc();
                if (c2 < cur - 1e-12) { cur = c2; imp = true; }
                else { t.gidx[g] = o; apply_group(P, t, g); }
            }
        }
        if (!imp) break;
    }
}

// Strong mode: K analytical seeds, packed and ranked by penalty-weighted cost.
std::vector<State> make_seeds(const Problem &P, const Config &cfg, Rng &rng) {
    std::vector<State> seeds;
    std::vector<std::pair<double, int>> rank;
    std::vector<std::pair<int64_t, int64_t>> sc;
    Packer K;
    K.init(P.n);
    int ns = std::max(1, (int)cfg.s_seeds);
    for (int v = 0; v < ns; ++v) {
        GlobalPlacer G(P, rng);
        G.rand_axis = v > 0;
        G.run((int)cfg.s_gp_iters, cfg.s_gp_anchor);
        State s;
        std::string err;
        if (!sp_from_gp(P, s, G, cfg.s_grow > 0, cfg.s_geo_prio > 0, err)) { std::fprintf(stderr, "WARN seed %d: %s\n", v, err.c_str()); continue; }
        if (!invariant_ok(P, s)) { std::fprintf(stderr, "WARN seed %d breaks the invariant\n", v); continue; }
        if (cfg.s_shape_desc > 0) shape_descent(P, s, K);
        K.pack(P, s);
        sc.push_back({K.ovfx + K.ovfy + K.viol, K.wl2});
        seeds.push_back(s);
    }
    if (seeds.empty()) return seeds;
    int64_t wmin = INT64_MAX;
    for (auto &q : sc) wmin = std::min(wmin, q.second);
    for (size_t i = 0; i < seeds.size(); ++i)
        rank.push_back({(double)sc[i].second / std::max<int64_t>(1, wmin) + 10.0 * sc[i].first / (P.W + P.H), (int)i});
    std::sort(rank.begin(), rank.end());
    std::vector<State> out;
    for (auto &q : rank) out.push_back(seeds[q.second]);
    std::fprintf(stderr, "STATS best_seed pen=%lld wl2=%lld\n", (long long)sc[rank[0].second].first,
                 (long long)sc[rank[0].second].second);
    return out;
}

// ------------------------------------------------------------------ homotopy (Phase F, strong)
//
// Ported from solver/test_homotopy.cpp (the diagnostic harness). The relaxed
// instance at dimension scale s: fixed dims max(1, floor(true * sqrt s)), soft
// dims = the true shape at the group's width index scaled by sqrt s (rounded),
// preplaced blocks / outline / terminals untouched. Shape-list indices are the
// same at every s, so an SP + width indices carry over between scales.

// Relaxed Problem at scale s (homotopy phase 1). The forbidden-relation table
// stays the TRUE one, so the relaxed anneal never adds a relation that is
// impossible at s = 1.
Problem relax_problem(const Problem &P, double s) {
    Problem R = P;
    const double r = std::sqrt(s);
    for (int i = 0; i < P.n; ++i) {
        if (P.kind[i] != K_FIXED) continue;
        R.fw[i] = std::max(1, (int)std::floor(P.fw[i] * r));
        R.fh[i] = std::max(1, (int)std::floor(P.fh[i] * r));
    }
    for (auto &g : R.gshape) {
        g.A = std::max<int64_t>(1, std::llround(g.A * s));
        for (int &w : g.widths) w = std::max(1, (int)std::lround(w * r));
        for (int &h : g.heights) h = std::max(1, (int)std::lround(h * r));
    }
    return R;
}

void scale_dims(const Problem &P, State &S, double s) {
    S.w.assign(P.n, 0);
    S.h.assign(P.n, 0);
    const double r = std::sqrt(s);
    const bool one = s >= 1.0;
    for (int i = 0; i < P.n; ++i) {
        if (P.kind[i] == K_PRE) { S.w[i] = P.fw[i]; S.h[i] = P.fh[i]; }
        else if (P.kind[i] == K_FIXED) {
            S.w[i] = one ? P.fw[i] : std::max(1, (int)std::floor(P.fw[i] * r));
            S.h[i] = one ? P.fh[i] : std::max(1, (int)std::floor(P.fh[i] * r));
        }
    }
    for (int g = 0; g < (int)P.gmem.size(); ++g) {
        int w = P.gshape[g].widths[S.gidx[g]], h = P.gshape[g].heights[S.gidx[g]];
        if (!one) { w = std::max(1, (int)std::lround(w * r)); h = std::max(1, (int)std::lround(h * r)); }
        for (int b : P.gmem[g]) { S.w[b] = w; S.h[b] = h; }
    }
}

struct HEval {
    bool feas = false;
    int64_t ovf = 0, viol = 0, forb = 0, wl2 = 0;
    bool inv = false;
};

// Pack S at scale s (true problem P): feasible = inside the outline, pins
// satisfied, invariant, no forbidden relation; at s = 1 also the gate mirror.
HEval eval_at(const Problem &P, State &S, Packer &K, double s) {
    scale_dims(P, S, s);
    HEval e;
    e.inv = invariant_ok(P, S);
    e.forb = forb_count(P, S);
    K.pack(P, S);
    e.ovf = K.ovfx + K.ovfy;
    e.viol = K.viol;
    e.wl2 = K.wl2;
    e.feas = e.inv && e.forb == 0 && K.feasible();
    if (e.feas && s >= 1.0) e.feas = gate_check(P, to_rects(K, S)).ok;
    return e;
}

// Largest grid s in [lo, 1] (grid 0.001) with eval_at feasible; -1 if none
// (searches down to 0.80). Assumes rough monotonicity in s.
double h_bisect(const Problem &P, State &S, Packer &K, double lo) {
    auto ok = [&](int k) { return eval_at(P, S, K, k / 1000.0).feas; };
    int a = (int)std::lround(lo * 1000), b = 1000;
    if (ok(b)) return 1.0;
    if (!ok(a)) {
        int lo2 = 800;
        if (a <= lo2 || !ok(lo2)) return -1;
        b = a; a = lo2;
    }
    while (b - a > 1) {
        int m = (a + b) / 2;
        if (ok(m)) a = m; else b = m;
    }
    return a / 1000.0;
}

// Corrector SA at target st (floor sf: a fallback state must stay feasible
// there): cost wl2/wl_norm + 30 (ovf + viol)/(W + H) + forb at low T (1%
// increase accepted with p = 0.1); returns at the first state feasible at
// st, else keeps the lowest-cost visited state feasible at sf (else S).
bool h_corrector(const Problem &P, State &S, Packer &K, double st, double sf, double budget, double wl_norm, Rng &rng,
                 int64_t &moves) {
    const int n = P.n;
    const double lam = 30.0, lamf = 1.0;
    auto cost_of = [&](const HEval &e) {
        return (double)e.wl2 / wl_norm + lam * (double)(e.ovf + e.viol) / (P.W + P.H) + lamf * (double)e.forb;
    };
    const Clock::time_point c0 = Clock::now();
    State cur = S;
    HEval ec = eval_at(P, cur, K, st);
    if (ec.feas) { S = cur; return true; }
    double cc = cost_of(ec);
    const double T = 0.01 * cc / std::log(10.0);
    State best = S, nx;
    double bc = 1e300;
    const int ng = (int)P.gmem.size();
    for (int64_t m = 0;; ++m) {
        if ((m & 255) == 0 && std::chrono::duration<double>(Clock::now() - c0).count() > budget) break;
        ++moves;
        nx = cur;
        int mv = rng.below(ng > 0 ? 6 : 5);
        if (mv == 0) {
            int a = rng.below(n), b = rng.below(n);
            if (a == b) continue;
            int which = rng.below(3);
            if (which != 1) std::swap(nx.gp[nx.pp[a]], nx.gp[nx.pp[b]]);
            if (which != 0) std::swap(nx.gm[nx.pm[a]], nx.gm[nx.pm[b]]);
        } else if (mv == 1 || mv == 2) {
            int a = rng.below(n);
            auto ins = [&](std::vector<int> &q, std::vector<int> &pos) {
                int i = pos[a], j = rng.below(n);
                if (i == j) return;
                q.erase(q.begin() + i);
                q.insert(q.begin() + j, a);
            };
            if (mv == 2) { ins(nx.gp, nx.pp); ins(nx.gm, nx.pm); }
            else if (rng.below(2)) ins(nx.gp, nx.pp);
            else ins(nx.gm, nx.pm);
        } else if (mv == 3 || mv == 4) {
            int i = rng.below(n - 1);
            if (mv == 3) std::swap(nx.gp[i], nx.gp[i + 1]);
            else std::swap(nx.gm[i], nx.gm[i + 1]);
        } else {
            int g = rng.below(ng), sz = (int)P.gshape[g].widths.size();
            if (sz < 2) continue;
            int k = 1 + rng.below(std::max(1, std::min(5, sz / 10)));
            int v = nx.gidx[g] + (rng.below(2) ? k : -k);
            if (v < 0 || v >= sz) continue;
            nx.gidx[g] = v;
        }
        nx.sync_pos();
        if (!invariant_ok(P, nx)) continue;
        HEval e = eval_at(P, nx, K, st);
        if (e.feas) { S = nx; return true; }
        double c = cost_of(e);
        if (c <= cc || rng.uni() < std::exp(-(c - cc) / T)) {
            cur = nx; cc = c;
            if (c < bc) {
                State t = nx;
                if (eval_at(P, t, K, sf).feas) { best = nx; bc = c; }
            }
        }
    }
    if (bc < 1e300) S = best;
    return false;
}

// Homotopy flow (strong, s_shrink > 0; docs in the file header). sol is the
// true-instance solver: it owns the output, ttf and the final polish. The
// phases are fractions of the budget left at the call (Rem = deadline - now):
// standalone (s_hyb = 0) that is the whole run, in the hybrid what is left
// after the step-1 flow's first feasible state (+ s_hyb_at).
void run_homotopy(const Problem &P, const Config &cfg, Solver &sol, const std::string &out) {
    const double D = sol.deadline, s0 = cfg.s_shrink, base = now_s(sol.t0), Rem = std::max(0.0, D - base);
    sol.flush_output(out);
    const Problem Pr = relax_problem(P, s0);
    Config cr = cfg;
    cr.s_polish = 0;
    if (cfg.s_rng_split > 0) cr.seed = cfg.seed ^ 0xA5A5A5A5ULL;   // own RNG streams (seeds, SA, GLS) for the relaxed run
    const int64_t m0 = sol.moves;
    std::vector<std::unique_ptr<Solver>> cands;
    auto new_cand = [&](uint64_t seed, Clock::time_point t0) {
        Config ci = cr;
        ci.seed = seed;
        cands.emplace_back(new Solver(Pr, ci, t0));
        Solver &R = *cands.back();
        R.no_gate = true;
        R.tru = &sol;
        R.tru_base = m0;
        R.pool = make_seeds(Pr, ci, R.rng);
        if (!R.pool.empty()) { R.S = R.pool[0]; R.pool_next = 1; }
        return &R;
    };
    const int race_k = (int)cfg.s_race_k;
    int best = 0;
    if (race_k >= 2 && cfg.s_race_t > 0) {
        // Phase 1 race (Phase G, s_race_k): k relaxed runs on their own RNG
        // streams and clocks, each paused after s_race_t s of its own anneal;
        // the s_race_keep best (by relaxed incumbent WL) resume, one after the
        // other, on schedules planned so that phase 1 ends at base + P1.
        const int keep = std::max(1, std::min(race_k, (int)cfg.s_race_keep));
        // The race never takes more than 70% of the homotopy's budget and phase 1
        // never more than 85% (a late first-feasible state leaves little: ttf 45 s
        // on hidden2_ws0.0450 gave k = 5, s_race_t = 16 a 132 s wall before this clamp).
        const double rt = std::min(cfg.s_race_t, 0.7 * Rem / race_k);
        const double P1 = std::min(0.85 * Rem, std::max(std::max(0.0, cfg.s_h1) * Rem, race_k * rt + std::max(0.0, cfg.s_race_min) * Rem));
        const double cont = std::max(0.0, P1 - race_k * rt) / keep;
        for (int i = 0; i < race_k; ++i) {
            Solver *R = new_cand(cr.seed + (uint64_t)i * 0x9E3779B97F4A7C15ULL, Clock::now());
            R->deadline = rt + cont;
            R->pause_at = rt;
            R->pause_hard = true;
            if (!R->pool.empty()) R->run_strong(out, true);
            if (R->pol_pending) R->polish_inc();
            sol.flush_output(out);
        }
        std::vector<int> order(race_k);
        for (int i = 0; i < race_k; ++i) order[i] = i;
        auto key = [&](int i) { const Solver &R = *cands[i]; return R.have_feasible ? (double)R.inc_wl2 : 1e300 + i; };
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return key(a) < key(b); });
        std::fprintf(stderr, "HOMO race k=%d keep=%d t=%.1f cont=%.1f wl2=", race_k, keep, rt, cont);
        for (int i = 0; i < race_k; ++i)
            std::fprintf(stderr, "%lld%s", cands[i]->have_feasible ? (long long)cands[i]->inc_wl2 : -1LL, i + 1 < race_k ? "," : "");
        std::fprintf(stderr, " at=%.2f\n", now_s(sol.t0));
        for (int j = 0; j < keep; ++j) {
            Solver &R = *cands[order[j]];
            if (!R.paused) continue;   // ended before its pause point (nothing to resume)
            R.t0 += Clock::now() - R.pause_wall;
            R.pause_at = 1e18;
            R.run_strong(out, true);
            if (R.pol_pending) R.polish_inc();
            sol.flush_output(out);
        }
        best = order[0];
        for (int j = 1; j < keep; ++j)
            if (key(order[j]) < key(best)) best = order[j];
    } else {
        Solver *R = new_cand(cr.seed, sol.t0);
        R->deadline = base + std::max(0.0, cfg.s_h1) * Rem;
        // Phase 1: the strong flow on the relaxed instance.
        if (!R->pool.empty()) R->run_strong(out, true);
    }
    Solver &R = *cands[best];
    sol.moves = m0;
    for (auto &c : cands) sol.moves += c->moves;
    sol.flush_output(out);
    std::fprintf(stderr, "HOMO phase1 relaxed_feasible=%d relaxed_wl2=%lld true_feasible=%d t=%.2f\n", R.have_feasible ? 1 : 0,
                 (long long)R.inc_wl2, sol.have_feasible ? 1 : 0, now_s(sol.t0));
    // Phase 2: ramp s0 -> 1 from the relaxed packed incumbent (or no ramp).
    bool handed = false;
    if (R.have_feasible) {
        State S = R.inc_S;
        Packer K;
        K.init(P.n);
        const double t_end = std::min(D, base + (std::max(0.0, cfg.s_h1) + std::max(0.0, cfg.s_h2)) * Rem);
        Rng hr(cfg.seed * 0x2545F4914F6CDD1DULL + 99);
        const double wl_norm = std::max<double>(1.0, (double)R.inc_wl2);
        double s = cfg.s_h_mode > 0 ? -1 : h_bisect(P, S, K, s0);
        const double s_start = s;
        int steps = 0;
        int64_t hmoves = 0;
        while (cfg.s_h_mode <= 0 && s < 1.0 && now_s(sol.t0) < t_end) {
            double st = s < 0 ? s0 : std::min(1.0, s + 0.005), sf = s < 0 ? s0 : s;
            h_corrector(P, S, K, st, sf, std::min(1.0, t_end - now_s(sol.t0)), wl_norm, hr, hmoves);
            ++steps;
            double sn = h_bisect(P, S, K, s < 0 ? s0 : s);
            const bool stalled = !(sn > s);
            if (sn > s) s = sn;
            if (stalled && s >= 0.95) break;
        }
        sol.moves += hmoves;
        std::fprintf(stderr, "HOMO phase2 s_start=%.3f s=%.3f steps=%d moves=%lld t=%.2f\n", s_start, s, steps, (long long)hmoves,
                     now_s(sol.t0));
        if (s >= 1.0) {
            eval_at(P, S, K, 1.0);   // true dims
            sol.K.pack(P, S);
            if (sol.K.feasible()) {
                sol.S = S;
                sol.consider(out);
                handed = true;
            }
        }
        if (!handed) {
            // Absorber: GLS from the packing at s = 1 (mode 0) or from the
            // relaxed placement re-centred at true shapes (mode 1).
            double bud = std::min(std::min(cfg.s_h_gls, 0.15 * Rem), D - now_s(sol.t0) - 1.0);
            if (sol.have_feasible && cfg.s_h_gls_feas >= 0) bud = std::min(bud, cfg.s_h_gls_feas);
            std::vector<lab3::Rect> r0;
            std::vector<int> g0;
            if (cfg.s_h_mode > 0) { r0 = R.inc_rects; g0 = R.inc_S.gidx; }
            else { scale_dims(P, S, 1.0); K.pack(P, S); r0 = to_rects(K, S); g0 = S.gidx; }
            const State keep = sol.S;
            if (!r0.empty() && sol.gls_try(r0, g0, bud, out, "absorb")) handed = true;
            else sol.S = keep;
        }
    }
    sol.flush_output(out);
    std::fprintf(stderr, "HOMO handoff=%d true_feasible=%d t=%.2f\n", handed ? 1 : 0, sol.have_feasible ? 1 : 0, now_s(sol.t0));
    // Phase 3: post-feasible anneal from the best packed true incumbent, or
    // the normal true-instance flow when nothing legal came out.
    if (sol.have_feasible) {
        sol.S = sol.inc_S;
        sol.gls_handoff = true;
        sol.run_strong(out, false);
        return;
    }
    sol.pool = make_seeds(P, cfg, sol.rng);
    if (sol.pool.empty()) return;
    sol.S = sol.pool[0];
    sol.pool_next = 1;
    sol.run_strong(out, true);
}

// Restart portfolio (strong, s_rs; docs in the file header). s1 is the
// initialized first run. Leaves the process result in s1 (best_rects /
// best_wl2 / moves; ttf stays run 1's, the first feasible time of the process).
void run_portfolio(const Problem &P, const Config &cfg, Solver &s1, const std::string &out) {
    const double D = s1.deadline, split = std::max(0.0, std::min(1.0, cfg.s_rs_frac)) * D;
    Config c2 = cfg;
    c2.seed = cfg.seed * 0x9E3779B97F4A7C15ULL + 0x632BE59BD9B4E019ULL;
    std::unique_ptr<Solver> s2;
    auto start2 = [&](Clock::time_point t0b) {
        s2.reset(new Solver(P, c2, t0b));
        s2->pool = make_seeds(P, c2, s2->rng);
        if (s2->pool.empty()) return false;
        s2->S = s2->pool[0];
        s2->pool_next = 1;
        s2->rival = &s1;
        return true;
    };
    auto flush = [&](Solver &s) {
        if (s.pol_pending) s.polish_inc();
        if (s.have_feasible && s.may_write()) write_output(P, s.best_rects, out);
    };
    auto merge = [&](Solver &w) {
        s1.moves += s2->moves;
        if (&w != &s1) { s1.best_rects = w.best_rects; s1.best_wl2 = w.best_wl2; }
    };
    if (cfg.s_rs < 1.5) {   // 1: run 1 to the split point, then run 2 to the deadline; best of both
        s1.deadline = split;
        s1.pause_at = split;
        s1.ext_deadline = D;
        s1.run_strong(out, true);
        if (!s1.paused) return;
        s1.on_incumbent(true);
        s1.final_done = true;
        flush(s1);
        std::fprintf(stderr, "PORTFOLIO run1 wl2=%lld t=%.2f\n", (long long)s1.best_wl2, now_s(s1.t0));
        if (!start2(s1.t0)) return;
        s2->deadline = D;
        s2->run_strong(out, true);
        if (s2->have_feasible) s2->on_incumbent(true);
        bool w2 = s2->have_feasible && s2->best_wl2 < s1.best_wl2;
        std::fprintf(stderr, "PORTFOLIO run2 feasible=%d wl2=%lld winner=%d\n", s2->have_feasible ? 1 : 0,
                     (long long)s2->best_wl2, w2 ? 2 : 1);
        merge(w2 ? *s2 : s1);
        return;
    }
    // 2: race run 1 and run 2 for `split` s each on schedules planned to D - split,
    // then resume the better one (its clock shifted by the time it was paused).
    const double Dv = D - split;
    s1.deadline = Dv;
    s1.pause_at = split;
    s1.ext_deadline = D;
    s1.run_strong(out, true);
    if (!s1.paused) return;
    flush(s1);
    if (!start2(Clock::now())) { s1.t0 += Clock::now() - s1.pause_wall; s1.pause_at = 1e18; s1.run_strong(out, true); return; }
    s2->deadline = Dv;
    s2->pause_at = split;
    s2->pause_hard = true;
    s2->run_strong(out, true);
    if (s2->pol_pending) s2->polish_inc();
    bool w2 = s2->have_feasible && s2->best_wl2 < s1.best_wl2 && s2->rs.valid;
    std::fprintf(stderr, "PORTFOLIO race run1 wl2=%lld run2 feasible=%d wl2=%lld winner=%d\n", (long long)s1.best_wl2,
                 s2->have_feasible ? 1 : 0, (long long)s2->best_wl2, w2 ? 2 : 1);
    Solver &w = w2 ? *s2 : s1;
    if (w2) flush(w);
    w.t0 += Clock::now() - w.pause_wall;
    w.pause_at = 1e18;
    w.run_strong(out, true);
    if (w.have_feasible) w.on_incumbent(true);
    w.final_done = true;
    merge(w);
    s1.final_done = true;
}

int usage(const char *argv0) {
    std::fprintf(stderr, "Usage: %s <input> <output> [--time SEC=24] [--seed N=1] [--mode weak|strong] "
                         "[--init-from PLACEMENT] [--pack-only] [--set key=value]...\n", argv0);
    return 2;
}

}  // namespace

int main(int argc, char **argv) {
    Clock::time_point t0 = Clock::now();
    if (argc < 3) return usage(argv[0]);
    std::string in = argv[1], out = argv[2];
    Config cfg;
    for (int i = 3; i < argc; ++i) {
        std::string a = argv[i];
        auto val = [&]() -> const char * { return i + 1 < argc ? argv[++i] : nullptr; };
        if (a == "--time") { const char *v = val(); if (!v) return usage(argv[0]); cfg.time = std::atof(v); }
        else if (a == "--seed") { const char *v = val(); if (!v) return usage(argv[0]); cfg.seed = std::strtoull(v, nullptr, 10); }
        else if (a == "--mode") {
            const char *v = val();
            if (!v) return usage(argv[0]);
            std::string m = v;
            if (m == "weak") cfg.strong = false;
            else if (m == "strong") cfg.strong = true;
            else return usage(argv[0]);
        }
        else if (a == "--init-from") { const char *v = val(); if (!v) return usage(argv[0]); cfg.init_from = v; }
        else if (a == "--pack-only") cfg.pack_only = true;
        else if (a == "--set") { const char *v = val(); if (!v || !cfg.set(v)) return usage(argv[0]); }
        else return usage(argv[0]);
    }

    if (cfg.strong && cfg.s_wl0 > 0) { cfg.s_seeds = 0; cfg.s_polish = 0; }   // sanity run: no WL-driven step at all
    Problem P;
    std::string err;
    if (!build_problem(in, P, err)) { std::fprintf(stderr, "Input error: %s\n", err.c_str()); return 2; }

    if (cfg.strong) build_forbidden(P);
    Solver sol(P, cfg, t0);
    State &S = sol.S;
    const bool homo = cfg.strong && cfg.s_shrink > 0 && cfg.s_shrink < 1 && cfg.init_from.empty() && cfg.s_seeds > 0 &&
                      cfg.s_wl0 <= 0 && cfg.s_rs <= 0 && !cfg.pack_only && cfg.s_hyb <= 0;
    const bool hyb = cfg.strong && cfg.s_hyb > 0 && cfg.s_shrink > 0 && cfg.s_shrink < 1 && cfg.init_from.empty() &&
                     cfg.s_seeds > 0 && cfg.s_wl0 <= 0 && cfg.s_rs <= 0 && !cfg.pack_only;
    if (!cfg.init_from.empty()) {
        if (!initial_from_placement(P, S, cfg.init_from, sol.rng, err)) {
            std::fprintf(stderr, "init-from failed: %s\n", err.c_str());
            std::fprintf(stderr, "RESULT feasible=0 ttf=-1 wl=0.0 boundary=0/%zu mib=0/%zu moves=0\n",
                         P.inst.boundaries.size(), P.inst.mibs.size());
            return 1;
        }
    } else if (homo) {
        // Homotopy: seeds are built per phase in run_homotopy(). S gets a
        // valid (random, near-square) state so that nothing ever packs or
        // reads an empty SP (separate Rng: sol.rng's stream is unchanged).
        S.gidx.assign(P.gmem.size(), 0);
        for (int g = 0; g < (int)P.gmem.size(); ++g) {
            const auto &ws = P.gshape[g].widths;
            int target = (int)std::lround(std::sqrt((double)P.gshape[g].A));
            S.gidx[g] = std::min((int)(std::lower_bound(ws.begin(), ws.end(), target) - ws.begin()), (int)ws.size() - 1);
        }
        init_dims(P, S);
        Rng hr(cfg.seed ^ 0x9E3779B97F4A7C15ULL);
        if (!random_initial_sp(P, S, hr, err)) { std::fprintf(stderr, "initial SP failed: %s\n", err.c_str()); return 1; }
    } else if (cfg.strong && cfg.s_seeds > 0) {
        double tg = now_s(t0);
        sol.pool = make_seeds(P, cfg, sol.rng);
        std::fprintf(stderr, "STATS seeds=%zu seed_time=%.3f\n", sol.pool.size(), now_s(t0) - tg);
        if (sol.pool.empty()) { std::fprintf(stderr, "no analytical seed\n"); return 1; }
        size_t first = std::min(sol.pool.size() - 1, (size_t)std::max(0.0, cfg.s_first));
        S = sol.pool[first];
        sol.pool_next = first + 1;
    } else {
        S.gidx.assign(P.gmem.size(), 0);
        for (int g = 0; g < (int)P.gmem.size(); ++g) {  // start near square
            const auto &ws = P.gshape[g].widths;
            int target = (int)std::lround(std::sqrt((double)P.gshape[g].A));
            S.gidx[g] = (int)(std::lower_bound(ws.begin(), ws.end(), target) - ws.begin());
            if (S.gidx[g] >= (int)ws.size()) S.gidx[g] = (int)ws.size() - 1;
        }
        init_dims(P, S);
        if (!random_initial_sp(P, S, sol.rng, err)) {
            std::fprintf(stderr, "initial SP failed: %s\n", err.c_str());
            return 1;
        }
    }
    bool inv = homo || invariant_ok(P, S);
    if (cfg.strong && !homo)
        std::fprintf(stderr, "INIT forbidden_pairs=%zu forbidden_now=%lld\n", P.forb.size(), (long long)forb_count(P, S));
    std::fprintf(stderr, "INIT invariant=%d n=%d groups=%d mode=%s\n", inv ? 1 : 0, P.n, (int)P.gmem.size(),
                 cfg.strong ? "strong" : "weak");

    if (cfg.pack_only) {
        sol.K.pack(P, S);
        sol.consider(out);
    } else if (cfg.strong) {
        if (homo) run_homotopy(P, cfg, sol, out);
        else if (hyb) {
            // Hybrid (Phase F): the step-1 flow (GLS-first, SA) secures the
            // incumbent; s_hyb_at s after its first feasible state the rest of
            // the budget goes to the homotopy. The homotopy only offers
            // candidates through consider() (feasible + gate mirror + lower
            // WL), and the switch happens only once feasible, so the hybrid is
            // never less feasible than the step-1 flow on the same run prefix.
            sol.hyb_at = std::max(0.0, cfg.s_hyb_at);
            sol.run_strong(out, true);
            if (sol.hyb_go) run_homotopy(P, cfg, sol, out);
        }
        else if (cfg.s_rs > 0 && cfg.init_from.empty() && cfg.s_seeds > 0 && cfg.s_wl0 <= 0) run_portfolio(P, cfg, sol, out);
        else sol.run_strong(out, cfg.init_from.empty() && cfg.s_seeds > 0);
    } else {
        sol.run(out);
    }
    if (cfg.strong && sol.have_feasible && !sol.final_done) sol.on_incumbent(true);

    std::vector<lab3::Rect> r = sol.have_feasible ? sol.best_rects
                              : (!sol.fallback_rects.empty() ? sol.fallback_rects : to_rects(sol.K, S));
    if (!write_output(P, r, out)) { std::fprintf(stderr, "cannot write %s\n", out.c_str()); return 1; }
    GateInfo gi = gate_check(P, r);
    if (!sol.have_feasible && !gi.ok)
        std::fprintf(stderr, "INFEASIBLE best penalty=%lld (%s) ovf=%lld+%lld viol=%lld\n", (long long)sol.best_pen,
                     gi.why.c_str(), (long long)sol.K.ovfx, (long long)sol.K.ovfy, (long long)sol.K.viol);
    std::fprintf(stderr, "RESULT feasible=%d ttf=%.3f wl=%s boundary=%d/%zu mib=%d/%zu moves=%lld\n",
                 gi.ok ? 1 : 0, sol.have_feasible ? sol.ttf : -1.0,
                 lab3::format_half(lab3::wirelength2(P.inst, r)).c_str(), gi.ok ? gi.bsat : 0,
                 P.inst.boundaries.size(), gi.ok ? gi.msat : 0, P.inst.mibs.size(), (long long)sol.moves);
    return 0;
}
