// test_oracle.cpp -- InsertionOracle (oracle.h) checks and go/no-go stats.
//
//   test_oracle check <case.in> <placement.out> [n_sp=3] [seed=1]
//       Correctness against brute force: SPs extracted from the placement with
//       n_sp random tie-break seeds; every block removed in turn; every cell
//       (i, j) inserted into a real State, packed by Packer (+ flush).
//         - legal verdict: oracle (mask && fit) == invariant_ok && feasible,
//           and gate_check of the packed shape-legal placement passes;
//         - invariant mask == invariant_ok (non-preplaced b);
//         - windows: (wx, wy) feasible, (wx+1, wy) and (wx, wy+1) infeasible,
//           b packed at (xlo, ylo) (W-w / H-h when b is RIGHT / TOP);
//         - illegal soft cells: on a 2% sample, every legal shape is packed
//           and must be infeasible.
//   test_oracle gng <case.in> <placement.out> [seed=1] [topk=8]
//       Go/no-go statistics on a finished placement (see report at the end).
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "gp.h"
#include "oracle.h"

using namespace fp;

static int g_shapes = 0;   // env ORACLE_SHAPES: extra soft shapes per evaluated cell

static double now_us() {
    return std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static bool load(const std::string &in, const std::string &out, Problem &P, State &S, uint64_t seed) {
    std::string err;
    if (!build_problem(in, P, err)) { std::fprintf(stderr, "problem: %s\n", err.c_str()); return false; }
    build_forbidden(P);
    Rng rng(seed);
    if (!initial_from_placement(P, S, out, rng, err)) { std::fprintf(stderr, "placement: %s\n", err.c_str()); return false; }
    return true;
}

static bool state_legal(const Problem &P, const State &S, Packer &K) {
    K.pack(P, S);
    return K.feasible() && invariant_ok(P, S);
}

// ------------------------------------------------------------------ check

struct Counts {
    long long cells = 0, legal = 0, verdict_bad = 0, inv_bad = 0, gate_bad = 0;
    long long probe_lo_bad = 0, probe_hix_bad = 0, probe_hiy_bad = 0, pos_bad = 0, unverif = 0;
    long long shape_sample = 0, shape_bad = 0, geo_only = 0, orig_bad = 0, reduced_bad = 0;
};

static int run_check(const std::string &in, const std::string &out, int nsp, uint64_t seed0) {
    Counts C;
    double t_oracle = 0;
    long long n_oracle = 0;
    for (int s = 0; s < nsp; ++s) {
        Problem P;
        State S;
        if (!load(in, out, P, S, seed0 + s)) return 2;
        Packer K;
        K.init(P.n);
        if (!state_legal(P, S, K) || !gate_check(P, to_rects(K, S)).ok) {
            std::printf("SP %d: base SP not legal, skipped\n", s);
            continue;
        }
        Rng rng(1000 + seed0 + s);
        InsertionOracle O;
        State T = S;
        for (int b = 0; b < P.n; ++b) {
            double t0 = now_us();
            O.build(P, S, b);
            t_oracle += now_us() - t0;
            ++n_oracle;
            if (!O.reduced_ok) { ++C.reduced_bad; continue; }
            const int c0 = O.idx(O.ib, O.jb);
            if (!(O.inv[c0] && O.fits(P, c0, S.w[b], S.h[b]))) ++C.orig_bad;
            const bool pin = P.kind[b] == K_PRE;
            for (int i = 0; i <= O.m; ++i) {
                for (int j = 0; j <= O.m; ++j) {
                    const int c = O.idx(i, j);
                    ++C.cells;
                    O.insert(i, j, T.gp, T.gm);
                    T.sync_pos();
                    if (i == O.ib && j == O.jb && (T.gp != S.gp || T.gm != S.gm)) ++C.verdict_bad;
                    const bool invb = invariant_ok(P, T);
                    if (!pin && invb != (bool)O.inv[c]) ++C.inv_bad;
                    // Verdict with a real legal shape.
                    const bool ol = O.legal(P, S, c);
                    C.legal += ol;
                    int w = S.w[b], h = S.h[b];
                    if (ol && O.free_shape) {
                        int k = O.pick_shape(P, S, c);
                        w = P.gshape[P.group[b]].widths[k];
                        h = P.gshape[P.group[b]].heights[k];
                    }
                    T.w[b] = w; T.h[b] = h;
                    K.pack(P, T);
                    const bool bl = K.feasible() && invb;
                    if (bl != ol) ++C.verdict_bad;
                    if (bl) {
                        GateInfo gi = gate_check(P, to_rects(K, T));
                        if (!gi.ok) {
                            ++C.gate_bad;
                            if (C.gate_bad <= 5) std::printf("  gate_check fail b=%d cell=(%d,%d) %s\n", b, i, j, gi.why.c_str());
                        }
                    }
                    if (invb) {
                        // Packer is conservative: an SP pair c left-of pin p
                        // must have x_c + w_c <= px even when the two are
                        // y-disjoint, so geometry can pass where Packer fails.
                        bool go = mcf::gate_ok(P.inst, to_rects(K, T), to_rects(K, T));
                        if (bl && !go) ++C.gate_bad;
                        if (!bl && go) ++C.geo_only;
                    }
                    // Illegal soft (geometry): sample, pack every shape.
                    if (O.free_shape && !O.geo_ok(P, S, c) && rng.below(50) == 0) {
                        ++C.shape_sample;
                        const ShapeRange &g = P.gshape[P.group[b]];
                        for (size_t k = 0; k < g.widths.size(); ++k) {
                            T.w[b] = g.widths[k]; T.h[b] = g.heights[k];
                            K.pack(P, T);
                            if (K.feasible()) { ++C.shape_bad; break; }
                        }
                    }
                    // Window probes (movable b).
                    if (!pin) {
                        const int wx = O.xhi[c] - O.xlo[c], wy = O.yhi[c] - O.ylo[c];
                        auto feas = [&](int pw, int ph) {
                            T.w[b] = pw; T.h[b] = ph;
                            K.pack(P, T);
                            return K.feasible();
                        };
                        if (wx >= 1 && wy >= 1) {
                            if (!feas(wx, wy)) ++C.probe_lo_bad;
                            else {
                                int ex = (P.bmask[b] & lab3::EDGE_RIGHT) ? P.W - wx : O.xlo[c];
                                int ey = (P.bmask[b] & lab3::EDGE_TOP) ? P.H - wy : O.ylo[c];
                                if (K.x[b] != ex || K.y[b] != ey) ++C.pos_bad;
                            }
                            if (feas(wx + 1, wy)) ++C.probe_hix_bad;
                            if (feas(wx, wy + 1)) ++C.probe_hiy_bad;
                        } else {
                            // A window is 0 (windows are never negative: every
                            // left-set block precedes every right-set block).
                            // Zero-extent probes isolate the other axis; Packer
                            // accepts w or h = 0 (gate_check is not used here).
                            ++C.unverif;
                            if (!feas(wx, wy)) ++C.probe_lo_bad;
                            if (feas(wx + 1, wy)) ++C.probe_hix_bad;
                            if (feas(wx, wy + 1)) ++C.probe_hiy_bad;
                            if (wx < 0 || wy < 0) ++C.probe_lo_bad;
                        }
                    }
                    T.w[b] = S.w[b]; T.h[b] = S.h[b];
                }
            }
        }
    }
    long long bad = C.verdict_bad + C.inv_bad + C.gate_bad + C.probe_lo_bad + C.probe_hix_bad + C.probe_hiy_bad +
                    C.pos_bad + C.shape_bad + C.orig_bad + C.reduced_bad;
    std::printf("CHECK %s cells=%lld legal=%lld verdict_bad=%lld inv_bad=%lld gate_bad=%lld probe_lo_bad=%lld "
                "probe_hix_bad=%lld probe_hiy_bad=%lld pos_bad=%lld shape_sample=%lld shape_bad=%lld orig_bad=%lld "
                "reduced_bad=%lld zero_window_cells=%lld packer_reject_geo_ok=%lld oracle_us=%.1f %s\n",
                in.c_str(), C.cells, C.legal, C.verdict_bad, C.inv_bad, C.gate_bad, C.probe_lo_bad, C.probe_hix_bad,
                C.probe_hiy_bad, C.pos_bad, C.shape_sample, C.shape_bad, C.orig_bad, C.reduced_bad, C.unverif, C.geo_only,
                n_oracle ? t_oracle / n_oracle : 0.0, bad ? "FAIL" : "OK");
    return bad ? 1 : 0;
}

// ------------------------------------------------------------------ go/no-go

// 2 * (weighted L1 WL of the nets touching b) with b's doubled centre (cx2, cy2)
// and every other block at `r`.
static int64_t wl2_of(const Problem &P, const std::vector<std::vector<std::pair<int, int>>> &adjb,
                      const std::vector<std::vector<std::array<int, 3>>> &adjt, const std::vector<lab3::Rect> &r,
                      int b, int64_t cx2, int64_t cy2) {
    int64_t s = 0;
    for (auto &e : adjb[b]) {
        const lab3::Rect &q = r[e.first];
        s += (int64_t)e.second * (std::llabs(cx2 - (2 * q.x + q.w)) + std::llabs(cy2 - (2 * q.y + q.h)));
    }
    for (auto &t : adjt[b]) s += (int64_t)t[2] * (std::llabs(cx2 - t[0]) + std::llabs(cy2 - t[1]));
    (void)P;
    return s;
}

static int64_t wmedian(std::vector<std::pair<int64_t, int>> v) {
    std::sort(v.begin(), v.end());
    int64_t tot = 0, acc = 0;
    for (auto &p : v) tot += p.second;
    for (auto &p : v) { acc += p.second; if (2 * acc >= tot) return p.first; }
    return v.empty() ? 0 : v.back().first;
}

static int run_gng(const std::string &in, const std::string &out, uint64_t seed, int topk) {
    Problem P;
    State S;
    if (!load(in, out, P, S, seed)) return 2;
    const int n = P.n;
    Packer K;
    K.init(n);
    if (!state_legal(P, S, K)) { std::printf("base SP not legal\n"); return 2; }
    std::vector<lab3::Rect> filer;
    std::string err;
    read_placement(P, out, filer, err);
    const int64_t file_wl2 = lab3::wirelength2(P.inst, filer);
    std::vector<lab3::Rect> base_asap = to_rects(K, S);
    const int64_t base_asap2 = K.wl2;
    mcf::ConstraintGraphs sg;
    sp_graphs(S.gp, S.pm, sg);
    mcf::PolishInfo pi;
    std::vector<lab3::Rect> base = mcf::polish_topology(P.inst, base_asap, sg, &pi);
    const int64_t base_pol2 = pi.wl2_out;
    std::printf("CASE %s n=%d file_wl=%.1f base_asap=%.1f base_pol=%.1f polok=%d\n", in.c_str(), n, file_wl2 / 2.0,
                base_asap2 / 2.0, base_pol2 / 2.0, (int)pi.ok);

    std::vector<std::vector<std::pair<int, int>>> adjb(n);
    std::vector<std::vector<std::array<int, 3>>> adjt(n);
    for (const BBConn &c : P.bb) { adjb[c.a].push_back({c.b, c.w}); adjb[c.b].push_back({c.a, c.w}); }
    for (const BTConn &c : P.bt) adjt[c.a].push_back({c.tx2, c.ty2, c.w});

    InsertionOracle O;
    State T = S;
    std::vector<double> t_or, t_eval, gains_pol, gains_asap;
    int movable = 0, with_legal = 0, with_geo = 0, impr_pol = 0, impr_asap = 0, impr_proxy = 0;
    long long sum_legal = 0, sum_geo = 0;
    for (int b = 0; b < n; ++b) {
        double t0 = now_us();
        O.build(P, S, b);
        t_or.push_back(now_us() - t0);
        if (P.kind[b] == K_PRE) continue;
        ++movable;
        double t1 = now_us();
        // Net-weighted median target (doubled centre).
        std::vector<std::pair<int64_t, int>> vx, vy;
        for (auto &e : adjb[b]) {
            const lab3::Rect &q = base[e.first];
            vx.push_back({2 * q.x + q.w, e.second});
            vy.push_back({2 * q.y + q.h, e.second});
        }
        for (auto &t : adjt[b]) { vx.push_back({t[0], t[2]}); vy.push_back({t[1], t[2]}); }
        const int64_t tx2 = wmedian(vx), ty2 = wmedian(vy);
        const lab3::Rect &rb = base[b];
        const int64_t old_part = wl2_of(P, adjb, adjt, base, b, 2 * rb.x + rb.w, 2 * rb.y + rb.h);
        struct Cand { int64_t gain; int c, w, h; };
        std::vector<Cand> cand;
        int nl = 0, ng = 0;
        for (int i = 0; i <= O.m; ++i) {
            for (int j = 0; j <= O.m; ++j) {
                if (i == O.ib && j == O.jb) continue;
                const int c = O.idx(i, j);
                if (!O.geo_ok(P, S, c)) continue;
                ++ng;
                if (!O.inv[c]) continue;
                ++nl;
                int w = S.w[b], h = S.h[b];
                if (O.free_shape) {
                    int k = O.pick_shape(P, S, c);
                    w = P.gshape[P.group[b]].widths[k];
                    h = P.gshape[P.group[b]].heights[k];
                }
                auto place = [](int64_t t2, int sz, int lo, int hi) {
                    int64_t v = (t2 - sz) / 2;
                    if (v > hi - sz) v = hi - sz;
                    if (v < lo) v = lo;
                    return v;
                };
                int64_t x = place(tx2, w, O.xlo[c], O.xhi[c]), y = place(ty2, h, O.ylo[c], O.yhi[c]);
                if (P.bmask[b] & lab3::EDGE_LEFT) x = 0;
                if (P.bmask[b] & lab3::EDGE_RIGHT) x = P.W - w;
                if (P.bmask[b] & lab3::EDGE_BOTTOM) y = 0;
                if (P.bmask[b] & lab3::EDGE_TOP) y = P.H - h;
                int64_t np = wl2_of(P, adjb, adjt, base, b, 2 * x + w, 2 * y + h);
                cand.push_back({old_part - np, c, w, h});
            }
        }
        sum_legal += nl; sum_geo += ng;
        with_legal += nl > 0; with_geo += ng > 0;
        std::sort(cand.begin(), cand.end(), [](const Cand &a, const Cand &c) { return a.gain > c.gain; });
        int64_t best_proxy = cand.empty() ? 0 : cand[0].gain;
        int64_t best_pol = INT64_MIN, best_asap = INT64_MIN;
        for (int t = 0; t < (int)cand.size() && t < topk; ++t) {
            const Cand &cd = cand[t];
            O.insert(cd.c / (O.m + 1), cd.c % (O.m + 1), T.gp, T.gm);
            T.sync_pos();
            // Shapes tried: the picked one; with ORACLE_SHAPES=k also up to k
            // more fitting widths spread over the fitting range (soft, non-MIB).
            std::vector<std::pair<int, int>> shp{{cd.w, cd.h}};
            if (g_shapes > 0 && O.free_shape) {
                const ShapeRange &g = P.gshape[P.group[b]];
                std::vector<int> fit;
                for (int k = 0; k < (int)g.widths.size(); ++k)
                    if (g.widths[k] <= O.xhi[cd.c] - O.xlo[cd.c] && g.heights[k] <= O.yhi[cd.c] - O.ylo[cd.c]) fit.push_back(k);
                for (int q = 0; q < g_shapes && !fit.empty(); ++q) {
                    int k = fit[(size_t)((fit.size() - 1) * (double)q / std::max(1, g_shapes - 1))];
                    if (g.widths[k] != cd.w) shp.push_back({g.widths[k], g.heights[k]});
                }
            }
            for (auto &sh : shp) {
                T.w[b] = sh.first; T.h[b] = sh.second;
                K.pack(P, T);
                std::vector<lab3::Rect> r = to_rects(K, T);
                if (!(K.feasible() && invariant_ok(P, T) && gate_check(P, r).ok)) {
                    std::printf("ERROR b=%d cell legal by oracle but pack illegal\n", b);
                    continue;
                }
                best_asap = std::max(best_asap, base_asap2 - K.wl2);
                mcf::ConstraintGraphs g2;
                sp_graphs(T.gp, T.pm, g2);
                mcf::PolishInfo p2;
                mcf::polish_topology(P.inst, r, g2, &p2);
                best_pol = std::max(best_pol, base_pol2 - p2.wl2_out);
            }
        }
        T.w[b] = S.w[b]; T.h[b] = S.h[b];
        t_eval.push_back(now_us() - t1);
        if (best_proxy > 0) ++impr_proxy;
        if (best_pol > 0) { ++impr_pol; gains_pol.push_back(100.0 * best_pol / base_pol2); }
        if (best_asap > 0) { ++impr_asap; gains_asap.push_back(100.0 * best_asap / base_asap2); }
        std::printf("BLK %d kind=%d mib=%d bnd=%d legal_nonorig=%d geo_nonorig=%d proxy=%.1f asap_gain=%.1f pol_gain=%.1f "
                    "oracle_us=%.1f eval_ms=%.1f\n",
                    b, P.kind[b], (int)(P.kind[b] == K_SOFT && !O.free_shape), P.bmask[b], nl, ng, best_proxy / 2.0,
                    best_asap == INT64_MIN ? 0.0 : best_asap / 2.0, best_pol == INT64_MIN ? 0.0 : best_pol / 2.0,
                    t_or.back(), (now_us() - t1) / 1000.0);
    }
    auto med = [](std::vector<double> v) {
        if (v.empty()) return 0.0;
        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
    };
    double mean_or = 0;
    for (double v : t_or) mean_or += v;
    mean_or /= std::max<size_t>(1, t_or.size());
    std::printf("SUMMARY %s topk=%d shapes=%d n=%d movable=%d with_legal_nonorig=%.3f with_geo_nonorig=%.3f mean_legal_cells=%.1f "
                "mean_geo_cells=%.1f proxy_improving=%.3f asap_improving=%.3f pol_improving=%.3f "
                "median_pol_gain_pct=%.3f median_asap_gain_pct=%.3f oracle_us_median=%.1f oracle_us_mean=%.1f "
                "eval_ms_median=%.1f\n",
                in.c_str(), topk, g_shapes, n, movable, (double)with_legal / movable, (double)with_geo / movable,
                (double)sum_legal / movable, (double)sum_geo / movable, (double)impr_proxy / movable,
                (double)impr_asap / movable, (double)impr_pol / movable, med(gains_pol), med(gains_asap), med(t_or),
                mean_or, med(t_eval) / 1000.0);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: %s check|gng <case.in> <placement.out> [...]\n", argv[0]);
        return 2;
    }
    if (const char *e = std::getenv("ORACLE_SHAPES")) g_shapes = std::atoi(e);
    std::string mode = argv[1];
    if (mode == "check")
        return run_check(argv[2], argv[3], argc > 4 ? std::atoi(argv[4]) : 3, argc > 5 ? std::atoll(argv[5]) : 1);
    if (mode == "gng")
        return run_gng(argv[2], argv[3], argc > 4 ? std::atoll(argv[4]) : 1, argc > 5 ? std::atoi(argv[5]) : 8);
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}
