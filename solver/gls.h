// gls.h -- coordinate-space overlap removal with guided local search (GLS).
//
// Fallback legalizer for the strong solver: when SP simulated annealing
// stalls before its first feasible packing, take the best (overflowing)
// packed placement, clamp it into the outline and minimize the total
// pairwise overlap directly on block coordinates. Once the overlap is zero
// the placement is legal; extract_sp (extract.h) turns it back into an SP
// that packs legal, and the normal SA continues from there.
//
// Invariants kept by every move (so "zero overlap" == "legal placement"):
//   * every block inside the outline,
//   * preplaced blocks at their fixed position (never moved),
//   * boundary blocks on their required edge(s),
//   * soft shapes from the legal width list of their shape group; an MIB
//     group is one shape group, so all members always share one shape.
//
// Cost of a pair (u = sqrt(W*H/n), the mean block side):
//   ov(i,j) = dx*dy/u^2 + 0.1*min(dx,dy)/u   (dx, dy > 0 the overlap extents)
// SA minimizes the GLS-weighted sum  sum_{i<j} wgt(i,j) * ov(i,j).  Every
// gls_every steps while infeasible, wgt += gls_inc on every currently
// overlapping pair: a pair that keeps overlapping becomes more expensive,
// which pushes the search out of the local minimum it sits in.
//
// Moves on a block i (90% picked among blocks that overlap something):
//   shape  (soft)  new width index of i's group within +-2.5% (10%: +-15%)
//                  of the list, each member anchored at its centre, its
//                  lower-left or its upper-right corner;
//   swap           exchange the centres of i and another movable block;
//   snap           x or y of i to an edge of a block j (abut left/right,
//                  align low/high edges), to the outline, or a small shift;
//                  j is one of i's overlap partners 65% of the time;
//   push           move i along one axis to clear a partner (or a small
//                  step) and push every block it then overlaps to abut it,
//                  recursively; rejected if a pushed block would leave the
//                  outline, leave its boundary edge or move a preplaced one.
// Changed blocks are clamped into the outline / onto their boundary edges.
// Temperature: T0 * (T1/T0)^phase with phase cycling in [0,1) every
// t_period seconds (wall clock), as a sequence of short anneals.
//
// WL bias (wlw): while infeasible the SA cost also carries
//   wlw(t) * dWL / (u * #connections),  wlw(t) = wlw * 0.5^(t / wlw_half),
// (optionally fading out; GLS weights grow meanwhile, so the overlap term
// dominates over time anyway). Measured on dev_tight (test_gls, seed start,
// 16 cases x 5 seeds, legal within 30 s / WL of the legal result over the
// construction WL): wlw 0: 99% / 3.2x; wlw 0.1: 99% / 2.1x (default);
// wlw 0.3: 91% / 1.5x; wlw 0.3 halved every 5 s: 94-98% / 1.9x.
// Downstream (what matters): GLS (30 s) then 70 s of post-feasible SP-SA from
// the hand-off (fp --init-from), dev_tight 16 + hidden4 ws0.05/0.06 4 cases x
// 5 seeds, verifier-checked, legal / final WL over construction, median over
// legal runs: wlw 0: 100% / 1.47x; wlw 0.1: 96% / 1.39x (default, best);
// wlw 0.3 halved every 5 s: 96% / 1.41x; wlw 0.1 + 10 s WL phase: 96% / 1.47x;
// + 20 s WL phase: 95% / 1.42x (plain fp 110 s: 77% / 1.40x). The lower
// hand-off WL of the WL phase does not survive the SA; keep wl_time = 0 here.
// Optional WL phase (wl_time > 0): after the first legal state, SA on WL with
// every overlap-creating move rejected; the best legal placement is returned.
// gls_wl_polish runs only this phase on an already legal placement.
//
// Cost of a step: O(|changed| * n) overlap evaluations against the pair
// matrix (n <= 200 here; ~0.5-1M steps/s, so no grid); see pair_scan().
// Timing is wall-clock based (temperature cycle, WL decay, deadline), so runs
// with the same Rng are not bit-reproducible.
#ifndef FP_GLS_H
#define FP_GLS_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <vector>

#include "gp.h"

namespace fp {

struct GlsOpts {
    double T0 = 0.015, T1 = 2e-5;   // cyclic temperature range (cost units: u^2 area)
    double t_period = 3.0;          // seconds per temperature cycle
    int64_t gls_every = 100000;     // steps between GLS weight bumps
    double gls_inc = 0.25;          // weight increment on overlapping pairs
    double p_pick_ovl = 0.9;        // probability to move a block that overlaps
    int p_shape = 28, p_swap = 6, p_push = 12;   // move mix in percent (rest: snap)
    double shift_frac = 0.06;       // small shift: +-shift_frac * u
    double wl_time = 0;             // after legal: seconds of overlap-free WL descent (0: off)
    double wl_T = 0.3;              // WL phase start temperature, x mean connection WL
    double wlw = 0.1;               // WL bias while infeasible: wlw * dWL / (u * #connections)
    double wlw_half = 0;            // ... halved every wlw_half seconds (0: constant)
    bool verbose = false;
};

struct GlsStats {
    bool legal = false;
    double ttl = -1;                // seconds from the call to the first legal state
    int64_t steps = 0, accepted = 0, bumps = 0;
    double min_overlap = 0;         // smallest raw overlap cost seen
    int64_t wl2_legal = 0, wl2_out = 0;
};

struct GlsLegalizer {
    const Problem &P;
    Rng &rng;
    int n;
    double unit, inv_u, inv_u2;
    std::vector<int> X, Y, Wd, Ht;        // current placement
    std::vector<int> NX, NY, NW, NH;      // tentative (equal to current outside a move)
    std::vector<int> gidx, ngidx;
    std::vector<double> pr, wg;           // pair overlap cost / GLS weight, n*n
    std::vector<int> novl;                // overlap partners per block
    std::vector<int> olist, opos;         // blocks with novl > 0 (index set)
    std::vector<int> mov;                 // movable (non-preplaced) blocks
    int64_t npairs = 0;                   // overlapping pairs
    double total = 0;                     // raw sum of pr (for stats only)
    std::vector<int> ch;                  // changed blocks of the current move
    std::vector<char> inch;
    std::vector<int> todo;

    GlsLegalizer(const Problem &p, Rng &r) : P(p), rng(r), n(p.n) {
        unit = std::sqrt((double)P.W * P.H / std::max(1, n));
        inv_u = 1.0 / unit;
        inv_u2 = inv_u * inv_u;
    }

    double ov(int i, int j) const {   // tentative arrays
        int dx = std::min(NX[i] + NW[i], NX[j] + NW[j]) - std::max(NX[i], NX[j]);
        if (dx <= 0) return 0;
        int dy = std::min(NY[i] + NH[i], NY[j] + NH[j]) - std::max(NY[i], NY[j]);
        if (dy <= 0) return 0;
        return (double)dx * dy * inv_u2 + 0.1 * std::min(dx, dy) * inv_u;
    }

    void clamp(int i) {
        if (P.kind[i] == K_PRE) { NX[i] = P.px[i]; NY[i] = P.py[i]; return; }
        NX[i] = std::max(0, std::min(P.W - NW[i], NX[i]));
        NY[i] = std::max(0, std::min(P.H - NH[i], NY[i]));
        int k = P.bmask[i];
        if (k & lab3::EDGE_LEFT) NX[i] = 0;
        if (k & lab3::EDGE_RIGHT) NX[i] = P.W - NW[i];
        if (k & lab3::EDGE_BOTTOM) NY[i] = 0;
        if (k & lab3::EDGE_TOP) NY[i] = P.H - NH[i];
    }
    bool placed_ok(int i) const {   // outline / boundary / preplaced, no clamping
        if (NX[i] < 0 || NY[i] < 0 || NX[i] + NW[i] > P.W || NY[i] + NH[i] > P.H) return false;
        if (P.kind[i] == K_PRE) return NX[i] == P.px[i] && NY[i] == P.py[i];
        int k = P.bmask[i];
        if ((k & lab3::EDGE_LEFT) && NX[i] != 0) return false;
        if ((k & lab3::EDGE_RIGHT) && NX[i] + NW[i] != P.W) return false;
        if ((k & lab3::EDGE_BOTTOM) && NY[i] != 0) return false;
        if ((k & lab3::EDGE_TOP) && NY[i] + NH[i] != P.H) return false;
        return true;
    }

    void oset(int i, int d) {   // novl[i] += d, keep olist
        int was = novl[i];
        novl[i] += d;
        if (was == 0 && novl[i] > 0) { opos[i] = (int)olist.size(); olist.push_back(i); }
        else if (was > 0 && novl[i] == 0) {
            int last = olist.back(); olist[opos[i]] = last; opos[last] = opos[i]; olist.pop_back(); opos[i] = -1;
        }
    }
    void set_pair(int j, int k, double z) {
        double &a = pr[(size_t)j * n + k];
        if ((a > 0) != (z > 0)) {
            int d = z > 0 ? 1 : -1;
            npairs += d; oset(j, d); oset(k, d);
        }
        total += z - a;
        a = z;
        pr[(size_t)k * n + j] = z;
    }

    void init(const std::vector<lab3::Rect> &r, const std::vector<int> &g) {
        X.resize(n); Y.resize(n); Wd.resize(n); Ht.resize(n);
        gidx = g;
        for (int i = 0; i < n; ++i) {
            if (P.kind[i] == K_SOFT) {
                const ShapeRange &s = P.gshape[P.group[i]];
                Wd[i] = s.widths[gidx[P.group[i]]]; Ht[i] = s.heights[gidx[P.group[i]]];
            } else { Wd[i] = P.fw[i]; Ht[i] = P.fh[i]; }
            // keep the centre of the given rectangle
            X[i] = (int)std::lround(r[i].x + r[i].w / 2.0 - Wd[i] / 2.0);
            Y[i] = (int)std::lround(r[i].y + r[i].h / 2.0 - Ht[i] / 2.0);
        }
        NX = X; NY = Y; NW = Wd; NH = Ht;
        for (int i = 0; i < n; ++i) clamp(i);
        X = NX; Y = NY;
        ngidx = gidx;
        pr.assign((size_t)n * n, 0);
        wg.assign((size_t)n * n, 1.0);
        novl.assign(n, 0); opos.assign(n, -1); olist.clear();
        npairs = 0; total = 0;
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) { double z = ov(i, j); if (z > 0) set_pair(i, j, z); }
        mov.clear();
        for (int i = 0; i < n; ++i) if (P.kind[i] != K_PRE) mov.push_back(i);
        inch.assign(n, 0);
        cb.assign(n, {}); ct.assign(n, {});
        for (int e = 0; e < (int)P.bb.size(); ++e) { cb[P.bb[e].a].push_back(e); if (P.bb[e].b != P.bb[e].a) cb[P.bb[e].b].push_back(e); }
        for (int e = 0; e < (int)P.bt.size(); ++e) ct[P.bt[e].a].push_back(e);
    }

    void add_ch(int j) { if (!inch[j]) { inch[j] = 1; ch.push_back(j); } }
    void revert() {
        for (int j : ch) { NX[j] = X[j]; NY[j] = Y[j]; NW[j] = Wd[j]; NH[j] = Ht[j]; inch[j] = 0; }
        ch.clear();
    }
    // Weighted delta of the tentative move: every pair touching a changed block once.
    double pair_scan(double &raw) const {
        double d = 0;
        raw = 0;
        for (int j : ch) {
            const double *prj = &pr[(size_t)j * n], *wgj = &wg[(size_t)j * n];
            for (int k = 0; k < n; ++k) {
                if (k == j || (inch[k] && k > j)) continue;
                double z = ov(j, k), e = z - prj[k];
                if (e != 0) { d += e * wgj[k]; raw += e; }
            }
        }
        return d;
    }
    void commit() {
        for (int j : ch) for (int k = 0; k < n; ++k) {
            if (k == j || (inch[k] && k > j)) continue;
            double z = ov(j, k);
            if (z != pr[(size_t)j * n + k]) set_pair(j, k, z);
        }
        for (int j : ch) { X[j] = NX[j]; Y[j] = NY[j]; Wd[j] = NW[j]; Ht[j] = NH[j]; inch[j] = 0; }
        ch.clear();
        gidx = ngidx;
    }

    int partner(int i, int tries) {   // a block overlapping i (or random)
        int j = rng.below(n);
        for (int t = 0; t < tries && (j == i || pr[(size_t)i * n + j] == 0); ++t) j = rng.below(n);
        return j;
    }

    // Push chain along axis (0: x, 1: y) by delta from root. False: infeasible.
    bool push(int root, int axis, int delta) {
        if (delta == 0) return false;
        int dir = delta > 0 ? 1 : -1;
        add_ch(root);
        (axis ? NY : NX)[root] += delta;
        todo.assign(1, root);
        for (size_t z = 0; z < todo.size(); ++z) {
            if (todo.size() > (size_t)4 * n) return false;
            int i = todo[z];
            if (!placed_ok(i)) return false;
            for (int j = 0; j < n; ++j) {
                if (j == i || ov(i, j) == 0) continue;
                if (P.kind[j] == K_PRE) return false;
                int nv = axis == 0 ? (dir > 0 ? NX[i] + NW[i] : NX[i] - NW[j]) : (dir > 0 ? NY[i] + NH[i] : NY[i] - NH[j]);
                int &c = (axis ? NY : NX)[j];
                if (dir > 0 ? nv <= c : nv >= c) continue;   // never pull back
                c = nv;
                add_ch(j);
                todo.push_back(j);
            }
        }
        return true;
    }

    // Build a tentative move; false: nothing to evaluate.
    bool propose() {
        int i;
        if (npairs > 0 && rng.uni() < O.p_pick_ovl) {
            i = olist[rng.below((int)olist.size())];
            for (int t = 0; t < 8 && P.kind[i] == K_PRE; ++t) i = olist[rng.below((int)olist.size())];
            if (P.kind[i] == K_PRE) i = mov[rng.below((int)mov.size())];
        } else {
            i = mov[rng.below((int)mov.size())];
        }
        int a = rng.below(100);
        if (a < O.p_shape && P.kind[i] == K_SOFT) {
            int g = P.group[i];
            int m = (int)P.gshape[g].widths.size();
            if (m < 2) return false;
            int rad = std::max(2, (int)(m * (rng.below(10) == 0 ? 0.15 : 0.025)));
            int ni = std::max(0, std::min(m - 1, gidx[g] + rng.below(2 * rad + 1) - rad));
            if (ni == gidx[g]) return false;
            ngidx[g] = ni;
            int nw = P.gshape[g].widths[ni], nh = P.gshape[g].heights[ni];
            for (int j : P.gmem[g]) {
                add_ch(j);
                int mode = rng.below(3);
                if (mode == 0) { NX[j] += (Wd[j] - nw) / 2; NY[j] += (Ht[j] - nh) / 2; }
                else if (mode == 1) { NX[j] += Wd[j] - nw; NY[j] += Ht[j] - nh; }
                NW[j] = nw; NH[j] = nh;
                clamp(j);
            }
            return true;
        }
        if (a < O.p_shape + O.p_swap) {
            int j = mov[rng.below((int)mov.size())];
            if (j == i) return false;
            add_ch(i); add_ch(j);
            int cxi = 2 * X[i] + Wd[i], cyi = 2 * Y[i] + Ht[i], cxj = 2 * X[j] + Wd[j], cyj = 2 * Y[j] + Ht[j];
            NX[i] = (cxj - Wd[i]) / 2; NY[i] = (cyj - Ht[i]) / 2;
            NX[j] = (cxi - Wd[j]) / 2; NY[j] = (cyi - Ht[j]) / 2;
            clamp(i); clamp(j);
            return true;
        }
        if (a >= 100 - O.p_push) {
            int j = partner(i, 40);
            int axis = rng.below(2), dir = rng.below(2) ? 1 : -1;
            int delta = dir * (1 + rng.below(std::max(2, (int)(unit * 0.1))));
            if (j != i && pr[(size_t)i * n + j] > 0) {
                if (axis == 0) delta = dir > 0 ? X[j] + Wd[j] - X[i] : X[j] - Wd[i] - X[i];
                else delta = dir > 0 ? Y[j] + Ht[j] - Y[i] : Y[j] - Ht[i] - Y[i];
            }
            return push(i, axis, delta);
        }
        int j = rng.uni() < 0.65 ? partner(i, 30) : rng.below(n);
        add_ch(i);
        int axis = rng.below(2), mode = rng.below(6);
        int &c = axis ? NY[i] : NX[i];
        const int pj = axis ? Y[j] : X[j], sj = axis ? Ht[j] : Wd[j], si = axis ? Ht[i] : Wd[i], L = axis ? P.H : P.W;
        switch (mode) {
            case 0: c = pj + sj; break;
            case 1: c = pj - si; break;
            case 2: c = pj; break;
            case 3: c = pj + sj - si; break;
            case 4: c += (int)std::lround((rng.uni() - 0.5) * 2 * O.shift_frac * unit); break;
            default: c = rng.below(2) ? 0 : L - si; break;
        }
        clamp(i);
        return true;
    }

    GlsOpts O;

    // ---- optional WL phase (legal -> legal): moves that create any overlap are
    // rejected; SA on the doubled-centre WL of the connections of changed blocks.
    std::vector<std::vector<int>> cb;   // per block: indices into P.bb
    std::vector<std::vector<int>> ct;   // per block: indices into P.bt
    int64_t conn_wl2(const BBConn &c) const {
        int64_t dx = (2LL * NX[c.a] + NW[c.a]) - (2LL * NX[c.b] + NW[c.b]);
        int64_t dy = (2LL * NY[c.a] + NH[c.a]) - (2LL * NY[c.b] + NH[c.b]);
        return c.w * ((dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy));
    }
    int64_t term_wl2(const BTConn &c) const {
        int64_t dx = (2LL * NX[c.a] + NW[c.a]) - c.tx2, dy = (2LL * NY[c.a] + NH[c.a]) - c.ty2;
        return c.w * ((dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy));
    }
    int64_t wl_delta() {   // tentative minus current, connections of changed blocks
        int64_t d = 0;
        auto cur = [&](auto f) {   // evaluate f on the current placement
            for (int j : ch) { std::swap(NX[j], X[j]); std::swap(NY[j], Y[j]); std::swap(NW[j], Wd[j]); std::swap(NH[j], Ht[j]); }
            int64_t v = f();
            for (int j : ch) { std::swap(NX[j], X[j]); std::swap(NY[j], Y[j]); std::swap(NW[j], Wd[j]); std::swap(NH[j], Ht[j]); }
            return v;
        };
        auto sum = [&]() {
            int64_t s = 0;
            for (int j : ch) {
                for (int e : cb[j]) { const BBConn &c = P.bb[e]; int o = c.a == j ? c.b : c.a; if (!inch[o] || o < j) s += conn_wl2(c); }
                for (int e : ct[j]) s += term_wl2(P.bt[e]);
            }
            return s;
        };
        d = sum() - cur(sum);
        return d;
    }
    void wl_phase(double seconds, GlsStats &st) {
        using Clock = std::chrono::steady_clock;
        auto t0 = Clock::now();
        int ncon = std::max<int>(1, (int)(P.bb.size() + P.bt.size()));
        int64_t wl = 0;
        for (const BBConn &c : P.bb) wl += conn_wl2(c);
        for (const BTConn &c : P.bt) wl += term_wl2(c);
        int64_t bestwl = wl;
        std::vector<int> bX = X, bY = Y, bg = gidx;
        double T = 1, now = 0;
        int64_t steps = 0;
        while (true) {
            if ((steps & 1023) == 0) {
                now = std::chrono::duration<double>(Clock::now() - t0).count();
                if (now >= seconds) break;
                double ph = std::fmod(now, 1.0), scale = (double)wl / ncon;
                T = O.wl_T * scale * std::pow(0.01, ph) * (1.0 - now / seconds) + 1e-9;
            }
            ++steps;
            if (!propose()) { revert(); ngidx = gidx; continue; }
            double raw;
            pair_scan(raw);
            bool ok = raw <= 0;
            if (ok) for (int j : ch) if (!placed_ok(j)) { ok = false; break; }
            if (!ok) { revert(); ngidx = gidx; continue; }
            int64_t d = wl_delta();
            if (d <= 0 || rng.uni() < std::exp(-(double)d / T)) {
                commit();
                wl += d;
                if (wl < bestwl) { bestwl = wl; bX = X; bY = Y; bg = gidx; }
            } else { revert(); ngidx = gidx; }
        }
        // restore the best legal placement
        NX = X = bX; NY = Y = bY; gidx = ngidx = bg;
        for (int i = 0; i < n; ++i) {
            if (P.kind[i] == K_SOFT) { NW[i] = Wd[i] = P.gshape[P.group[i]].widths[gidx[P.group[i]]]; NH[i] = Ht[i] = P.gshape[P.group[i]].heights[gidx[P.group[i]]]; }
        }
        st.steps += steps;
    }

    // Returns true when legal. r / g: start (in) and result (out; the legal
    // placement, or the least-overlap one seen).
    bool run(std::vector<lab3::Rect> &r, std::vector<int> &g, double seconds, GlsStats &st, const GlsOpts &opts) {
        using Clock = std::chrono::steady_clock;
        O = opts;
        auto t0 = Clock::now();
        auto el = [&]() { return std::chrono::duration<double>(Clock::now() - t0).count(); };
        init(r, g);
        std::vector<int> bX = X, bY = Y, bg = gidx;
        double best = total;
        double T = O.T0, now = 0, wls = 0;
        st = GlsStats();
        const double wls0 = O.wlw / (2.0 * unit * std::max<size_t>(1, P.bb.size() + P.bt.size()));   // per doubled WL unit
        while (npairs > 0 && !mov.empty()) {
            if ((st.steps & 1023) == 0) {
                now = el();
                if (now >= seconds) break;
                double ph = std::fmod(now, O.t_period) / O.t_period;
                T = O.T0 * std::pow(O.T1 / O.T0, ph);
                wls = O.wlw_half > 0 ? wls0 * std::pow(0.5, now / O.wlw_half) : wls0;
            }
            ++st.steps;
            if (st.steps % O.gls_every == 0) {
                ++st.bumps;
                for (int j : olist) for (int k = 0; k < n; ++k)
                    if (k > j && pr[(size_t)j * n + k] > 0) { wg[(size_t)j * n + k] += O.gls_inc; wg[(size_t)k * n + j] += O.gls_inc; }
                if (O.verbose)
                    std::fprintf(stderr, "GLS t=%.2f steps=%lld pairs=%lld total=%.4f best=%.4f\n", el(), (long long)st.steps,
                                 (long long)npairs, total, best);
            }
            if (!propose()) { revert(); ngidx = gidx; continue; }
            double raw, d = pair_scan(raw);
            if (wls > 0) d += wls * (double)wl_delta();
            if (d <= 0 || rng.uni() < std::exp(-d / T)) {
                commit();
                ++st.accepted;
                if (total < best - 1e-12 || npairs == 0) { best = total; bX = X; bY = Y; bg = gidx; }
            } else {
                revert();
                ngidx = gidx;
            }
        }
        st.min_overlap = npairs == 0 ? 0 : best;
        st.legal = npairs == 0;
        if (st.legal) {
            st.ttl = el();
            if (O.wl_time > 0) {
                O.p_pick_ovl = 0;   // no overlaps left: uniform block choice
                st.wl2_legal = 0;
                for (const BBConn &c : P.bb) st.wl2_legal += conn_wl2(c);
                for (const BTConn &c : P.bt) st.wl2_legal += term_wl2(c);
                wl_phase(std::min(O.wl_time, std::max(0.0, seconds - el())), st);
            }
            bX = X; bY = Y; bg = gidx;
        }
        g = bg;
        r.assign(n, lab3::Rect());
        for (int i = 0; i < n; ++i) {
            int w = P.kind[i] == K_SOFT ? P.gshape[P.group[i]].widths[bg[P.group[i]]] : P.fw[i];
            int h = P.kind[i] == K_SOFT ? P.gshape[P.group[i]].heights[bg[P.group[i]]] : P.fh[i];
            r[i] = lab3::Rect{bX[i], bY[i], w, h};
        }
        st.wl2_out = lab3::wirelength2(P.inst, r);
        if (!(st.legal && O.wl_time > 0)) st.wl2_legal = st.wl2_out;
        return st.legal;
    }
};

// One-call API. `seconds`: wall budget for this call.
inline bool gls_legalize(const Problem &P, std::vector<lab3::Rect> &r, std::vector<int> &gidx, double seconds, Rng &rng,
                         GlsStats &st, const GlsOpts &opts = GlsOpts()) {
    GlsLegalizer L(P, rng);
    return L.run(r, gidx, seconds, st, opts);
}

// Post-legal polish: only the WL phase, started from a legal placement r
// (soft shapes: gidx per shape group must index the exact widths of r; a
// height above hmin(width) is shrunk around the centre, which keeps it legal).
// Returns false (r untouched) when the start is not legal under the gls
// invariants; else r / gidx get the best legal placement seen (never worse).
inline bool gls_wl_polish(const Problem &P, std::vector<lab3::Rect> &r, std::vector<int> &gidx, double seconds, Rng &rng,
                          GlsStats &st, const GlsOpts &opts = GlsOpts()) {
    GlsLegalizer L(P, rng);
    L.O = opts;
    L.init(r, gidx);
    st = GlsStats();
    if (L.npairs != 0) return false;
    for (int i = 0; i < P.n; ++i) if (!L.placed_ok(i)) return false;
    L.O.p_pick_ovl = 0;
    st.legal = true;
    st.ttl = 0;
    for (const BBConn &c : P.bb) st.wl2_legal += L.conn_wl2(c);
    for (const BTConn &c : P.bt) st.wl2_legal += L.term_wl2(c);
    L.wl_phase(seconds, st);
    gidx = L.gidx;
    r.assign(P.n, lab3::Rect());
    for (int i = 0; i < P.n; ++i) r[i] = lab3::Rect{L.X[i], L.Y[i], L.Wd[i], L.Ht[i]};
    st.wl2_out = lab3::wirelength2(P.inst, r);
    return true;
}

}  // namespace fp

#endif  // FP_GLS_H
