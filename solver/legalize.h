// legalize.h -- topology-preserving legalizer for analytical placements.
//
//   legalize(P, cx, cy, gidx, rects, g, stats)
//
// Input: block centres (cx, cy; preplaced centres are ignored and taken from
// the input, boundary blocks are snapped to their edges) and a soft shape per
// shape group (gidx, index into P.gshape[g].widths). Output: a legal placement
// (gate checks 3-8, every boundary edge touched, every MIB group one shape)
// with exact minimum WL for its H/V constraint graphs, or failure.
//
// 1. Relations. Every pair gets one axis: H if its normalized x displacement
//    |dcx| / ((w_i + w_j)/2) is >= the y one, else V. The direction on an axis
//    always follows a per-axis key order (centre, with LEFT/BOTTOM blocks
//    before and RIGHT/TOP blocks after everything), so both graphs are DAGs
//    by construction and the key order is a topological order. A relation
//    that is impossible outright (an edge into a LEFT block, out of a RIGHT
//    block, a block that does not fit beside a preplaced block, two preplaced
//    blocks against their geometry) takes the other axis if that one is not
//    also impossible.
// 2. Feasibility. Per axis, ASAP longest path (not clamped) with bounds
//    lo/hi (preplaced and boundary pinned, others [0, L - size]) and ALAP.
//    excess = sum of ASAP - hi > 0 over nodes; tns = sum over constraint
//    edges i -> j of ASAP_i + s_i - ALAP_j > 0 (sees parallel critical paths).
//    Feasible <=> excess == 0.
// 3. Repair (greedy descent on excess + tns): per axis, trace the critical
//    chain of the worst violator and evaluate candidates: flip a chain pair to
//    the other axis; relocate a chain block out of its band (mode 0/1/2);
//    teleport a chain block's centre by +-0.6/1.2 of its size and re-derive
//    its relations; narrow a soft chain block (MIB groups share the shape);
//    plus the top gflip_k negative-slack pairs whose flip has slack on the
//    other axis. The best improving candidate is applied; stop when none
//    improves or after max_iter iterations.
// 4. Exact WL LP (mcf::polish_topology) on the transitively reduced graphs,
//    starting from the ASAP placement (which is legal once excess == 0).
//
// STATUS (measured, see the test_legalize report): converges on near-packings
// (construction centres, also jittered by +-40% / +-80% of block size) but
// stalls on fp's GlobalPlacer leaf centres on every preset, so it is not
// wired into fp.cpp. LG_DEBUG=1 dumps the critical chains where it stalls.
#ifndef FP_LEGALIZE_H
#define FP_LEGALIZE_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <string>
#include <vector>

#include "gp.h"

namespace fp {

struct LegalizeOpts {
    int max_iter = 0;          // 0: 10 * n
    bool shape = true;         // allow soft narrowing candidates
    bool flip = true;          // allow pair flip candidates
    bool reloc = true;         // allow band relocation candidates
    bool gflip = true;         // global negative-slack flip candidates
    bool teleport = true;      // move a chain block's centre and re-derive its relations
    int gflip_k = 20;          // ... evaluated per iteration
    int spread_iters = 0;      // overlap-removal iterations on the centres first
    bool lp = true;            // exact WL LP on the final graphs
};

struct LegalizeStats {
    bool ok = false;
    std::string why;
    int iters = 0, flips = 0, shapes = 0, relocs = 0, teleports = 0;
    int64_t excess0 = 0, excess = 0;   // initial / final total excess
    int64_t wl2_asap = 0, wl2 = 0;     // 2*WL of the ASAP placement / the output
    bool lp_ok = false;
    double ms = 0, ms_lp = 0;
};

class Legalizer {
public:
    Legalizer(const Problem &P) : P_(P), n_(P.n) {}

    bool run(const std::vector<double> &cx_in, const std::vector<double> &cy_in, std::vector<int> &gidx,
             std::vector<lab3::Rect> &out, mcf::ConstraintGraphs *gout, LegalizeStats &st,
             const LegalizeOpts &o = LegalizeOpts()) {
        auto t0 = std::chrono::steady_clock::now();
        st = LegalizeStats();
        const int n = n_;
        gidx_ = gidx;
        w_.assign(n, 0); h_.assign(n, 0);
        for (int i = 0; i < n; ++i) if (P_.kind[i] != K_SOFT) { w_[i] = P_.fw[i]; h_[i] = P_.fh[i]; }
        for (int g = 0; g < (int)P_.gmem.size(); ++g) set_group(g, gidx_[g]);
        // Centres: preplaced from the input, boundary blocks snapped to their edges.
        cx_ = cx_in; cy_ = cy_in;
        for (int i = 0; i < n; ++i) {
            if (P_.kind[i] == K_PRE) { cx_[i] = P_.px[i] + P_.fw[i] / 2.0; cy_[i] = P_.py[i] + P_.fh[i] / 2.0; continue; }
            int k = P_.bmask[i];
            if (k & lab3::EDGE_LEFT) cx_[i] = w_[i] / 2.0;
            if (k & lab3::EDGE_RIGHT) cx_[i] = P_.W - w_[i] / 2.0;
            if (k & lab3::EDGE_BOTTOM) cy_[i] = h_[i] / 2.0;
            if (k & lab3::EDGE_TOP) cy_[i] = P_.H - h_[i] / 2.0;
        }
        if (o.spread_iters > 0) spread_overlap(o.spread_iters);
        // Key orders (topological orders of both graphs).
        for (int a = 0; a < 2; ++a) {
            key_[a].resize(n);
            for (int i = 0; i < n; ++i) key_[a][i] = make_key(a, i);
            ord_[a].resize(n);
            std::iota(ord_[a].begin(), ord_[a].end(), 0);
            const std::vector<double> &key = key_[a];
            std::sort(ord_[a].begin(), ord_[a].end(), [&](int p, int q) { return key[p] != key[q] ? key[p] < key[q] : p < q; });
            rk_[a].assign(n, 0);
            for (int r = 0; r < n; ++r) rk_[a][ord_[a][r]] = r;
        }
        // Initial relations.
        ax_.assign((size_t)n * n, 0);
        flipped_.assign((size_t)n * n, 0);
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) classify(i, j);
        int max_iter = o.max_iter > 0 ? o.max_iter : 10 * n;
        Eval cur = evaluate(true);
        st.excess0 = cur.excess;
        std::vector<uint8_t> ax_save;
        std::vector<int> ord_save, rk_save;
        while (cur.excess > 0 && st.iters < max_iter) {
            ++st.iters;
            // Chains: the worst violator on each axis.
            Eval best = cur;
            int bkind = -1, bi = -1, bj = -1, bidx = -1, bax = -1;
            double tpx = 0, tpy = 0;
            for (int va = 0; va < 2; ++va) {
                int vv = -1;
                int64_t vmax = 0;
                for (int i = 0; i < n; ++i)
                    if (exc_[va][i] > vmax) { vmax = exc_[va][i]; vv = i; }
                if (vv < 0) continue;
                std::vector<int> chain;
                for (int u = vv; u >= 0; u = pred_[va][u]) chain.push_back(u);
                std::reverse(chain.begin(), chain.end());
                // first violator on the chain: the local excess there
                int64_t E = vmax;
                for (int u : chain) if (exc_[va][u] > 0) { E = exc_[va][u]; break; }
                if (o.flip) {
                    for (size_t k = 0; k + 1 < chain.size(); ++k) {
                        int i = chain[k], j = chain[k + 1];
                        if (flipped_[(size_t)i * n + j] >= 3 || impossible(i, j, 1 - va)) continue;
                        if (P_.kind[i] == K_PRE && P_.kind[j] == K_PRE) continue;
                        toggle(i, j);
                        Eval e = evaluate(false);
                        toggle(i, j);
                        if (e.cost() < best.cost()) { best = e; bkind = 0; bi = i; bj = j; }
                    }
                }
                if (o.reloc) {
                    for (int u : chain) {
                        if (P_.kind[u] == K_PRE) continue;
                        for (int mode = 0; mode < 3; ++mode) {
                            ax_save = ax_; ord_save = ord_[1 - va]; rk_save = rk_[1 - va];
                            bool ok = relocate(u, va, mode);
                            Eval e;
                            if (ok) e = evaluate(false);
                            ax_.swap(ax_save); ord_[1 - va].swap(ord_save); rk_[1 - va].swap(rk_save);
                            if (ok && e.cost() < best.cost()) { best = e; bkind = 2; bi = u; bax = va; bidx = mode; }
                        }
                    }
                }
                if (o.teleport) {
                    static const double fo[5] = {0, -0.6, 0.6, -1.2, 1.2};
                    for (int u : chain) {
                        if (P_.kind[u] == K_PRE) continue;
                        double ox = cx_[u], oy = cy_[u];
                        for (int fx = 0; fx < 5; ++fx)
                            for (int fy = 0; fy < 5; ++fy) {
                                if (!fx && !fy) continue;
                                ax_save = ax_;
                                std::vector<int> o0 = ord_[0], o1 = ord_[1], r0 = rk_[0], r1 = rk_[1];
                                double k0 = key_[0][u], k1 = key_[1][u];
                                teleport(u, ox + fo[fx] * w_[u], oy + fo[fy] * h_[u]);
                                Eval e = evaluate(false);
                                ax_.swap(ax_save); ord_[0].swap(o0); ord_[1].swap(o1); rk_[0].swap(r0); rk_[1].swap(r1);
                                cx_[u] = ox; cy_[u] = oy; key_[0][u] = k0; key_[1][u] = k1;
                                if (e.cost() < best.cost()) { best = e; bkind = 3; bi = u; tpx = ox + fo[fx] * w_[u]; tpy = oy + fo[fy] * h_[u]; }
                            }
                    }
                }
                if (o.shape) {
                    std::vector<int> seen;
                    for (int u : chain) {
                        if (P_.kind[u] != K_SOFT) continue;
                        int g = P_.group[u];
                        if (std::find(seen.begin(), seen.end(), g) != seen.end()) continue;
                        seen.push_back(g);
                        int old = gidx_[g];
                        int64_t ds[4] = {E, (E + 1) / 2, (E + 3) / 4, E < 10 ? 1 : 0};
                        int lastc = -1;
                        for (int64_t d : ds) {
                            if (d <= 0) continue;
                            int c = narrow_index(g, va, d);
                            if (c < 0 || c == old || c == lastc) continue;
                            lastc = c;
                            set_group(g, c);
                            Eval e = evaluate(false);
                            set_group(g, old);
                            if (e.cost() < best.cost()) { best = e; bkind = 1; bi = g; bidx = c; }
                        }
                    }
                }
            }
            if (o.gflip) {
                // Every negative-slack constraint pair whose flip to the other axis
                // has non-negative slack there, most negative first; top K evaluated.
                std::vector<std::pair<int64_t, std::pair<int, int>>> cand;
                for (int a = 0; a < 2; ++a) {
                    const std::vector<int> &sz = a == 0 ? w_ : h_, &sb = a == 0 ? h_ : w_;
                    const int b = 1 - a;
                    for (int i = 0; i < n; ++i)
                        for (int j = 0; j < n; ++j) {
                            if (i == j || ax_[(size_t)i * n + j] != a || rk_[a][i] > rk_[a][j]) continue;
                            int64_t sl = la_[a][j] - e_[a][i] - sz[i];
                            if (sl >= 0) continue;
                            if (P_.kind[i] == K_PRE && P_.kind[j] == K_PRE) continue;
                            if (flipped_[(size_t)i * n + j] >= 3 || impossible(i, j, b)) continue;
                            int u = rk_[b][i] < rk_[b][j] ? i : j, v = u == i ? j : i;
                            int64_t m = la_[b][v] - e_[b][u] - sb[u];
                            if (m < 0) continue;
                            cand.push_back({sl, {i, j}});
                        }
                }
                std::sort(cand.begin(), cand.end());
                int K = std::min((int)cand.size(), o.gflip_k);
                for (int c = 0; c < K; ++c) {
                    int i = cand[c].second.first, j = cand[c].second.second;
                    toggle(i, j);
                    Eval e = evaluate(false);
                    toggle(i, j);
                    if (e.cost() < best.cost()) { best = e; bkind = 0; bi = i; bj = j; }
                }
            }
            if (bkind < 0) {
                if (std::getenv("LG_DEBUG")) dump_stuck();
                break;
            }
            if (bkind == 0) { toggle(bi, bj); ++flipped_[(size_t)bi * n + bj]; ++flipped_[(size_t)bj * n + bi]; ++st.flips; }
            else if (bkind == 1) { set_group(bi, bidx); ++st.shapes; }
            else if (bkind == 2) { relocate(bi, bax, bidx); ++st.relocs; }
            else { teleport(bi, tpx, tpy); ++st.teleports; }
            cur = evaluate(true);
        }
        st.excess = cur.excess;
        auto ms_since = [](std::chrono::steady_clock::time_point t) {
            return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count();
        };
        if (cur.excess > 0) { st.why = "excess"; st.ms = ms_since(t0); return false; }
        // ASAP placement (feasible: unclamped == clamped).
        out.assign(n, lab3::Rect());
        for (int i = 0; i < n; ++i) out[i] = lab3::Rect{e_[0][i], e_[1][i], w_[i], h_[i]};
        st.wl2_asap = st.wl2 = lab3::wirelength2(P_.inst, out);
        mcf::ConstraintGraphs g;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                int a = ax_[(size_t)i * n + j];
                if (rk_[a][i] < rk_[a][j]) (a == 0 ? g.h : g.v).emplace_back(i, j);
            }
        bool cyc = false;
        g.h = mcf::transitive_reduction(n, g.h, &cyc);
        if (cyc) { st.why = "cycle h"; st.ms = ms_since(t0); return false; }
        g.v = mcf::transitive_reduction(n, g.v, &cyc);
        if (cyc) { st.why = "cycle v"; st.ms = ms_since(t0); return false; }
        if (o.lp) {
            auto t1 = std::chrono::steady_clock::now();
            mcf::PolishInfo info;
            std::vector<lab3::Rect> r = mcf::polish_topology(P_.inst, out, g, &info);
            st.ms_lp = ms_since(t1);
            st.lp_ok = info.ok;
            if (info.ok) { out = r; st.wl2 = info.wl2_out; }
        }
        if (gout) *gout = g;
        gidx = gidx_;
        st.ok = true;
        st.ms = ms_since(t0);
        return true;
    }

    // Optional leaf regions (GlobalPlacer leafw/leafh): movable blocks use
    // them instead of their shapes to normalize displacements.
    void set_leaf(const std::vector<double> *lw, const std::vector<double> *lh) { lw_ = lw; lh_ = lh; }

private:
    const Problem &P_;
    const std::vector<double> *lw_ = nullptr, *lh_ = nullptr;
    double ew(int i) const { return lw_ && P_.kind[i] != K_PRE ? (*lw_)[i] : w_[i]; }
    double eh(int i) const { return lh_ && P_.kind[i] != K_PRE ? (*lh_)[i] : h_[i]; }
    int n_;
    std::vector<int> gidx_, w_, h_;
    std::vector<double> cx_, cy_;
    std::vector<int> ord_[2], rk_[2];
    std::vector<uint8_t> ax_, flipped_;
    std::vector<int64_t> e_[2], exc_[2];
    std::vector<int> pred_[2];
    std::vector<int64_t> la_[2];

    // Pairwise overlap removal on the centres (shapes fixed): each
    // overlapping pair is pushed apart along the axis of smaller normalized
    // overlap, split by movability (preplaced and edge-pinned coordinates do
    // not move), then clamped into the outline.
    void spread_overlap(int iters) {
        const int n = n_;
        std::vector<char> mx(n), my(n);
        for (int i = 0; i < n; ++i) {
            int k = P_.bmask[i];
            mx[i] = P_.kind[i] != K_PRE && !(k & (lab3::EDGE_LEFT | lab3::EDGE_RIGHT));
            my[i] = P_.kind[i] != K_PRE && !(k & (lab3::EDGE_BOTTOM | lab3::EDGE_TOP));
        }
        std::vector<double> dx(n), dy(n);
        for (int it = 0; it < iters; ++it) {
            std::fill(dx.begin(), dx.end(), 0.0);
            std::fill(dy.begin(), dy.end(), 0.0);
            double tot = 0;
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j) {
                    double ox = (w_[i] + w_[j]) / 2.0 - std::fabs(cx_[i] - cx_[j]);
                    double oy = (h_[i] + h_[j]) / 2.0 - std::fabs(cy_[i] - cy_[j]);
                    if (ox <= 0 || oy <= 0) continue;
                    tot += ox * oy;
                    bool hx = ox / (w_[i] + w_[j]) <= oy / (h_[i] + h_[j]);
                    if (hx && !mx[i] && !mx[j]) hx = false;
                    if (!hx && !my[i] && !my[j]) hx = true;
                    if (hx) {
                        double sgn = cx_[i] < cx_[j] || (cx_[i] == cx_[j] && i < j) ? -1 : 1;
                        double fi = mx[i] ? (mx[j] ? 0.5 : 1.0) : 0.0, fj = mx[j] ? 1.0 - fi : 0.0;
                        dx[i] += sgn * ox * fi; dx[j] -= sgn * ox * fj;
                    } else {
                        double sgn = cy_[i] < cy_[j] || (cy_[i] == cy_[j] && i < j) ? -1 : 1;
                        double fi = my[i] ? (my[j] ? 0.5 : 1.0) : 0.0, fj = my[j] ? 1.0 - fi : 0.0;
                        dy[i] += sgn * oy * fi; dy[j] -= sgn * oy * fj;
                    }
                }
            if (tot == 0) break;
            for (int i = 0; i < n; ++i) {
                if (mx[i]) cx_[i] = std::min(P_.W - w_[i] / 2.0, std::max(w_[i] / 2.0, cx_[i] + o_step_ * dx[i]));
                if (my[i]) cy_[i] = std::min(P_.H - h_[i] / 2.0, std::max(h_[i] / 2.0, cy_[i] + o_step_ * dy[i]));
            }
        }
    }
    double o_step_ = 0.5;
    std::vector<double> key_[2];

    double make_key(int a, int i) const {
        int k = P_.bmask[i];
        double c = a == 0 ? cx_[i] : cy_[i];
        int lowE = a == 0 ? lab3::EDGE_LEFT : lab3::EDGE_BOTTOM, highE = a == 0 ? lab3::EDGE_RIGHT : lab3::EDGE_TOP;
        if (P_.kind[i] != K_PRE && (k & lowE)) c -= 1e12;
        if (P_.kind[i] != K_PRE && (k & highE)) c += 1e12;
        return c;
    }
    // Axis of pair (i, j) from the current centres (see header, step 1).
    void classify(int i, int j) {
        const int n = n_;
        int a;
        if (P_.kind[i] == K_PRE && P_.kind[j] == K_PRE) {
            int64_t gx = std::max(P_.px[i], P_.px[j]) - std::min(P_.px[i] + P_.fw[i], P_.px[j] + P_.fw[j]);
            int64_t gy = std::max(P_.py[i], P_.py[j]) - std::min(P_.py[i] + P_.fh[i], P_.py[j] + P_.fh[j]);
            a = gx >= gy ? 0 : 1;
        } else {
            double sh = std::fabs(cx_[i] - cx_[j]) / std::max(1.0, (ew(i) + ew(j)) / 2.0);
            double sv = std::fabs(cy_[i] - cy_[j]) / std::max(1.0, (eh(i) + eh(j)) / 2.0);
            a = sh >= sv ? 0 : 1;
            if (impossible(i, j, a) && !impossible(i, j, 1 - a)) a = 1 - a;
        }
        ax_[(size_t)i * n + j] = ax_[(size_t)j * n + i] = (uint8_t)a;
    }
    // Re-sort u in the key order of axis a after its key changed.
    void rekey(int u, int a) {
        const int n = n_;
        std::vector<int> &ord = ord_[a];
        const std::vector<double> &key = key_[a];
        ord.erase(ord.begin() + rk_[a][u]);
        auto less = [&](int p, int q) { return key[p] != key[q] ? key[p] < key[q] : p < q; };
        ord.insert(std::lower_bound(ord.begin(), ord.end(), u, less), u);
        for (int r = 0; r < n; ++r) rk_[a][ord[r]] = r;
    }
    // Move movable block u's centre to (px, py) (edge-pinned coordinates
    // stay) and re-derive all of u's relations.
    void teleport(int u, double px, double py) {
        int k = P_.bmask[u];
        if (!(k & (lab3::EDGE_LEFT | lab3::EDGE_RIGHT))) cx_[u] = std::min(P_.W - w_[u] / 2.0, std::max(w_[u] / 2.0, px));
        if (!(k & (lab3::EDGE_BOTTOM | lab3::EDGE_TOP))) cy_[u] = std::min(P_.H - h_[u] / 2.0, std::max(h_[u] / 2.0, py));
        for (int a = 0; a < 2; ++a) { key_[a][u] = make_key(a, u); rekey(u, a); }
        for (int j = 0; j < n_; ++j) if (j != u) classify(u, j);
    }

    void dump_stuck() const {
        const int n = n_;
        for (int a = 0; a < 2; ++a) {
            int vv = -1; int64_t vmax = 0, tot = 0; int cnt = 0;
            for (int i = 0; i < n; ++i) { if (exc_[a][i] > 0) { ++cnt; tot += exc_[a][i]; } if (exc_[a][i] > vmax) { vmax = exc_[a][i]; vv = i; } }
            std::fprintf(stderr, "STUCK axis=%d violators=%d sum=%lld worst=%d (%lld)\n", a, cnt, (long long)tot, vv, (long long)vmax);
            if (vv < 0) continue;
            std::vector<int> chain;
            for (int u = vv; u >= 0; u = pred_[a][u]) chain.push_back(u);
            std::reverse(chain.begin(), chain.end());
            for (int u : chain) {
                int64_t lo, hi; bounds(a, u, lo, hi);
                int64_t nsl_o = la_[1 - a][u] - e_[1 - a][u];
                std::fprintf(stderr, "  %3d kind=%d bm=%d w=%4d h=%4d c=(%5.0f,%5.0f) e=%5lld lo=%5lld hi=%5lld  other: e=%5lld slack=%lld\n", u, P_.kind[u], P_.bmask[u], w_[u], h_[u],
                             cx_[u], cy_[u], (long long)e_[a][u], (long long)lo, (long long)hi, (long long)e_[1 - a][u], (long long)nsl_o);
            }
        }
    }

    struct Eval {
        int64_t excess = 0, tns = 0;
        int64_t cost() const { return excess + tns; }
    };

    void set_group(int g, int idx) {
        gidx_[g] = idx;
        int w = P_.gshape[g].widths[idx], h = P_.gshape[g].heights[idx];
        for (int b : P_.gmem[g]) { w_[b] = w; h_[b] = h; }
    }
    void toggle(int i, int j) {
        uint8_t &a = ax_[(size_t)i * n_ + j];
        a ^= 1;
        ax_[(size_t)j * n_ + i] = a;
    }

    // Smallest change of group g's index that reduces its extent on axis a
    // by at least d (or as much as possible); -1 if none.
    int narrow_index(int g, int a, int64_t d) const {
        const ShapeRange &s = P_.gshape[g];
        int k = gidx_[g], m = (int)s.widths.size();
        if (a == 0) {   // narrower: lower index
            int64_t target = s.widths[k] - d;
            int c = k;
            while (c > 0 && s.widths[c] > target) --c;
            return c == k ? -1 : c;
        }
        int64_t target = s.heights[k] - d;   // lower: higher index
        int c = k;
        while (c + 1 < m && s.heights[c] > target) ++c;
        return c == k || s.heights[c] >= s.heights[k] ? -1 : c;
    }

    // Relation of pair (i, j) on axis a (direction by key order) impossible
    // outright for the current shapes.
    bool impossible(int i, int j, int a) const {
        int u = rk_[a][i] < rk_[a][j] ? i : j, v = u == i ? j : i;   // u before v
        int lowE = a == 0 ? lab3::EDGE_LEFT : lab3::EDGE_BOTTOM, highE = a == 0 ? lab3::EDGE_RIGHT : lab3::EDGE_TOP;
        bool pu = P_.kind[u] == K_PRE, pv = P_.kind[v] == K_PRE;
        if (!pv && (P_.bmask[v] & lowE)) return true;
        if (!pu && (P_.bmask[u] & highE)) return true;
        int64_t L = a == 0 ? P_.W : P_.H;
        auto pos = [&](int b) -> int64_t { return a == 0 ? P_.px[b] : P_.py[b]; };
        auto sz = [&](int b) -> int64_t { return a == 0 ? w_[b] : h_[b]; };
        if (pu && pv) return pos(u) + sz(u) > pos(v);
        if (pv) return sz(u) > pos(v);
        if (pu) return sz(v) > L - pos(u) - sz(u);
        return false;
    }

    void bounds(int a, int i, int64_t &lo, int64_t &hi) const {
        int64_t L = a == 0 ? P_.W : P_.H, s = a == 0 ? w_[i] : h_[i];
        if (P_.kind[i] == K_PRE) { lo = hi = a == 0 ? P_.px[i] : P_.py[i]; return; }
        lo = 0; hi = L - s;
        int k = P_.bmask[i];
        if (k & (a == 0 ? lab3::EDGE_LEFT : lab3::EDGE_BOTTOM)) hi = 0;
        if (k & (a == 0 ? lab3::EDGE_RIGHT : lab3::EDGE_TOP)) lo = L - s;
    }

    // ASAP (forward, raw: not clamped at hi) and ALAP (backward, raw: not
    // clamped at lo) on both axes. excess = sum over nodes of ASAP - hi > 0;
    // tns = sum over constraint edges i -> j of (ASAP_i + s_i - ALAP_j) > 0,
    // which sees parallel critical paths. With keep, fills e_/exc_/pred_.
    Eval evaluate(bool keep) {
        const int n = n_;
        Eval ev;
        thread_local std::vector<int64_t> e, ex, la, lo, hi;
        thread_local std::vector<int> pr;
        for (int a = 0; a < 2; ++a) {
            e.assign(n, 0); ex.assign(n, 0); pr.assign(n, -1); la.assign(n, 0); lo.assign(n, 0); hi.assign(n, 0);
            const std::vector<int> &ord = ord_[a];
            const std::vector<int> &sz = a == 0 ? w_ : h_;
            for (int r = 0; r < n; ++r) {
                int j = ord[r];
                bounds(a, j, lo[j], hi[j]);
                int64_t best = lo[j];
                int bp = -1;
                const uint8_t *row = &ax_[(size_t)j * n];
                for (int q = 0; q < r; ++q) {
                    int i = ord[q];
                    if (row[i] != a) continue;
                    int64_t v = e[i] + sz[i];
                    if (v > best) { best = v; bp = i; }
                }
                if (best > hi[j]) { ex[j] = best - hi[j]; ev.excess += best - hi[j]; }
                e[j] = best; pr[j] = bp;
            }
            if (ev.excess > 0 || keep) {
                for (int r = n - 1; r >= 0; --r) {
                    int i = ord[r];
                    int64_t best = hi[i];
                    const uint8_t *row = &ax_[(size_t)i * n];
                    const int64_t ei = e[i] + sz[i];
                    for (int q = r + 1; q < n; ++q) {
                        int j = ord[q];
                        if (row[j] != a) continue;
                        if (ei > la[j]) ev.tns += ei - la[j];
                        int64_t v = la[j] - sz[i];
                        if (v < best) best = v;
                    }
                    la[i] = best;
                }
            }
            if (keep) { e_[a] = e; exc_[a] = ex; pred_[a] = pr; la_[a] = la; }
        }
        return ev;
    }

    // Move u out of its band on axis a (the blocks related to u on a whose
    // extent on the other axis b overlaps u's current ASAP extent): every
    // band relation becomes a b-relation. mode 0 keeps u's key on b, mode 1
    // moves u after all band blocks in the b order (above / right), mode 2
    // before them. Returns false if a resulting relation is impossible.
    bool relocate(int u, int a, int mode) {
        const int n = n_, b = 1 - a;
        int lowE = b == 0 ? lab3::EDGE_LEFT : lab3::EDGE_BOTTOM, highE = b == 0 ? lab3::EDGE_RIGHT : lab3::EDGE_TOP;
        if (mode && (P_.bmask[u] & (lowE | highE))) return false;
        const std::vector<int> &sz = b == 0 ? w_ : h_;
        const std::vector<int64_t> &p = e_[b];
        std::vector<int> band;
        for (int k = 0; k < n; ++k) {
            if (k == u || ax_[(size_t)u * n + k] != a) continue;
            if (std::min(p[u] + sz[u], p[k] + sz[k]) > std::max(p[u], p[k])) band.push_back(k);
        }
        if (band.empty()) return false;
        if (mode) {
            int tgt = mode == 1 ? -1 : n;
            for (int k : band) tgt = mode == 1 ? std::max(tgt, rk_[b][k]) : std::min(tgt, rk_[b][k]);
            std::vector<int> &ord = ord_[b];
            int from = rk_[b][u];
            ord.erase(ord.begin() + from);
            if (mode == 1) { if (tgt > from) --tgt; ord.insert(ord.begin() + tgt + 1, u); }
            else { if (tgt > from) --tgt; ord.insert(ord.begin() + tgt, u); }
            for (int r = 0; r < n; ++r) rk_[b][ord[r]] = r;
        }
        for (int k : band) {
            ax_[(size_t)u * n + k] = ax_[(size_t)k * n + u] = (uint8_t)b;
        }
        for (int k = 0; k < n; ++k)
            if (k != u && ax_[(size_t)u * n + k] == b && impossible(u, k, b)) return false;
        return true;
    }
};

inline bool legalize(const Problem &P, const std::vector<double> &cx, const std::vector<double> &cy,
                     std::vector<int> &gidx, std::vector<lab3::Rect> &out, mcf::ConstraintGraphs *g,
                     LegalizeStats &st, const LegalizeOpts &o = LegalizeOpts(),
                     const std::vector<double> *leafw = nullptr, const std::vector<double> *leafh = nullptr) {
    Legalizer L(P);
    L.set_leaf(leafw, leafh);
    return L.run(cx, cy, gidx, out, g, st, o);
}

}  // namespace fp

#endif  // FP_LEGALIZE_H
