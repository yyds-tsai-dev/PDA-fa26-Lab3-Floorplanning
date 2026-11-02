// gp.h -- shared problem model, packer, gate mirror and analytical global
// placer of the strong solver, extracted verbatim from fp.cpp (only change:
// everything is `inline` in namespace fp so the header can be included from
// several translation units). Used by test_legalize.cpp / legalize.h now;
// fp.cpp is meant to include it later and drop its own copies of:
//   Rng, K_*, BBConn, BTConn, Problem, sp_rel, build_forbidden, forb_count,
//   build_problem, State, sp_graphs, apply_group, init_dims, g_rej_* / g_t*
//   stats globals, invariant_ok, Packer, GateInfo, gate_check, to_rects,
//   write_output, sp_from_rects, random_initial_sp, read_placement,
//   initial_from_placement, GlobalPlacer, sp_from_gp, shape_descent.
#ifndef FP_GP_H
#define FP_GP_H

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "extract.h"
#include "lab3_io.h"
#include "mcf/wl_lp.h"
#include "shape.h"

namespace fp {

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

inline void build_forbidden(Problem &P) {
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

inline int64_t forb_count(const Problem &P, const struct State &S);

inline bool build_problem(const std::string &path, Problem &P, std::string &err) {
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

    void sync_pos() {
        for (int i = 0; i < (int)gp.size(); ++i) pp[gp[i]] = i;
        for (int i = 0; i < (int)gm.size(); ++i) pm[gm[i]] = i;
    }
};

inline int64_t forb_count(const Problem &P, const State &S) {
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
inline void sp_graphs(const std::vector<int> &gp, const std::vector<int> &pm, mcf::ConstraintGraphs &g) {
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

inline void apply_group(const Problem &P, State &S, int g) {
    int w = P.gshape[g].widths[S.gidx[g]], h = P.gshape[g].heights[S.gidx[g]];
    for (int b : P.gmem[g]) { S.w[b] = w; S.h[b] = h; }
}

inline void init_dims(const Problem &P, State &S) {
    S.w.assign(P.n, 0);
    S.h.assign(P.n, 0);
    for (int i = 0; i < P.n; ++i) if (P.kind[i] != K_SOFT) { S.w[i] = P.fw[i]; S.h[i] = P.fh[i]; }
    for (int g = 0; g < (int)P.gmem.size(); ++g) apply_group(P, S, g);
}

// Boundary source/sink invariant + preplaced-pair geometric consistency.
inline int64_t g_rej_boundary = 0, g_rej_pre = 0;  // rejection causes (stats only)
inline bool invariant_ok(const Problem &P, const State &S) {
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
inline GateInfo gate_check(const Problem &P, const std::vector<lab3::Rect> &r) {
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

inline std::vector<lab3::Rect> to_rects(const Packer &K, const State &S) {
    std::vector<lab3::Rect> r(K.n);
    for (int i = 0; i < K.n; ++i) r[i] = lab3::Rect{K.x[i], K.y[i], S.w[i], S.h[i]};
    return r;
}

inline bool write_output(const Problem &P, const std::vector<lab3::Rect> &r, const std::string &path) {
    std::string tmp = path + ".tmp";
    {
        std::ofstream f(tmp);
        if (!f) return false;
        f << "WIRELENGTH " << lab3::format_half(lab3::wirelength2(P.inst, r)) << "\n";
        f << "BLOCKS " << P.n << "\n";
        for (int i = 0; i < P.n; ++i)
            f << P.inst.blocks[i].name << ' ' << r[i].x << ' ' << r[i].y << ' ' << r[i].w << ' ' << r[i].h << "\n";
        if (!f) return false;
    }
    return std::rename(tmp.c_str(), path.c_str()) == 0;
}

// ------------------------------------------------------------------ initial SP

// SP from a placement (all blocks). Soft (w,h) from the placement set the
// group width index (nearest legal width; MIB groups take the first member).
inline bool sp_from_rects(const Problem &P, State &S, const std::vector<fp::XRect> &r, Rng &rng, std::string &err) {
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
inline bool random_initial_sp(const Problem &P, State &S, Rng &rng, std::string &err) {
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

inline bool read_placement(const Problem &P, const std::string &path, std::vector<lab3::Rect> &r, std::string &err) {
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

inline bool initial_from_placement(const Problem &P, State &S, const std::string &path, Rng &rng, std::string &err) {
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

    // ---- obstacle-aware bisection (obstacle = true; default path unchanged)
    //
    // Same slicing recursion, but every cut is chosen from candidates that
    // cross no preplaced block inside the region: for each split k of the
    // blocks (sorted by coordinate on the cut axis) the cut that gives side A
    // exactly its share of the free area, plus every preplaced edge in the
    // region (with the k that best balances utilization). Candidates on both
    // axes are scored by
    //   capacity mismatch + displacement of blocks out of their side
    //   + misfit (a side narrower than the narrowest legal shape of a block in
    //     it, or a region lower than one) + a small penalty for the short axis.
    // Regions still holding a preplaced block keep being cut (k may be 0 or m:
    // a side without blocks is dropped), so leaves never contain a preplaced
    // block when a guillotine separation exists. Boundary blocks are forced
    // to the side touching their edge, so they end in a leaf on that edge
    // (corner blocks in the corner leaf).
    bool obstacle = false;
    double ob_wcap = 4.0, ob_wdisp = 1.0, ob_wfit = 20.0, ob_waxis = 0.05, ob_wutil = 1.0, ob_umax = 0.97;
    std::vector<double> minw, minh;

    void init_min_dims() {
        minw.assign(n, 0); minh.assign(n, 0);
        for (int i = 0; i < n; ++i) {
            if (P.kind[i] == K_SOFT) {
                const fp::ShapeRange &g = P.gshape[P.group[i]];
                minw[i] = g.widths.front();
                minh[i] = *std::min_element(g.heights.begin(), g.heights.end());
            } else { minw[i] = P.fw[i]; minh[i] = P.fh[i]; }
        }
    }

    // Free area of [x0,x1]x[y0,y1] minus the given preplaced blocks.
    double free_area_in(const std::vector<int> &pres, double x0, double x1, double y0, double y1) const {
        double a = (x1 - x0) * (y1 - y0);
        for (int p : pres) {
            double ox = std::min(x1, (double)P.px[p] + P.fw[p]) - std::max(x0, (double)P.px[p]);
            double oy = std::min(y1, (double)P.py[p] + P.fh[p]) - std::max(y0, (double)P.py[p]);
            if (ox > 0 && oy > 0) a -= ox * oy;
        }
        return std::max(a, 0.0);
    }

    // Slicing tree of the last obstacle bisection: leaf (blk >= 0), cut node
    // (vert, a = low side, b = high side) or opaque (blk = -2: fallback
    // bisection, not sized). free = the region holds no preplaced block.
    struct ObNode { int blk = -1, vert = 0, a = -1, b = -1; bool free = false; double x0, x1, y0, y1; };
    std::vector<ObNode> obn;
    int ob_root = -1;
    int new_node(double x0, double x1, double y0, double y1, bool fr) {
        ObNode nd; nd.x0 = x0; nd.x1 = x1; nd.y0 = y0; nd.y1 = y1; nd.free = fr;
        obn.push_back(nd);
        return (int)obn.size() - 1;
    }

    int bisect_ob(std::vector<int> &ids, double x0, double x1, double y0, double y1, std::vector<int> &sp,
                  std::vector<int> &sm) {
        if (ids.empty()) return -1;
        std::vector<int> pres;
        for (int p : P.pre) {
            double ox = std::min(x1, (double)P.px[p] + P.fw[p]) - std::max(x0, (double)P.px[p]);
            double oy = std::min(y1, (double)P.py[p] + P.fh[p]) - std::max(y0, (double)P.py[p]);
            if (ox > 1e-9 && oy > 1e-9) pres.push_back(p);
        }
        const int m = (int)ids.size();
        if (m == 1 && pres.empty()) {
            int b = ids[0];
            tx[b] = (x0 + x1) / 2; ty[b] = (y0 + y1) / 2; leafw[b] = x1 - x0; leafh[b] = y1 - y0;
            sp.push_back(b); sm.push_back(b);
            int nd = new_node(x0, x1, y0, y1, true);
            obn[nd].blk = b;
            return nd;
        }
        if (m == 1) {   // one block beside preplaced blocks: largest fitting maximal strip
            int b = ids[0];
            double sx0 = x0, sx1 = x1, sy0 = y0, sy1 = y1, bestA = -1, bfit = -1;
            double L = x1, R = x0, Bo = y1, T = y0;   // extent of the preplaced blocks in the region
            for (int p : pres) {
                L = std::min(L, std::max(x0, (double)P.px[p])); R = std::max(R, std::min(x1, (double)P.px[p] + P.fw[p]));
                Bo = std::min(Bo, std::max(y0, (double)P.py[p])); T = std::max(T, std::min(y1, (double)P.py[p] + P.fh[p]));
            }
            const double cand[4][4] = {{x0, L, y0, y1}, {R, x1, y0, y1}, {x0, x1, y0, Bo}, {x0, x1, T, y1}};
            for (auto &q : cand) {
                double w = q[1] - q[0], h = q[3] - q[2];
                if (w <= 0 || h <= 0) continue;
                double fit = (w >= minw[b] && h >= minh[b] && w * h >= area[b]) ? 1 : 0;
                if (fit > bfit || (fit == bfit && w * h > bestA)) { bfit = fit; bestA = w * h; sx0 = q[0]; sx1 = q[1]; sy0 = q[2]; sy1 = q[3]; }
            }
            tx[b] = (sx0 + sx1) / 2; ty[b] = (sy0 + sy1) / 2; leafw[b] = sx1 - sx0; leafh[b] = sy1 - sy0;
            sp.push_back(b); sm.push_back(b);
            int nd = new_node(sx0, sx1, sy0, sy1, true);
            obn[nd].blk = b;
            return nd;
        }
        double tot = 0;
        for (int b : ids) tot += area[b];
        const double ftot = free_area_in(pres, x0, x1, y0, y1);
        const double wr = x1 - x0, hr = y1 - y0;
        struct Cand { double cost; int vert, k; double cut; };
        Cand best{1e300, -1, -1, 0};
        std::vector<int> ord[2];
        for (int vert = 0; vert < 2; ++vert) {   // vert = 1: vertical cut line (split along x)
            const std::vector<double> &c = vert ? x : y;
            const double lo = vert ? x0 : y0, hi = vert ? x1 : y1, ext = hi - lo, oext = vert ? hr : wr;
            const int lowE = vert ? lab3::EDGE_LEFT : lab3::EDGE_BOTTOM, highE = vert ? lab3::EDGE_RIGHT : lab3::EDGE_TOP;
            const bool atlo = lo <= 0.5, athi = hi >= (vert ? P.W : P.H) - 0.5;
            std::vector<int> &o = ord[vert];
            o = ids;
            auto key = [&](int b) {
                double v = c[b];
                if (atlo && (P.bmask[b] & lowE)) v -= 1e12;
                if (athi && (P.bmask[b] & highE)) v += 1e12;
                return v;
            };
            std::sort(o.begin(), o.end(), [&](int a, int b) { double ka = key(a), kb = key(b); return ka != kb ? ka < kb : a < b; });
            int kmin = 0, kmax = m;
            for (int b : o) { if (atlo && (P.bmask[b] & lowE)) ++kmin; if (athi && (P.bmask[b] & highE)) --kmax; }
            if (kmin > kmax) continue;
            // prefix area, prefix max of min-dim along the axis, suffix max; max other-axis min-dim
            std::vector<double> pa(m + 1, 0), pmx(m + 1, 0), smx(m + 2, 0);
            const std::vector<double> &md = vert ? minw : minh, &od = vert ? minh : minw;
            double omax = 0;
            for (int i = 0; i < m; ++i) { pa[i + 1] = pa[i] + area[o[i]]; pmx[i + 1] = std::max(pmx[i], md[o[i]]); omax = std::max(omax, od[o[i]]); }
            for (int i = m - 1; i >= 0; --i) smx[i] = std::max(smx[i + 1], md[o[i]]);
            const double omis = std::max(0.0, omax - oext) / std::max(omax, 1.0);
            auto fA = [&](double cut) { return vert ? free_area_in(pres, x0, cut, y0, y1) : free_area_in(pres, x0, x1, y0, cut); };
            auto valid = [&](double cut) {
                if (cut <= lo + 1e-6 || cut >= hi - 1e-6) return false;
                for (int p : pres) {
                    double a = vert ? P.px[p] : P.py[p], b = a + (vert ? P.fw[p] : P.fh[p]);
                    if (cut > a + 1e-6 && cut < b - 1e-6) return false;
                }
                return true;
            };
            auto score = [&](double cut, int k, double FA) {
                double FB = ftot - FA, acc = pa[k], rest = tot - acc;
                if ((acc > 0 && FA <= 1e-9) || (rest > 0 && FB <= 1e-9)) return 1e300;
                double u = tot / std::max(ftot, 1e-9);
                double cap = (std::fabs(acc - u * FA) + std::fabs(rest - u * FB)) / tot;
                double disp = 0;
                for (int i = 0; i < k; ++i) { double kk = c[o[i]]; if (kk > cut) disp += area[o[i]] * (kk - cut); }
                for (int i = k; i < m; ++i) { double kk = c[o[i]]; if (kk < cut) disp += area[o[i]] * (cut - kk); }
                disp /= tot * ext;
                double mis = 0;
                if (k > 0) mis += std::max(0.0, pmx[k] - (cut - lo)) / pmx[k];
                if (k < m) mis += std::max(0.0, smx[k] - (hi - cut)) / smx[k];
                mis += omis;
                if (k == 0) mis += FA / ftot;   // free area stranded without blocks
                if (k == m) mis += FB / ftot;
                if (k > 0) mis += ob_wutil * std::max(0.0, acc / (ob_umax * FA) - 1);
                if (k < m) mis += ob_wutil * std::max(0.0, rest / (ob_umax * FB) - 1);
                bool pref = vert ? wr >= hr : hr > wr;
                return ob_wcap * cap + ob_wdisp * disp + ob_wfit * mis + (pref ? 0.0 : ob_waxis);
            };
            auto consider = [&](double cut, int k, double FA) {
                double s = score(cut, k, FA);
                if (s < best.cost) best = {s, vert, k, cut};
            };
            // exact-capacity cuts per split k
            for (int k = std::max(kmin, 1); k <= std::min(kmax, m - 1); ++k) {
                double frac = pa[k] / tot, a = lo, b = hi;
                if (pres.empty()) {
                    double cut = lo + frac * ext;
                    if (valid(cut)) consider(cut, k, frac * ftot);
                    continue;
                }
                for (int it = 0; it < 40; ++it) {
                    double mid = (a + b) / 2;
                    if (fA(mid) < frac * ftot) a = mid; else b = mid;
                }
                double cut = (a + b) / 2;
                if (valid(cut)) consider(cut, k, fA(cut));
            }
            // fit cuts: side A exactly as wide as its widest block, or side B
            for (int k = std::max(kmin, 1); k <= std::min(kmax, m - 1); ++k) {
                for (double cut : {lo + pmx[k], hi - smx[k]})
                    if (valid(cut)) consider(cut, k, fA(cut));
            }
            // preplaced edges, k balancing utilization
            for (int p : pres) {
                double e0 = vert ? P.px[p] : P.py[p], e1 = e0 + (vert ? P.fw[p] : P.fh[p]);
                for (double cut : {e0, e1}) {
                    if (!valid(cut)) continue;
                    double FA = fA(cut);
                    double want = tot * FA / std::max(ftot, 1e-9);
                    int k = (int)(std::lower_bound(pa.begin(), pa.end(), want) - pa.begin());
                    for (int kk : {k - 1, k}) {
                        if (kk < kmin || kk > kmax) continue;
                        consider(cut, kk, FA);
                    }
                }
            }
        }
        if (best.vert < 0) {   // no guillotine cut avoids the preplaced blocks: original rule
            int nd = new_node(x0, x1, y0, y1, false);
            obn[nd].blk = -2;
            if (m == 1) {
                int b = ids[0];
                tx[b] = (x0 + x1) / 2; ty[b] = (y0 + y1) / 2; leafw[b] = x1 - x0; leafh[b] = y1 - y0;
                sp.push_back(b); sm.push_back(b);
                return nd;
            }
            bisect(ids, x0, x1, y0, y1, sp, sm);
            return nd;
        }
        const std::vector<int> &o = ord[best.vert];
        std::vector<int> A(o.begin(), o.begin() + best.k), B(o.begin() + best.k, o.end());
        std::vector<int> pa, ma, pb, mb;
        int nd = new_node(x0, x1, y0, y1, pres.empty()), na, nb;
        if (best.vert) {
            na = bisect_ob(A, x0, best.cut, y0, y1, pa, ma);
            nb = bisect_ob(B, best.cut, x1, y0, y1, pb, mb);
            sp.insert(sp.end(), pa.begin(), pa.end()); sp.insert(sp.end(), pb.begin(), pb.end());
            sm.insert(sm.end(), ma.begin(), ma.end()); sm.insert(sm.end(), mb.begin(), mb.end());
        } else {
            na = bisect_ob(A, x0, x1, y0, best.cut, pa, ma);
            nb = bisect_ob(B, x0, x1, best.cut, y1, pb, mb);
            sp.insert(sp.end(), pb.begin(), pb.end()); sp.insert(sp.end(), pa.begin(), pa.end());
            sm.insert(sm.end(), ma.begin(), ma.end()); sm.insert(sm.end(), mb.begin(), mb.end());
        }
        if (na < 0 || nb < 0) return na < 0 ? nb : na;   // one side empty: no cut node
        obn[nd].vert = best.vert; obn[nd].a = na; obn[nd].b = nb;
        return nd;
    }

    // ---- shape-curve sizing of the pres-free subtrees (obstacle mode)
    //
    // Every maximal pres-free region of the slicing tree is sized bottom-up
    // (Stockmeyer shape curves: soft leaves all legal shapes, fixed leaves
    // their one shape; a vertical cut adds widths, a horizontal one heights;
    // points beyond the region are dropped, curves thinned to 256 points),
    // then realized top-down inside the region: the root point with the
    // most relative slack, child boxes split in proportion to their widths
    // (heights), each block centred in its box (boundary blocks against
    // their edge). If no point fits, the least-overflowing one is used.
    // MIB groups: a second pass forces every member to the group's median
    // chosen width. Results: sx, sy (centres), sgidx (width index per group).
    struct SPt { int w, h, a, b; };
    std::vector<std::vector<SPt>> curv;
    std::vector<double> sx, sy;
    std::vector<int> sgidx, sbw, forced;
    int ob_nfree = 0, ob_nmiss = 0;

    static void pareto(std::vector<SPt> &c) {   // sorted by w asc; keep strictly decreasing h
        std::vector<SPt> o;
        for (const SPt &p : c) {
            if (!o.empty() && p.w == o.back().w) { if (p.h < o.back().h) o.back() = p; continue; }
            if (o.empty() || p.h < o.back().h) o.push_back(p);
        }
        c.swap(o);
    }
    void prune(std::vector<SPt> &c, double Rw, double Rh) {
        std::vector<SPt> o;
        int bi = -1;
        double bo = 1e300;
        for (int i = 0; i < (int)c.size(); ++i) {
            if (c[i].w <= Rw + 1e-9 && c[i].h <= Rh + 1e-9) o.push_back(c[i]);
            double ov = std::max(c[i].w / Rw, c[i].h / Rh);
            if (ov < bo) { bo = ov; bi = i; }
        }
        if (o.empty() && bi >= 0) o.push_back(c[bi]);
        if (o.size() > 256) {
            std::vector<SPt> t;
            for (int k = 0; k < 256; ++k) t.push_back(o[(size_t)k * (o.size() - 1) / 255]);
            o.swap(t);
        }
        c.swap(o);
    }
    void size_node(int nd, double Rw, double Rh) {
        ObNode &N = obn[nd];
        std::vector<SPt> &c = curv[nd];
        c.clear();
        if (N.blk >= 0) {
            int b = N.blk;
            if (P.kind[b] == K_SOFT) {
                int g = P.group[b];
                const fp::ShapeRange &sr = P.gshape[g];
                if (forced[g] >= 0) c.push_back({sr.widths[forced[g]], sr.heights[forced[g]], forced[g], -1});
                else for (int i = 0; i < (int)sr.widths.size(); ++i) c.push_back({sr.widths[i], sr.heights[i], i, -1});
            } else {
                c.push_back({P.fw[b], P.fh[b], -1, -1});
            }
            pareto(c);
            prune(c, Rw, Rh);
            return;
        }
        size_node(N.a, Rw, Rh);
        size_node(N.b, Rw, Rh);
        const std::vector<SPt> &A = curv[N.a], &B = curv[N.b];
        size_t i = 0, j = 0;
        if (N.vert) {   // side by side: w = wA + wB, h = max
            while (i < A.size() && j < B.size()) {
                c.push_back({A[i].w + B[j].w, std::max(A[i].h, B[j].h), (int)i, (int)j});
                if (A[i].h > B[j].h) ++i; else if (B[j].h > A[i].h) ++j; else { ++i; ++j; }
            }
        } else {        // stacked: w = max, h = hA + hB
            while (i < A.size() && j < B.size()) {
                c.push_back({std::max(A[i].w, B[j].w), A[i].h + B[j].h, (int)i, (int)j});
                if (A[i].w < B[j].w) ++i; else if (B[j].w < A[i].w) ++j; else { ++i; ++j; }
            }
        }
        std::stable_sort(c.begin(), c.end(), [](const SPt &p, const SPt &q) { return p.w != q.w ? p.w < q.w : p.h < q.h; });
        pareto(c);
        prune(c, Rw, Rh);
    }
    void realize(int nd, int pi, double x0, double y0, double bw, double bh) {
        const ObNode &N = obn[nd];
        const SPt &p = curv[nd][pi];
        if (N.blk >= 0) {
            int b = N.blk, k = P.bmask[b];
            double cx = x0 + bw / 2, cy = y0 + bh / 2;
            if (k & lab3::EDGE_LEFT) cx = p.w / 2.0;
            if (k & lab3::EDGE_RIGHT) cx = P.W - p.w / 2.0;
            if (k & lab3::EDGE_BOTTOM) cy = p.h / 2.0;
            if (k & lab3::EDGE_TOP) cy = P.H - p.h / 2.0;
            sx[b] = cx; sy[b] = cy;
            sbw[b] = P.kind[b] == K_SOFT ? p.a : -1;
            return;
        }
        const SPt &pa = curv[N.a][p.a], &pb = curv[N.b][p.b];
        if (N.vert) {
            double wa = bw * pa.w / std::max(1, pa.w + pb.w);
            realize(N.a, p.a, x0, y0, wa, bh);
            realize(N.b, p.b, x0 + wa, y0, bw - wa, bh);
        } else {
            double ha = bh * pa.h / std::max(1, pa.h + pb.h);
            realize(N.a, p.a, x0, y0, bw, ha);
            realize(N.b, p.b, x0, y0 + ha, bw, bh - ha);
        }
    }
    void size_walk(int nd) {
        if (nd < 0) return;
        const ObNode &N = obn[nd];
        if (N.blk == -2) return;
        if (!N.free) { size_walk(N.a); size_walk(N.b); return; }
        double Rw = N.x1 - N.x0, Rh = N.y1 - N.y0;
        size_node(nd, Rw, Rh);
        const std::vector<SPt> &c = curv[nd];
        int bi = 0;
        double bs = -1e300;
        for (int i = 0; i < (int)c.size(); ++i) {
            double sl = std::min((Rw - c[i].w) / Rw, (Rh - c[i].h) / Rh);
            if (sl > bs) { bs = sl; bi = i; }
        }
        ++ob_nfree;
        if (bs < 0) ++ob_nmiss;
        if (std::getenv("OB_DEBUG")) {
            int cnt = 0; double ar = 0; int nfix = 0;
            std::vector<int> st{nd};
            while (!st.empty()) { int q = st.back(); st.pop_back(); if (obn[q].blk >= 0) { ++cnt; ar += area[obn[q].blk]; nfix += P.kind[obn[q].blk] == K_FIXED; } else if (obn[q].blk == -1) { st.push_back(obn[q].a); st.push_back(obn[q].b); } }
            std::fprintf(stderr, "FREE m=%d fixed=%d R=%.0fx%.0f util=%.3f best=%dx%d slack=%.3f\n", cnt, nfix, Rw, Rh, ar / (Rw * Rh), c[bi].w, c[bi].h, bs);
        }
        realize(nd, bi, N.x0, N.y0, Rw, Rh);
    }
    void size_tree() {
        sx = tx; sy = ty;
        sbw.assign(n, -1);
        forced.assign(P.gmem.size(), -1);
        curv.assign(obn.size(), {});
        ob_nfree = ob_nmiss = 0;
        size_walk(ob_root);
        bool mib = false;
        for (int g = 0; g < (int)P.gmem.size(); ++g) {
            if (P.gmem[g].size() < 2) continue;
            std::vector<int> v;
            for (int b : P.gmem[g]) if (sbw[b] >= 0) v.push_back(sbw[b]);
            if (v.empty()) continue;
            std::sort(v.begin(), v.end());
            forced[g] = v[v.size() / 2];
            mib = true;
        }
        if (mib) { ob_nfree = ob_nmiss = 0; size_walk(ob_root); }
        sgidx.assign(P.gmem.size(), -1);
        for (int g = 0; g < (int)P.gmem.size(); ++g) {
            if (forced[g] >= 0) { sgidx[g] = forced[g]; continue; }
            for (int b : P.gmem[g]) if (sbw[b] >= 0) sgidx[g] = sbw[b];
        }
    }

    void spread() {
        std::vector<int> ids;
        for (int i = 0; i < n; ++i) if (P.kind[i] != K_PRE) ids.push_back(i);
        seqp.clear(); seqm.clear();
        if (obstacle) {
            if ((int)minw.size() != n) init_min_dims();
            obn.clear();
            ob_root = bisect_ob(ids, 0, P.W, 0, P.H, seqp, seqm);
            return;
        }
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
        if (obstacle) size_tree();
    }
};

// SP + soft shapes from a GlobalPlacer result (see above).
inline bool sp_from_gp(const Problem &P, State &S, const GlobalPlacer &G, bool grow, bool geo_prio, std::string &err) {
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

// Soft shapes that fit the bisection leaves (obstacle mode): per shape group
// the width closest (log scale) to the leaf-aspect target among those whose
// (w, h) fits every member's leaf; if none fits, the one with the least
// total overhang.
inline void gp_fit_shapes(const Problem &P, const GlobalPlacer &G, std::vector<int> &gidx) {
    gidx.assign(P.gmem.size(), 0);
    for (int g = 0; g < (int)P.gmem.size(); ++g) {
        double lr = 0, lw = 1e300, lh = 1e300;
        for (int b : P.gmem[g]) {
            lr += std::log(std::max(1e-9, G.leafw[b]) / std::max(1e-9, G.leafh[b]));
            lw = std::min(lw, G.leafw[b]); lh = std::min(lh, G.leafh[b]);
        }
        lr /= P.gmem[g].size();
        const auto &ws = P.gshape[g].widths, &hs = P.gshape[g].heights;
        double target = std::sqrt((double)P.gshape[g].A * std::exp(lr));
        int best = -1;
        double bc = 1e300;
        for (int i = 0; i < (int)ws.size(); ++i) {
            double over = std::max(0.0, ws[i] - lw) + std::max(0.0, hs[i] - lh);
            double c = (over > 0 ? 1e9 + over : 0) + std::fabs(std::log(ws[i] / target));
            if (c < bc) { bc = c; best = i; }
        }
        gidx[g] = std::max(best, 0);
    }
}

// Greedy soft-shape descent on a fixed SP: per shape group try a few width
// steps and keep any that lowers the normalized overflow + violation.
inline void shape_descent(const Problem &P, State &t, Packer &K2) {
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

}  // namespace fp

#endif  // FP_GP_H
