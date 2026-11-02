// oracle.h -- exact insertion oracle for one block of a sequence pair
// (docs/research/adhd-solver-round2.md, shortlist item 1).
//
// Given a state S (legal under the solver's rules) and a block b, remove b
// from (Γ+, Γ-) and compute, for every insertion cell (i, j) -- b inserted in
// front of reduced position i of Γ+ and reduced position j of Γ- (i, j in
// [0, m], m = n - 1) -- the exact window that b's lower-left corner and far
// edge must stay in:
//
//   left  of b : reduced Γ+ < i and Γ- < j   -> xlo = max(0, max ASAP x-end)
//   right of b : reduced Γ+ >= i and Γ- >= j -> xhi = min(W, min ALAP x-start)
//   below b    : reduced Γ+ >= i and Γ- < j  -> ylo = max(0, max ASAP y-end)
//   above b    : reduced Γ+ < i and Γ- >= j  -> yhi = min(H, min ALAP y-start)
//
// ASAP / ALAP are longest paths in the reduced SP's constraint graphs with
// preplaced blocks pinned (ASAP = ALAP = px) and the outline as source/sink,
// exactly as Packer computes them. Each SP relation is transitive, so the
// blocks that b's position can push are exactly its right (above) set, and
// those that push b are exactly its left (below) set. Hence, with b of size
// (w, h), Packer on the new SP is feasible (no overflow, no pin violation)
// iff  xhi - xlo >= w  and  yhi - ylo >= h  (pinned b: xlo <= px and
// px + w <= xhi, same in y), and Packer puts b at (xlo, ylo) before flush.
//
// "Legal" = what the solver accepts: Packer.feasible() && invariant_ok().
// The invariant part (boundary source/sink rule) is a per-cell mask: b's own
// edges need the outward quadrant empty, and every other boundary block c
// forbids the corner rectangle of cells that put b on c's outward side.
// Under the invariant the RIGHT/TOP flush cannot create an overlap.
//
// Shapes: fixed, preplaced and MIB members keep their current (w, h) (an MIB
// shape change would move the other members). A non-MIB soft block may take
// any legal width with its minimum height (the solver's shape convention).
//
// Cost: O(n log n) sweeps + O(n^2) cell grids (four int arrays of (n)^2).
#ifndef FP_ORACLE_H
#define FP_ORACLE_H

#include <algorithm>
#include <climits>
#include <vector>

#include "gp.h"

namespace fp {

struct InsertionOracle {
    int n = 0, m = 0, b = -1;
    bool reduced_ok = false;       // reduced SP packs feasibly (else no cell is legal)
    bool free_shape = false;       // non-MIB soft: shape chosen per cell
    int ib = 0, jb = 0;            // original cell of b
    std::vector<int> rp, rm;       // reduced Γ+ / Γ- (block ids)
    std::vector<int> qp, qm;       // reduced positions per block (-1 for b)
    std::vector<int> xend, yend, alx, aly;   // per block: ASAP end, ALAP start
    // Per cell (i * (m + 1) + j):
    std::vector<int> xlo, xhi, ylo, yhi;
    std::vector<unsigned char> inv;          // invariant mask: 1 = allowed
    // Soft shape table (prefix minimum of heights by width).
    std::vector<int> sw, shmin;              // widths ascending, prefix-min heights

    int cells() const { return (m + 1) * (m + 1); }
    int idx(int i, int j) const { return i * (m + 1) + j; }

    // Geometric fit of a given (w, h) at cell c (no invariant).
    bool fits(const Problem &P, int c, int w, int h) const {
        if (!reduced_ok) return false;
        if (P.kind[b] == K_PRE)
            return xlo[c] <= P.px[b] && P.px[b] + w <= xhi[c] && ylo[c] <= P.py[b] && P.py[b] + h <= yhi[c];
        return xhi[c] - xlo[c] >= w && yhi[c] - ylo[c] >= h;
    }
    // Some allowed shape fits geometrically at cell c (no invariant).
    bool geo_ok(const Problem &P, const State &S, int c) const {
        if (!reduced_ok) return false;
        if (!free_shape) return fits(P, c, S.w[b], S.h[b]);
        int wx = xhi[c] - xlo[c], wy = yhi[c] - ylo[c];
        if (wx < sw.front() || wy <= 0) return false;
        int k = int(std::upper_bound(sw.begin(), sw.end(), wx) - sw.begin()) - 1;
        return shmin[k] <= wy;
    }
    bool legal(const Problem &P, const State &S, int c) const { return inv[c] && geo_ok(P, S, c); }

    // Shape index (into the group's width list) for a legal soft cell:
    // the current width if it fits, else the fitting width nearest to it.
    int pick_shape(const Problem &P, const State &S, int c) const {
        const ShapeRange &g = P.gshape[P.group[b]];
        int wx = xhi[c] - xlo[c], wy = yhi[c] - ylo[c];
        int best = -1, bd = INT_MAX;
        for (int k = 0; k < (int)g.widths.size(); ++k) {
            if (g.widths[k] > wx) break;
            if (g.heights[k] > wy) continue;
            int d = std::abs(g.widths[k] - S.w[b]);
            if (d < bd) { bd = d; best = k; }
        }
        return best;
    }

    // New SP with b inserted at cell (i, j) (positions not synced).
    void insert(int i, int j, std::vector<int> &gp, std::vector<int> &gm) const {
        gp.assign(rp.begin(), rp.begin() + i);
        gp.push_back(b);
        gp.insert(gp.end(), rp.begin() + i, rp.end());
        gm.assign(rm.begin(), rm.begin() + j);
        gm.push_back(b);
        gm.insert(gm.end(), rm.begin() + j, rm.end());
    }

    void build(const Problem &P, const State &S, int blk) {
        n = P.n;
        m = n - 1;
        b = blk;
        const int W = P.W, H = P.H;
        rp.clear(); rm.clear();
        for (int v : S.gp) if (v != b) rp.push_back(v);
        for (int v : S.gm) if (v != b) rm.push_back(v);
        ib = S.pp[b]; jb = S.pm[b];
        qp.assign(n, -1); qm.assign(n, -1);
        for (int k = 0; k < m; ++k) { qp[rp[k]] = k; qm[rm[k]] = k; }
        xend.assign(n, 0); yend.assign(n, 0); alx.assign(n, 0); aly.assign(n, 0);

        // Fenwick max over prefixes of reduced Γ- (1-based).
        std::vector<int> bit(m + 1);
        auto qmax = [&](int k) { int r = 0; for (; k > 0; k -= k & -k) r = std::max(r, bit[k]); return r; };
        auto umax = [&](int k, int v) { for (; k <= m; k += k & -k) bit[k] = std::max(bit[k], v); };
        // Suffix-min over Γ- positions via a Fenwick on reversed index.
        auto qmin = [&](int k) { int r = INT_MAX; for (; k > 0; k -= k & -k) r = std::min(r, bit[k]); return r; };
        auto umin = [&](int k, int v) { for (; k <= m; k += k & -k) bit[k] = std::min(bit[k], v); };

        reduced_ok = true;
        // ASAP x: predecessors = before in Γ+ and Γ-.
        std::fill(bit.begin(), bit.end(), 0);
        for (int k = 0; k < m; ++k) {
            int a = rp[k], q = qmax(qm[a]), xa = q;
            if (P.kind[a] == K_PRE) { xa = P.px[a]; if (q > xa) reduced_ok = false; }
            xend[a] = xa + S.w[a];
            if (xend[a] > W) reduced_ok = false;
            umax(qm[a] + 1, xend[a]);
        }
        // ASAP y: predecessors = after in Γ+, before in Γ-.
        std::fill(bit.begin(), bit.end(), 0);
        for (int k = m - 1; k >= 0; --k) {
            int a = rp[k], q = qmax(qm[a]), ya = q;
            if (P.kind[a] == K_PRE) { ya = P.py[a]; if (q > ya) reduced_ok = false; }
            yend[a] = ya + S.h[a];
            if (yend[a] > H) reduced_ok = false;
            umax(qm[a] + 1, yend[a]);
        }
        // ALAP x: successors = after in Γ+ and Γ- (reverse index m - qm).
        std::fill(bit.begin(), bit.end(), INT_MAX);
        for (int k = m - 1; k >= 0; --k) {
            int a = rp[k], q = std::min(W, qmin(m - qm[a] - 1));
            alx[a] = (P.kind[a] == K_PRE) ? P.px[a] : q - S.w[a];
            umin(m - qm[a], alx[a]);
        }
        // ALAP y: successors (above) = before in Γ+, after in Γ-.
        std::fill(bit.begin(), bit.end(), INT_MAX);
        for (int k = 0; k < m; ++k) {
            int a = rp[k], q = std::min(H, qmin(m - qm[a] - 1));
            aly[a] = (P.kind[a] == K_PRE) ? P.py[a] : q - S.h[a];
            umin(m - qm[a], aly[a]);
        }

        // Cell grids. Forward in i: xlo (Γ- < j), yhi (Γ- >= j).
        const int M1 = m + 1, NC = M1 * M1;
        xlo.resize(NC); xhi.resize(NC); ylo.resize(NC); yhi.resize(NC); inv.resize(NC);
        for (int j = 0; j <= m; ++j) { xlo[j] = 0; yhi[j] = H; }
        for (int i = 1; i <= m; ++i) {
            const int a = rp[i - 1], t = qm[a], ex = xend[a], ay = aly[a];
            const int *pxl = &xlo[(i - 1) * M1], *pyh = &yhi[(i - 1) * M1];
            int *cxl = &xlo[i * M1], *cyh = &yhi[i * M1];
            for (int j = 0; j <= t; ++j) { cxl[j] = pxl[j]; cyh[j] = std::min(pyh[j], ay); }
            for (int j = t + 1; j <= m; ++j) { cxl[j] = std::max(pxl[j], ex); cyh[j] = pyh[j]; }
        }
        // Backward in i: xhi (Γ- >= j), ylo (Γ- < j).
        for (int j = 0; j <= m; ++j) { xhi[m * M1 + j] = W; ylo[m * M1 + j] = 0; }
        for (int i = m - 1; i >= 0; --i) {
            const int a = rp[i], t = qm[a], ax = alx[a], ey = yend[a];
            const int *pxh = &xhi[(i + 1) * M1], *pyl = &ylo[(i + 1) * M1];
            int *cxh = &xhi[i * M1], *cyl = &ylo[i * M1];
            for (int j = 0; j <= t; ++j) { cxh[j] = std::min(pxh[j], ax); cyl[j] = pyl[j]; }
            for (int j = t + 1; j <= m; ++j) { cxh[j] = pxh[j]; cyl[j] = std::max(pyl[j], ey); }
        }

        // Invariant mask. Per i, the allowed j form an interval [jl, jh]:
        //   b LEFT   : no a with qp < i, qm < j       -> j <= min qm over qp < i
        //   b TOP    : no a with qp < i, qm >= j      -> j >  max qm over qp < i
        //   b RIGHT  : no a with qp >= i, qm >= j     -> j >  max qm over qp >= i
        //   b BOTTOM : no a with qp >= i, qm < j      -> j <= min qm over qp >= i
        //   c LEFT   forbids i <= qp[c], j <= qm[c]   -> j >  max qm[c] over LEFT c, qp >= i
        //   c TOP    forbids i <= qp[c], j >  qm[c]   -> j <= min qm[c] over TOP c, qp >= i
        //   c RIGHT  forbids i >  qp[c], j >  qm[c]   -> j <= min qm[c] over RIGHT c, qp < i
        //   c BOTTOM forbids i >  qp[c], j <= qm[c]   -> j >  max qm[c] over BOTTOM c, qp < i
        const int kb = P.bmask[b];
        std::vector<int> lo(M1, 0), hi(M1, m);
        {
            int pmin = m, pmax = -1, cr = m, cb = -1;   // prefix over qp < i
            for (int i = 0; i <= m; ++i) {
                if (i > 0) {
                    int a = rp[i - 1], k = P.bmask[a];
                    pmin = std::min(pmin, qm[a]); pmax = std::max(pmax, qm[a]);
                    if (k & lab3::EDGE_RIGHT) cr = std::min(cr, qm[a]);
                    if (k & lab3::EDGE_BOTTOM) cb = std::max(cb, qm[a]);
                }
                if (kb & lab3::EDGE_LEFT) hi[i] = std::min(hi[i], pmin);
                if (kb & lab3::EDGE_TOP) lo[i] = std::max(lo[i], pmax + 1);
                hi[i] = std::min(hi[i], cr);
                lo[i] = std::max(lo[i], cb + 1);
            }
            int smin = m, smax = -1, cl = -1, ct = m;     // suffix over qp >= i
            for (int i = m; i >= 0; --i) {
                if (i < m) {
                    int a = rp[i], k = P.bmask[a];
                    smin = std::min(smin, qm[a]); smax = std::max(smax, qm[a]);
                    if (k & lab3::EDGE_LEFT) cl = std::max(cl, qm[a]);
                    if (k & lab3::EDGE_TOP) ct = std::min(ct, qm[a]);
                }
                if (kb & lab3::EDGE_RIGHT) lo[i] = std::max(lo[i], smax + 1);
                if (kb & lab3::EDGE_BOTTOM) hi[i] = std::min(hi[i], smin);
                lo[i] = std::max(lo[i], cl + 1);
                hi[i] = std::min(hi[i], ct);
            }
        }
        for (int i = 0; i <= m; ++i) {
            unsigned char *r = &inv[i * M1];
            for (int j = 0; j <= m; ++j) r[j] = (j >= lo[i] && j <= hi[i]);
        }

        // Shape table.
        free_shape = P.kind[b] == K_SOFT && P.gmem[P.group[b]].size() == 1;
        if (free_shape) {
            const ShapeRange &g = P.gshape[P.group[b]];
            sw = g.widths;
            shmin.resize(sw.size());
            int mn = INT_MAX;
            for (size_t k = 0; k < sw.size(); ++k) { mn = std::min(mn, g.heights[k]); shmin[k] = mn; }
        }
    }
};

}  // namespace fp

#endif  // FP_ORACLE_H
