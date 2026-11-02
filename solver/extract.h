// extract.h -- sequence pair from a non-overlapping placement.
//
// SP convention (Murata), used everywhere in the solver:
//   a before b in Γ+ and before b in Γ-  =>  a left-of b   (H relation)
//   a after  b in Γ+ and before b in Γ-  =>  a below b     (V relation)
//
// Every pair of non-overlapping rectangles forces orders in at most two
// sequences (touching counts as disjoint, as in the verifier):
//   a x-before b, y-overlap          : a left-of b  -> Γ+: a<b, Γ-: a<b
//   a y-before b, x-overlap          : a below b    -> Γ+: b<a, Γ-: a<b
//   a lower-left of b (both disjoint): left-of or below, both need Γ-: a<b
//   a upper-left of b (both disjoint): a left-of b or b below a, both need Γ+: a<b
// Any linear extension of each partial order then gives an SP whose every
// pairwise relation is geometrically true for the input placement. Packing
// that SP (longest path, pinned preplaced) therefore never moves a block past
// its input position, and an outline-edge block has nothing on its outer side,
// so the boundary source/sink invariant holds automatically.
#ifndef FP_EXTRACT_H
#define FP_EXTRACT_H

#include <cstdint>
#include <queue>
#include <string>
#include <utility>
#include <vector>

namespace fp {

struct XRect {
    int64_t x, y, w, h;
};

// prio_p / prio_m: tie-break keys for the topological sorts (smaller first).
// Returns false and sets err on an overlapping pair or a cycle.
inline bool extract_sp(const std::vector<XRect> &r, const std::vector<int64_t> &prio_p,
                       const std::vector<int64_t> &prio_m, std::vector<int> &gp, std::vector<int> &gm,
                       std::string &err) {
    const int n = (int)r.size();
    std::vector<std::vector<int>> ap(n), am(n);  // edges a -> b means a before b
    std::vector<int> inp(n, 0), inm(n, 0);
    for (int a = 0; a < n; ++a) {
        for (int b = a + 1; b < n; ++b) {
            const XRect &A = r[a], &B = r[b];
            int xr = 0, yr = 0;  // -1: a before b on that axis, +1: b before a, 0: overlap
            if (A.x + A.w <= B.x) xr = -1;
            else if (B.x + B.w <= A.x) xr = 1;
            if (A.y + A.h <= B.y) yr = -1;
            else if (B.y + B.h <= A.y) yr = 1;
            if (xr == 0 && yr == 0) {
                err = "overlapping rectangles " + std::to_string(a) + " " + std::to_string(b);
                return false;
            }
            // Normalize so that L is x-before-or-overlapping R.
            int L = a, R = b, ly = yr;
            if (xr == 1) { L = b; R = a; ly = -yr; }
            auto addp = [&](int u, int v) { ap[u].push_back(v); ++inp[v]; };
            auto addm = [&](int u, int v) { am[u].push_back(v); ++inm[v]; };
            if (xr != 0 && ly == 0) {            // L left-of R
                addp(L, R); addm(L, R);
            } else if (xr == 0) {                // vertical: yr decides
                int lo = (yr == -1) ? a : b, hi = (yr == -1) ? b : a;  // lo below hi
                addp(hi, lo); addm(lo, hi);
            } else if (ly == -1) {               // L lower-left of R
                addm(L, R);
            } else {                             // L upper-left of R
                addp(L, R);
            }
        }
    }
    auto topo = [&](std::vector<std::vector<int>> &adj, std::vector<int> &indeg,
                    const std::vector<int64_t> &prio, std::vector<int> &out, const char *nm) {
        typedef std::pair<int64_t, int> PI;
        std::priority_queue<PI, std::vector<PI>, std::greater<PI>> q;
        for (int i = 0; i < n; ++i) if (indeg[i] == 0) q.push(PI(prio[i], i));
        out.clear();
        while (!q.empty()) {
            int u = q.top().second;
            q.pop();
            out.push_back(u);
            for (int v : adj[u]) if (--indeg[v] == 0) q.push(PI(prio[v], v));
        }
        if ((int)out.size() != n) { err = std::string("cycle in ") + nm + " partial order"; return false; }
        return true;
    };
    if (!topo(ap, inp, prio_p, gp, "Gamma+")) return false;
    if (!topo(am, inm, prio_m, gm, "Gamma-")) return false;
    return true;
}

}  // namespace fp

#endif  // FP_EXTRACT_H
