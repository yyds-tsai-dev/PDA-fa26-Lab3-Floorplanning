// wl_lp.h -- exact fixed-topology wirelength LP via min-cost flow
// (header-only, C++17). Formulation: docs/research/lp-polish-results.md.
//
// Per axis, with integer lower-left coordinates x_i and a root r (x_r = 0):
//   min  sum_nets c * |2cx_a - 2cx_b|          (= 2 * WL on this axis)
//   s.t. x_j >= x_i + size_i                   per constraint edge i -> j
//        lo_i <= x_i <= hi_i
// A net between a and b, d = x_a - x_b, delta = size_a - size_b (fixed point
// at doubled coordinate p2: b = r, delta = size_a - p2) contributes
// c*|2d + delta| = c*|d + k1| + c*|d + k2| with k1 = floor(delta/2),
// k2 = delta - k1 (exact for integer d). Every net therefore always has two
// hinges, so the arc set does not depend on block sizes and a size change is
// a cost-only update (warm start).
//
// Dual min-cost circulation (all supplies 0, node n is r):
//   hinge c*|x_a - x_b + k| : a->b cost -k cap c,  b->a cost +k cap c
//   edge  x_j - x_i >= s_i  : i->j cost -s_i cap INF
//   x_i >= lo_i             : r->i cost -lo_i cap INF
//   x_i <= hi_i             : i->r cost +hi_i cap INF
// Optimal x_i = pi_r - pi_i (integral), objective = -flow cost.
//
// Also: H/V constraint graphs from a legal placement, bitset transitive
// reduction, and polish(): both axes on the placement's own topology.
#ifndef MCF_WL_LP_H
#define MCF_WL_LP_H

#include "network_simplex.h"
#include "lab3_io.h"

#include <algorithm>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

namespace mcf {

using Edge = std::pair<int, int>;   // (i, j): x_j >= x_i + size_i

// 2-pin net on one axis. b >= 0: block b. b < 0: fixed point whose doubled
// coordinate is p2 (a terminal at t gives p2 = 2*t).
struct AxisNet {
    int a = 0, b = -1;
    int64_t weight = 1;
    int64_t p2 = 0;
};

struct AxisProblem {
    std::vector<int64_t> size;   // block extent on this axis (w or h)
    std::vector<int64_t> lo, hi; // bounds on the lower-left coordinate
    std::vector<Edge> edges;     // full or transitively reduced, any order
    std::vector<AxisNet> nets;
};

struct AxisResult {
    bool ok = false;             // false: constraints infeasible (positive cycle)
    std::vector<int64_t> x;      // optimal lower-left coordinates
    int64_t obj2 = 0;            // sum c*|2cx_a - 2cx_b| = 2 * WL on this axis
    long long pivots = 0;
    bool warm = false;
};

inline int64_t floor_div2(int64_t v) { return v >= 0 ? v / 2 : -((-v + 1) / 2); }

// Reusable per-axis solver. Nets (and the block count) are structural; edges,
// sizes and bounds can change between solves with warm start:
//   - set_problem() with the same nets: sizes/bounds are cost updates, edges
//     are diffed (a removed edge's arc gets a slack cost making it a no-op,
//     a new edge revives such an arc, retargets a free one or appends one).
//   - the tree, flows and arc arrays are kept; only a changed net list or
//     block count rebuilds from scratch.
class AxisSolver {
public:
    void set_problem(const AxisProblem &p) {
        int n = (int)p.size.size();
        if (n != n_ || !same_nets(p.nets)) build(n, p.nets);
        size_ = p.size; lo_ = p.lo; hi_ = p.hi;
        set_edges(p.edges);
        update_costs();
    }

    AxisResult solve(bool warm = true) {
        AxisResult r;
        if (!warm || !ns_.has_tree()) make_crash();
        NetworkSimplex::Status st = ns_.solve(warm);
        r.pivots = ns_.pivots() - last_pivots_;
        last_pivots_ = ns_.pivots();
        r.warm = ns_.last_solve_was_warm();
        if (st != NetworkSimplex::OPTIMAL) return r;
        int root = n_;
        r.x.resize(n_);
        for (int i = 0; i < n_; ++i) r.x[i] = ns_.potential(root) - ns_.potential(i);
        r.obj2 = objective2(r.x);
        // Primal feasibility + strong duality (exact) certify optimality.
        r.ok = r.obj2 == -ns_.total_cost();
        for (int i = 0; i < n_ && r.ok; ++i) r.ok = lo_[i] <= r.x[i] && r.x[i] <= hi_[i];
        for (const Slot &s : slots_)
            if (s.alive && r.x[s.j] - r.x[s.i] < size_[s.i]) r.ok = false;
        return r;
    }

    // 2*WL on this axis for coordinates x under the current sizes.
    int64_t objective2(const std::vector<int64_t> &x) const {
        int64_t s = 0;
        for (const AxisNet &nt : nets_) {
            int64_t a2 = 2 * x[nt.a] + size_[nt.a];
            int64_t b2 = nt.b >= 0 ? 2 * x[nt.b] + size_[nt.b] : nt.p2;
            s += nt.weight * (a2 > b2 ? a2 - b2 : b2 - a2);
        }
        return s;
    }

    const NetworkSimplex &simplex() const { return ns_; }

private:
    int n_ = -1;
    std::vector<AxisNet> nets_;
    std::vector<int64_t> size_, lo_, hi_;
    NetworkSimplex ns_;
    int net_arc0_ = 0, bound_arc0_ = 0;
    // Constraint-edge arcs: slot -> (i, j, arc id, alive).
    struct Slot { int i, j, arc; bool alive; };
    std::vector<Slot> slots_;
    std::vector<int> slot_of_;          // i*n + j -> slot or -1
    std::vector<int> free_;             // dead slots
    std::vector<char> want_;
    long long last_pivots_ = 0;

    bool same_nets(const std::vector<AxisNet> &nets) const {
        if (nets.size() != nets_.size()) return false;
        for (size_t k = 0; k < nets.size(); ++k) {
            const AxisNet &a = nets[k], &b = nets_[k];
            if (a.a != b.a || a.b != b.b || a.weight != b.weight || a.p2 != b.p2) return false;
        }
        return true;
    }

    void build(int n, const std::vector<AxisNet> &nets) {
        n_ = n;
        nets_ = nets;
        ns_.reset(n + 1);
        ns_.reserve_arcs(4 * (int)nets.size() + 2 * n + 4 * n);
        int root = n;
        net_arc0_ = 0;
        for (const AxisNet &nt : nets_) {
            int b = nt.b >= 0 ? nt.b : root;
            for (int h = 0; h < 2; ++h) {
                ns_.add_arc(nt.a, b, nt.weight, 0);
                ns_.add_arc(b, nt.a, nt.weight, 0);
            }
        }
        bound_arc0_ = ns_.num_arcs();
        for (int i = 0; i < n; ++i) {
            ns_.add_arc(root, i, NetworkSimplex::INF, 0);
            ns_.add_arc(i, root, NetworkSimplex::INF, 0);
        }
        slots_.clear();
        free_.clear();
        slot_of_.assign((size_t)n * n, -1);
        last_pivots_ = ns_.pivots();
    }

    void set_edges(const std::vector<Edge> &edges) {
        for (Slot &s : slots_) s.alive = false;
        for (const Edge &e : edges) {
            int key = e.first * n_ + e.second;
            int sl = slot_of_[key];
            if (sl >= 0) { slots_[sl].alive = true; continue; }
            // retarget a dead slot whose arc carries no flow, else append
            int got = -1;
            while (!free_.empty() && got < 0) {
                int c = free_.back(); free_.pop_back();
                if (slots_[c].alive || slot_of_[slots_[c].i * n_ + slots_[c].j] != c) continue;
                if (ns_.is_basic(slots_[c].arc) || ns_.flow(slots_[c].arc) != 0) continue;
                slot_of_[slots_[c].i * n_ + slots_[c].j] = -1;
                ns_.set_endpoints(slots_[c].arc, e.first, e.second);
                got = c;
            }
            if (got < 0) {
                got = (int)slots_.size();
                slots_.push_back({0, 0, ns_.add_arc(e.first, e.second, NetworkSimplex::INF, 0), false});
            }
            slots_[got].i = e.first; slots_[got].j = e.second; slots_[got].alive = true;
            slot_of_[key] = got;
        }
        free_.clear();
        for (int s = 0; s < (int)slots_.size(); ++s)
            if (!slots_[s].alive) free_.push_back(s);
    }

    // ALAP crash basis: x_i = min(hi_i, min_{i->j} x_j - size_i); each node's
    // tree arc is its binding upper-bound arc i->r or edge arc i->j.
    void make_crash() {
        int n = n_;
        std::vector<int> start(n + 1, 0), adj, indeg(n, 0), order;
        for (const Slot &s : slots_) if (s.alive) { ++start[s.i + 1]; ++indeg[s.j]; }
        for (int i = 0; i < n; ++i) start[i + 1] += start[i];
        adj.resize(start[n]);
        {
            std::vector<int> fill(start.begin(), start.end() - 1);
            for (int k = 0; k < (int)slots_.size(); ++k)
                if (slots_[k].alive) adj[fill[slots_[k].i]++] = k;
        }
        order.reserve(n);
        for (int i = 0; i < n; ++i) if (indeg[i] == 0) order.push_back(i);
        for (size_t k = 0; k < order.size(); ++k)
            for (int a = start[order[k]]; a < start[order[k] + 1]; ++a)
                if (--indeg[slots_[adj[a]].j] == 0) order.push_back(slots_[adj[a]].j);
        if ((int)order.size() != n) { ns_.set_crash_tree({}); return; }
        std::vector<int64_t> x(n);
        std::vector<int> up(n + 1, -1);
        for (int k = n - 1; k >= 0; --k) {
            int i = order[k];
            x[i] = hi_[i];
            up[i] = bound_arc0_ + 2 * i + 1;          // i -> r
            for (int a = start[i]; a < start[i + 1]; ++a) {
                const Slot &s = slots_[adj[a]];
                if (x[s.j] - size_[i] < x[i]) { x[i] = x[s.j] - size_[i]; up[i] = s.arc; }
            }
        }
        ns_.set_crash_tree(up);                        // root r hangs from the artificial root
    }

    void update_costs() {
        int arc = net_arc0_;
        for (const AxisNet &nt : nets_) {
            int64_t delta = nt.b >= 0 ? size_[nt.a] - size_[nt.b] : size_[nt.a] - nt.p2;
            int64_t k1 = floor_div2(delta), k2 = delta - k1;
            ns_.set_cost(arc++, -k1); ns_.set_cost(arc++, k1);
            ns_.set_cost(arc++, -k2); ns_.set_cost(arc++, k2);
        }
        arc = bound_arc0_;
        int64_t mn = 0, mx = 0;
        for (int i = 0; i < n_; ++i) {
            ns_.set_cost(arc++, -lo_[i]);
            ns_.set_cost(arc++, hi_[i]);
            if (i == 0 || lo_[i] < mn) mn = lo_[i];
            if (i == 0 || hi_[i] > mx) mx = hi_[i];
        }
        // A dead edge arc i->j gets cost M: x_j - x_i >= -M always holds when
        // M >= max(hi) - min(lo), so it never constrains or carries flow.
        int64_t dead = (mx - mn) + 1;
        if (dead < 1) dead = 1;
        for (const Slot &s : slots_)
            ns_.set_cost(s.arc, s.alive ? -size_[s.i] : dead);
    }
};

// One-shot solve (cold).
inline AxisResult solve_axis(const AxisProblem &p) {
    AxisSolver s;
    s.set_problem(p);
    return s.solve(false);
}

// ------------------------------------------------------------ topology

struct ConstraintGraphs {
    std::vector<Edge> h, v;   // h: i left of j (x); v: i below j (y)
};

// For every pair: positive y-overlap -> H edge, positive x-overlap -> V edge,
// otherwise (diagonal; touching counts as disjoint) the larger-gap axis
// (H on ties), or diag_horiz(i, j) if given (called for i < j in index
// order). Edges follow the geometric order. Returns false on an overlap.
inline bool build_constraint_graphs(const std::vector<lab3::Rect> &r, ConstraintGraphs &g,
                                    const std::function<bool(int, int)> &diag_horiz = {}) {
    int n = (int)r.size();
    g.h.clear(); g.v.clear();
    g.h.reserve((size_t)n * (n - 1) / 2); g.v.reserve((size_t)n * (n - 1) / 2);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int64_t ox = std::min(r[i].x + r[i].w, r[j].x + r[j].w) - std::max(r[i].x, r[j].x);
            int64_t oy = std::min(r[i].y + r[i].h, r[j].y + r[j].h) - std::max(r[i].y, r[j].y);
            bool horiz;
            if (ox > 0 && oy > 0) return false;
            if (oy > 0) horiz = true;
            else if (ox > 0) horiz = false;
            else if (diag_horiz) horiz = diag_horiz(i, j);
            else horiz = -ox >= -oy;
            if (horiz) g.h.push_back(r[i].x + r[i].w <= r[j].x ? Edge(i, j) : Edge(j, i));
            else g.v.push_back(r[i].y + r[i].h <= r[j].y ? Edge(i, j) : Edge(j, i));
        }
    }
    return true;
}

// Transitive reduction of a DAG (bitset reachability). Duplicate edges are
// merged. Returns an empty list and sets *cyclic if the graph has a cycle.
inline std::vector<Edge> transitive_reduction(int n, const std::vector<Edge> &edges, bool *cyclic = nullptr) {
    // CSR adjacency + Kahn topological order
    std::vector<int> start(n + 1, 0), adj(edges.size()), indeg(n, 0);
    for (const Edge &e : edges) { ++start[e.first + 1]; ++indeg[e.second]; }
    for (int i = 0; i < n; ++i) start[i + 1] += start[i];
    {
        std::vector<int> fill(start.begin(), start.end() - 1);
        for (const Edge &e : edges) adj[fill[e.first]++] = e.second;
    }
    std::vector<int> topo, pos(n);
    topo.reserve(n);
    for (int i = 0; i < n; ++i) if (indeg[i] == 0) topo.push_back(i);
    for (size_t k = 0; k < topo.size(); ++k)
        for (int a = start[topo[k]]; a < start[topo[k] + 1]; ++a)
            if (--indeg[adj[a]] == 0) topo.push_back(adj[a]);
    if (cyclic) *cyclic = (int)topo.size() != n;
    if ((int)topo.size() != n) return {};
    for (int k = 0; k < n; ++k) pos[topo[k]] = k;
    // Bitsets indexed by topological position: succ[u] (direct successors),
    // reach[k] (nodes reachable from topo[k]). Successors are visited in
    // increasing position, so a successor reachable through an earlier one
    // is already in reach when it is examined.
    int words = (n + 63) / 64;
    std::vector<uint64_t> succ((size_t)n * words, 0), reach((size_t)n * words, 0);
    for (const Edge &e : edges) {
        int p = pos[e.second];
        succ[(size_t)pos[e.first] * words + (p >> 6)] |= (uint64_t)1 << (p & 63);
    }
    std::vector<Edge> res;
    for (int k = n - 1; k >= 0; --k) {
        uint64_t *ru = &reach[(size_t)k * words];
        const uint64_t *su = &succ[(size_t)k * words];
        for (int w = 0; w < words; ++w) {
            for (uint64_t bits = su[w]; bits; bits &= bits - 1) {
                int p = w * 64 + __builtin_ctzll(bits);
                if (ru[p >> 6] >> (p & 63) & 1) continue;
                res.emplace_back(topo[k], topo[p]);
                const uint64_t *rp = &reach[(size_t)p * words];
                for (int x = 0; x < words; ++x) ru[x] |= rp[x];
                ru[p >> 6] |= (uint64_t)1 << (p & 63);
            }
        }
    }
    return res;
}

// ------------------------------------------------------------ instance glue

// Nets of one axis (0 = x, 1 = y) from the instance connections.
inline std::vector<AxisNet> instance_nets(const lab3::Instance &inst, int axis) {
    std::vector<AxisNet> nets;
    nets.reserve(inst.conns.size());
    for (const lab3::Conn &c : inst.conns) {
        AxisNet nt;
        int a = c.a, b = c.b;
        if (a < 0) std::swap(a, b);
        nt.a = a; nt.b = b; nt.weight = c.weight;
        if (b < 0) {
            const lab3::Terminal &t = inst.terminals[-b - 1];
            nt.b = -1;
            nt.p2 = 2 * (axis == 0 ? t.x : t.y);
        }
        nets.push_back(nt);
    }
    return nets;
}

// Bounds of one axis for shapes `rects`: [0, L - size]; preplaced pinned to
// the input position; a boundary edge is pinned only if `rects` already
// touches it (so `rects` stays feasible and the boundary score is kept).
inline void instance_bounds(const lab3::Instance &inst, const std::vector<lab3::Rect> &rects, int axis,
                            std::vector<int64_t> &lo, std::vector<int64_t> &hi) {
    int n = (int)inst.blocks.size();
    int64_t L = axis == 0 ? inst.W : inst.H;
    lo.assign(n, 0); hi.assign(n, 0);
    for (int i = 0; i < n; ++i) {
        int64_t s = axis == 0 ? rects[i].w : rects[i].h;
        lo[i] = 0; hi[i] = L - s;
        const lab3::Block &b = inst.blocks[i];
        if (b.kind == lab3::Kind::Preplaced) lo[i] = hi[i] = axis == 0 ? b.x : b.y;
    }
    int lowE = axis == 0 ? lab3::EDGE_LEFT : lab3::EDGE_BOTTOM;
    int highE = axis == 0 ? lab3::EDGE_RIGHT : lab3::EDGE_TOP;
    for (const lab3::Boundary &bd : inst.boundaries) {
        int i = bd.block;
        int64_t p = axis == 0 ? rects[i].x : rects[i].y;
        int64_t s = axis == 0 ? rects[i].w : rects[i].h;
        if ((bd.edges & lowE) && p == 0) lo[i] = hi[i] = 0;
        if ((bd.edges & highE) && p + s == L) lo[i] = hi[i] = L - s;
    }
}

// Verifier gate items 3, 4 (shape/position), 7, 8 for `out` with the shapes
// of `in` (the soft-shape checks are unchanged because shapes are unchanged).
inline bool gate_ok(const lab3::Instance &inst, const std::vector<lab3::Rect> &in,
                    const std::vector<lab3::Rect> &out) {
    int n = (int)inst.blocks.size();
    if ((int)out.size() != n) return false;
    for (int i = 0; i < n; ++i) {
        const lab3::Rect &r = out[i];
        if (r.w != in[i].w || r.h != in[i].h || r.w <= 0 || r.h <= 0) return false;
        const lab3::Block &b = inst.blocks[i];
        if (b.kind == lab3::Kind::Preplaced && (r.x != b.x || r.y != b.y)) return false;
        if (r.x < 0 || r.y < 0 || r.x + r.w > inst.W || r.y + r.h > inst.H) return false;
    }
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            const lab3::Rect &a = out[i], &b = out[j];
            if (std::max(a.x, b.x) < std::min(a.x + a.w, b.x + b.w) &&
                std::max(a.y, b.y) < std::min(a.y + a.h, b.y + b.h)) return false;
        }
    return true;
}

struct PolishInfo {
    bool ok = false;        // true: result is the LP optimum and passed the gate check
    int64_t wl2_in = 0, wl2_out = 0;
    int edges_h = 0, edges_v = 0;   // constraint edges used (polish(): after reduction)
};

// Exact WL minimization on caller-supplied constraint graphs `g` (full or
// reduced; e.g. from a sequence pair). `rects` gives the shapes and should
// satisfy g; the result passes the gate, keeps every boundary edge `rects`
// touches and every MIB shape, and has WL <= WL(rects). If anything fails
// (illegal input, infeasible graphs, a check below), `rects` is returned
// unchanged with info->ok = false.
inline std::vector<lab3::Rect> polish_topology(const lab3::Instance &inst, const std::vector<lab3::Rect> &rects,
                                               const ConstraintGraphs &g, PolishInfo *info = nullptr) {
    PolishInfo loc;
    PolishInfo &pi = info ? *info : loc;
    pi = PolishInfo();
    pi.wl2_in = pi.wl2_out = lab3::wirelength2(inst, rects);
    int n = (int)inst.blocks.size();
    if ((int)rects.size() != n || !gate_ok(inst, rects, rects)) return rects;
    std::vector<lab3::Rect> out = rects;
    int64_t obj2 = 0;
    for (int axis = 0; axis < 2; ++axis) {
        AxisProblem p;
        p.size.resize(n);
        for (int i = 0; i < n; ++i) p.size[i] = axis == 0 ? rects[i].w : rects[i].h;
        instance_bounds(inst, rects, axis, p.lo, p.hi);
        p.edges = axis == 0 ? g.h : g.v;
        (axis == 0 ? pi.edges_h : pi.edges_v) = (int)p.edges.size();
        p.nets = instance_nets(inst, axis);
        AxisResult r = solve_axis(p);
        if (!r.ok) return rects;
        obj2 += r.obj2;
        for (int i = 0; i < n; ++i) (axis == 0 ? out[i].x : out[i].y) = r.x[i];
    }
    int64_t wl2 = lab3::wirelength2(inst, out);
    if (wl2 != obj2 || !gate_ok(inst, rects, out) || wl2 > pi.wl2_in) return rects;
    for (const lab3::Boundary &bd : inst.boundaries) {   // satisfied edges stay satisfied
        const lab3::Rect &a = rects[bd.block], &b = out[bd.block];
        if ((bd.edges & lab3::EDGE_LEFT) && a.x == 0 && b.x != 0) return rects;
        if ((bd.edges & lab3::EDGE_BOTTOM) && a.y == 0 && b.y != 0) return rects;
        if ((bd.edges & lab3::EDGE_RIGHT) && a.x + a.w == inst.W && b.x + b.w != inst.W) return rects;
        if ((bd.edges & lab3::EDGE_TOP) && a.y + a.h == inst.H && b.y + b.h != inst.H) return rects;
    }
    pi.ok = true;
    pi.wl2_out = wl2;
    return out;
}

// Exact WL minimization on the topology of the legal placement `rects`
// (H/V graphs by the geometry rule, larger-gap diagonals, transitively
// reduced), shapes kept. Same guarantees as polish_topology().
inline std::vector<lab3::Rect> polish(const lab3::Instance &inst, const std::vector<lab3::Rect> &rects,
                                      PolishInfo *info = nullptr) {
    ConstraintGraphs g;
    if ((int)rects.size() != (int)inst.blocks.size() || !build_constraint_graphs(rects, g)) {
        if (info) { *info = PolishInfo(); info->wl2_in = info->wl2_out = lab3::wirelength2(inst, rects); }
        return rects;
    }
    int n = (int)rects.size();
    g.h = transitive_reduction(n, g.h);
    g.v = transitive_reduction(n, g.v);
    return polish_topology(inst, rects, g, info);
}

}  // namespace mcf

#endif  // MCF_WL_LP_H
