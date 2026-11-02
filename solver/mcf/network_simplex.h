// network_simplex.h -- integer primal network simplex for min-cost flow
// (header-only, C++17, no dependencies).
//
// Algorithm follows LEMON's NetworkSimplex (Kiraly & Kovacs): artificial root
// with one artificial arc per node, strongly feasible spanning tree
// (Cunningham's leaving-arc rule), block-search pivot rule with block size
// max(10, sqrt(m)). The tree is stored as parent pointers plus child lists;
// a pivot re-roots the leaving subtree at the entering arc and walks it once
// to update depth and potentials.
//
// Problem:   min  sum_e cost_e * f_e
//            s.t. 0 <= f_e <= cap_e                  (cap_e >= INF: no bound)
//                 sum_out f - sum_in f = supply_u    for every node u
// Supplies must sum to 0 (equality form). Potentials pi satisfy, at optimum,
// rc_e = cost_e + pi[src] - pi[dst] >= 0 if f_e < cap_e and <= 0 if f_e > 0.
//
// Warm start: costs may be changed with set_cost() and new arcs appended with
// add_arc() after a solve; the next solve() keeps the spanning tree and flows,
// recomputes potentials from the tree and resumes pivoting. set_endpoints()
// retargets a non-basic zero-flow arc in place. Any change that invalidates
// the tree (capacity change on a basic/upper arc, supply change, new nodes)
// falls back to a cold start automatically.
#ifndef MCF_NETWORK_SIMPLEX_H
#define MCF_NETWORK_SIMPLEX_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace mcf {

class NetworkSimplex {
public:
    using i64 = int64_t;
    static constexpr i64 INF = std::numeric_limits<i64>::max() / 4;
    enum Status { OPTIMAL, INFEASIBLE, UNBOUNDED };

    explicit NetworkSimplex(int nodes = 0) { reset(nodes); }

    // Removes all arcs and supplies; n user nodes 0..n-1.
    void reset(int nodes) {
        n_ = nodes;
        src_.assign(n_, 0); dst_.assign(n_, 0); cap_.assign(n_, INF);
        cost_.assign(n_, 0); flow_.assign(n_, 0); state_.assign(n_, 0);
        supply_.assign(n_ + 1, 0);
        pi_.assign(n_ + 1, 0);
        parent_.assign(n_ + 1, -1); pred_.assign(n_ + 1, -1);
        dir_.assign(n_ + 1, 0); depth_.assign(n_ + 1, 0);
        child_.assign(n_ + 1, -1); next_.assign(n_ + 1, -1); prev_.assign(n_ + 1, -1);
        stack_.clear(); path_.clear(); crash_.clear();
        has_tree_ = false;
        pivots_ = 0;
        warm_used_ = false;
    }

    int num_nodes() const { return n_; }
    int num_arcs() const { return (int)src_.size() - n_; }

    // Returns the arc id (0-based, dense). cap >= INF means uncapacitated.
    int add_arc(int u, int v, i64 cap, i64 cost) {
        src_.push_back(u); dst_.push_back(v);
        cap_.push_back(cap >= INF ? INF : cap);
        cost_.push_back(cost); flow_.push_back(0); state_.push_back(STATE_LOWER);
        return (int)src_.size() - 1 - n_;
    }
    void reserve_arcs(int m) {
        size_t t = (size_t)n_ + (size_t)m;
        src_.reserve(t); dst_.reserve(t); cap_.reserve(t);
        cost_.reserve(t); flow_.reserve(t); state_.reserve(t);
    }

    void set_cost(int e, i64 c) { cost_[n_ + e] = c; }
    i64 cost(int e) const { return cost_[n_ + e]; }

    void set_cap(int e, i64 c) {
        int a = n_ + e;
        if (c >= INF) c = INF;
        if (c == cap_[a]) return;
        // Harmless only for a non-basic arc with zero flow.
        if (state_[a] != STATE_TREE && flow_[a] == 0) state_[a] = STATE_LOWER;
        else has_tree_ = false;
        cap_[a] = c;
    }
    i64 cap(int e) const { return cap_[n_ + e]; }

    // Retargets arc e. Keeps the warm start iff e is non-basic with zero flow;
    // returns whether it did.
    bool set_endpoints(int e, int u, int v) {
        int a = n_ + e;
        bool keep = state_[a] != STATE_TREE && flow_[a] == 0;
        if (keep) state_[a] = STATE_LOWER;
        else has_tree_ = false;
        src_[a] = u; dst_[a] = v;
        return keep;
    }
    bool is_basic(int e) const { return state_[n_ + e] == STATE_TREE; }

    // Optional initial basis for cold starts (all supplies zero only):
    // up_arc[u] is a user arc with source u whose target becomes u's tree
    // parent, or -1 to hang u from the artificial root. All these arcs start
    // with zero flow and point toward the root, so the tree is strongly
    // feasible. Ignored if it is not a forest. Example: an ALAP longest-path
    // tree of difference constraints gives a near-optimal dual start.
    void set_crash_tree(const std::vector<int> &up_arc) { crash_ = up_arc; }

    void set_supply(int u, i64 s) {
        if (supply_[u] != s) has_tree_ = false;
        supply_[u] = s;
    }

    // Solves; warm-starts from the previous tree when allowed and valid.
    Status solve(bool warm = true) {
        int all = (int)src_.size();
        i64 sum = 0;
        for (int u = 0; u < n_; ++u) sum += supply_[u];
        if (sum != 0) return INFEASIBLE;
        bool neg = false;
        for (int u = 0; u < n_; ++u) neg = neg || supply_[u] < 0;
        i64 art = neg ? art_cost() : 0;
        warm_used_ = warm && has_tree_ && art <= art_cost_ && warm_valid(all);
        if (!warm_used_) {
            art_cost_ = art;
            cold_init();
        } else {
            for (int u = 0; u < n_; ++u)
                if (supply_[u] < 0) cost_[u] = art_cost_;
            recompute_potentials();
        }
        has_tree_ = true;
        block_ = std::max(10, (int)(block_factor_ * std::sqrt((double)std::max(1, all - n_))));
        if (next_arc_ < n_ || next_arc_ >= all) next_arc_ = n_;
        int in_arc;
        while (find_entering(all, in_arc)) {
            ++pivots_;
            if (!pivot(in_arc)) { has_tree_ = false; return UNBOUNDED; }
        }
        for (int u = 0; u < n_; ++u)
            if (flow_[u] != 0) return INFEASIBLE;
        return OPTIMAL;
    }

    i64 flow(int e) const { return flow_[n_ + e]; }
    // Potential of user node u (pi of the artificial root is 0).
    i64 potential(int u) const { return pi_[u]; }
    i64 total_cost() const {
        i64 c = 0;
        for (size_t a = n_; a < src_.size(); ++a) c += flow_[a] * cost_[a];
        return c;
    }
    long long pivots() const { return pivots_; }
    // Block-search size = factor * sqrt(#arcs) (LEMON default 1.0).
    void set_block_factor(double f) { block_factor_ = f; }
    bool last_solve_was_warm() const { return warm_used_; }
    bool has_tree() const { return has_tree_; }

private:
    enum { STATE_UPPER = -1, STATE_TREE = 0, STATE_LOWER = 1 };
    enum { DIR_DOWN = -1, DIR_UP = 1 };

    // Arcs [0, n_) are artificial (arc u joins node u and the root n_);
    // user arc e is stored at n_ + e.
    int n_ = 0;
    std::vector<int> src_, dst_;
    std::vector<i64> cap_, cost_, flow_;
    std::vector<signed char> state_;
    std::vector<i64> supply_, pi_;
    std::vector<int> parent_, pred_, dir_, depth_;
    std::vector<int> child_, next_, prev_;     // child lists
    std::vector<int> stack_, path_, crash_;
    i64 art_cost_ = 0;
    int block_ = 10, next_arc_ = 0;
    double block_factor_ = 1.0;
    bool has_tree_ = false, warm_used_ = false;
    long long pivots_ = 0;

    i64 art_cost() const {
        i64 mx = 0;
        for (size_t a = n_; a < src_.size(); ++a) mx = std::max(mx, cost_[a] < 0 ? -cost_[a] : cost_[a]);
        return (mx + 1) * (i64)(n_ + 1);
    }

    bool warm_valid(int all) const {
        if ((int)parent_.size() != n_ + 1) return false;
        for (int a = n_; a < all; ++a) {
            if (state_[a] == STATE_LOWER && flow_[a] != 0) return false;
            if (state_[a] == STATE_UPPER && flow_[a] != cap_[a]) return false;
            if (flow_[a] < 0 || flow_[a] > cap_[a]) return false;
        }
        return true;
    }

    void detach(int u) {
        int p = parent_[u];
        if (p < 0) return;
        if (prev_[u] >= 0) next_[prev_[u]] = next_[u]; else child_[p] = next_[u];
        if (next_[u] >= 0) prev_[next_[u]] = prev_[u];
        next_[u] = prev_[u] = -1;
        parent_[u] = -1;
    }
    void attach(int u, int p) {
        parent_[u] = p;
        prev_[u] = -1;
        next_[u] = child_[p];
        if (child_[p] >= 0) prev_[child_[p]] = u;
        child_[p] = u;
    }

    void cold_init() {
        int root = n_;
        std::fill(child_.begin(), child_.end(), -1);
        for (size_t a = n_; a < src_.size(); ++a) { flow_[a] = 0; state_[a] = STATE_LOWER; }
        parent_[root] = -1; pred_[root] = -1; depth_[root] = 0; pi_[root] = 0;
        next_[root] = prev_[root] = -1;
        for (int u = 0; u < n_; ++u) {
            parent_[u] = -1;
            attach(u, root);
            pred_[u] = u;
            depth_[u] = 1;
            cap_[u] = INF;
            state_[u] = STATE_TREE;
            if (supply_[u] >= 0) {
                dir_[u] = DIR_UP; src_[u] = u; dst_[u] = root;
                flow_[u] = supply_[u]; cost_[u] = 0; pi_[u] = 0;
            } else {
                dir_[u] = DIR_DOWN; src_[u] = root; dst_[u] = u;
                flow_[u] = -supply_[u]; cost_[u] = art_cost_; pi_[u] = art_cost_;
            }
        }
        next_arc_ = n_;
        if (crash_ok()) {
            for (int u = 0; u < n_; ++u) {
                if (crash_[u] < 0) continue;
                int a = n_ + crash_[u];
                state_[u] = STATE_LOWER;            // artificial arc leaves the basis
                detach(u);
                attach(u, dst_[a]);
                pred_[u] = a; dir_[u] = DIR_UP; state_[a] = STATE_TREE;
            }
            recompute_potentials();
        }
    }

    bool crash_ok() {
        if ((int)crash_.size() != n_) return false;
        int m = (int)src_.size() - n_;
        for (int u = 0; u < n_; ++u) {
            if (supply_[u] != 0) return false;
            int e = crash_[u];
            if (e < 0) continue;
            if (e >= m || src_[n_ + e] != u || dst_[n_ + e] == u) return false;
        }
        // forest check: follow parent pointers with colours (0 new, 1 on path, 2 done)
        std::vector<signed char> col(n_, 0);
        for (int s = 0; s < n_; ++s) {
            int u = s;
            while (u >= 0 && u < n_ && col[u] == 0) {
                col[u] = 1;
                u = crash_[u] < 0 ? -1 : dst_[n_ + crash_[u]];
            }
            if (u >= 0 && u < n_ && col[u] == 1) return false;   // cycle
            for (u = s; u >= 0 && u < n_ && col[u] == 1; u = crash_[u] < 0 ? -1 : dst_[n_ + crash_[u]]) col[u] = 2;
        }
        return true;
    }

    // pi from the tree (rc = 0 on tree arcs), depth refreshed too.
    void recompute_potentials() {
        int root = n_;
        pi_[root] = 0; depth_[root] = 0;
        stack_.clear();
        for (int c = child_[root]; c >= 0; c = next_[c]) stack_.push_back(c);
        while (!stack_.empty()) {
            int u = stack_.back(); stack_.pop_back();
            int p = parent_[u], e = pred_[u];
            pi_[u] = dir_[u] == DIR_UP ? pi_[p] - cost_[e] : pi_[p] + cost_[e];
            depth_[u] = depth_[p] + 1;
            for (int c = child_[u]; c >= 0; c = next_[c]) stack_.push_back(c);
        }
    }

    bool find_entering(int all, int &in_arc) {
        i64 mn = 0;
        int cnt = block_, e;
        for (e = next_arc_; e < all; ++e) {
            i64 c = state_[e] * (cost_[e] + pi_[src_[e]] - pi_[dst_[e]]);
            if (c < mn) { mn = c; in_arc = e; }
            if (--cnt == 0) {
                if (mn < 0) goto found;
                cnt = block_;
            }
        }
        for (e = n_; e < next_arc_; ++e) {
            i64 c = state_[e] * (cost_[e] + pi_[src_[e]] - pi_[dst_[e]]);
            if (c < mn) { mn = c; in_arc = e; }
            if (--cnt == 0) {
                if (mn < 0) goto found;
                cnt = block_;
            }
        }
        if (mn >= 0) return false;
    found:
        next_arc_ = e;
        return true;
    }

    // One pivot. Returns false if the cycle has unbounded residual capacity.
    bool pivot(int in_arc) {
        int first, second;
        if (state_[in_arc] == STATE_LOWER) { first = src_[in_arc]; second = dst_[in_arc]; }
        else { first = dst_[in_arc]; second = src_[in_arc]; }
        // join node
        int a = first, b = second;
        while (a != b) {
            if (depth_[a] >= depth_[b]) a = parent_[a]; else b = parent_[b];
        }
        int join = a;
        // leaving arc (Cunningham: strict on first side, non-strict on second)
        i64 delta = cap_[in_arc];
        int result = 0, u_out = -1;
        for (int u = first; u != join; u = parent_[u]) {
            int e = pred_[u];
            i64 d = flow_[e];
            if (dir_[u] == DIR_DOWN) d = cap_[e] >= INF ? INF : cap_[e] - d;
            if (d < delta) { delta = d; u_out = u; result = 1; }
        }
        for (int u = second; u != join; u = parent_[u]) {
            int e = pred_[u];
            i64 d = flow_[e];
            if (dir_[u] == DIR_UP) d = cap_[e] >= INF ? INF : cap_[e] - d;
            if (d <= delta) { delta = d; u_out = u; result = 2; }
        }
        if (delta >= INF) return false;
        // flow change
        if (delta > 0) {
            i64 val = state_[in_arc] * delta;
            flow_[in_arc] += val;
            for (int u = src_[in_arc]; u != join; u = parent_[u]) flow_[pred_[u]] -= dir_[u] * val;
            for (int u = dst_[in_arc]; u != join; u = parent_[u]) flow_[pred_[u]] += dir_[u] * val;
        }
        if (result == 0) {
            state_[in_arc] = (signed char)-state_[in_arc];
            return true;
        }
        // tree update: subtree under u_out is re-rooted at u_in, hung from v_in
        int u_in = result == 1 ? first : second;
        int v_in = result == 1 ? second : first;
        int out_arc = pred_[u_out];
        state_[out_arc] = flow_[out_arc] == 0 ? STATE_LOWER : STATE_UPPER;
        state_[in_arc] = STATE_TREE;
        path_.clear();
        for (int u = u_in; ; u = parent_[u]) { path_.push_back(u); if (u == u_out) break; }
        int k = (int)path_.size();
        std::vector<int> &pp = path_;
        // save old pred/dir along the path before overwriting
        int prev_pred = pred_[pp[0]], prev_dir = dir_[pp[0]];
        for (int i = 0; i < k; ++i) detach(pp[i]);
        attach(pp[0], v_in);
        pred_[pp[0]] = in_arc;
        dir_[pp[0]] = src_[in_arc] == pp[0] ? DIR_UP : DIR_DOWN;
        for (int i = 1; i < k; ++i) {
            int u = pp[i];
            int sp = pred_[u], sd = dir_[u];
            attach(u, pp[i - 1]);
            pred_[u] = prev_pred;
            dir_[u] = -prev_dir;
            prev_pred = sp; prev_dir = sd;
        }
        // potentials and depth of the moved subtree
        int r = pp[0];
        i64 newpi = dir_[r] == DIR_UP ? pi_[v_in] - cost_[in_arc] : pi_[v_in] + cost_[in_arc];
        i64 sigma = newpi - pi_[r];
        stack_.clear();
        stack_.push_back(r);
        while (!stack_.empty()) {
            int u = stack_.back(); stack_.pop_back();
            pi_[u] += sigma;
            depth_[u] = depth_[parent_[u]] + 1;
            for (int c = child_[u]; c >= 0; c = next_[c]) stack_.push_back(c);
        }
        return true;
    }
};

}  // namespace mcf

#endif  // MCF_NETWORK_SIMPLEX_H
