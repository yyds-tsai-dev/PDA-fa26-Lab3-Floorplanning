// test_ns.cpp -- network simplex vs a successive-shortest-path oracle on
// random instances, plus an LP optimality certificate on every solve and
// warm-start re-solves vs cold solves.
#include "network_simplex.h"

#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

using mcf::NetworkSimplex;
using i64 = int64_t;

struct Arc { int u, v; i64 cap, cost; };
struct Inst { int n; std::vector<Arc> arcs; std::vector<i64> supply; };

static int g_fail = 0;
#define CHECK(c, ...) do { if (!(c)) { ++g_fail; std::printf("FAIL %s:%d: ", __FILE__, __LINE__); \
    std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

// Oracle: saturate negative arcs, then SSP with Bellman-Ford from a super
// source. INF caps are replaced by BIG (callers make sure no negative cycle
// consists of INF arcs only, so the optimum is unchanged). Returns false if
// infeasible.
static bool oracle(const Inst &in, i64 &best) {
    i64 big = 1;
    for (const Arc &a : in.arcs) if (a.cap < NetworkSimplex::INF) big += a.cap;
    for (i64 s : in.supply) big += s > 0 ? s : -s;
    int n = in.n, S = n, T = n + 1, N = n + 2;
    struct E { int to; i64 cap, cost; int rev; };
    std::vector<std::vector<E>> g(N);
    auto add = [&](int u, int v, i64 cap, i64 cost) {
        g[u].push_back({v, cap, cost, (int)g[v].size()});
        g[v].push_back({u, 0, -cost, (int)g[u].size() - 1});
    };
    std::vector<i64> ex(in.supply);
    i64 cost = 0;
    for (const Arc &a : in.arcs) {
        i64 cap = a.cap >= NetworkSimplex::INF ? big : a.cap;
        if (a.cost < 0) {  // saturate: reversed residual arc with positive cost
            cost += cap * a.cost;
            ex[a.u] -= cap; ex[a.v] += cap;
            add(a.v, a.u, cap, -a.cost);
        } else {
            add(a.u, a.v, cap, a.cost);
        }
    }
    i64 need = 0;
    for (int u = 0; u < n; ++u) {
        if (ex[u] > 0) { add(S, u, ex[u], 0); need += ex[u]; }
        else if (ex[u] < 0) add(u, T, -ex[u], 0);
    }
    const i64 BIGD = (i64)1 << 60;
    while (need > 0) {
        std::vector<i64> d(N, BIGD);
        std::vector<int> pv(N, -1), pe(N, -1);
        d[S] = 0;
        for (int it = 0; it < N; ++it) {
            bool ch = false;
            for (int u = 0; u < N; ++u) {
                if (d[u] == BIGD) continue;
                for (int k = 0; k < (int)g[u].size(); ++k) {
                    const E &e = g[u][k];
                    if (e.cap > 0 && d[u] + e.cost < d[e.to]) {
                        d[e.to] = d[u] + e.cost; pv[e.to] = u; pe[e.to] = k; ch = true;
                    }
                }
            }
            if (!ch) break;
        }
        if (d[T] == BIGD) return false;
        i64 f = need;
        for (int v = T; v != S; v = pv[v]) f = std::min(f, g[pv[v]][pe[v]].cap);
        for (int v = T; v != S; v = pv[v]) {
            E &e = g[pv[v]][pe[v]];
            e.cap -= f;
            g[v][e.rev].cap += f;
        }
        need -= f;
        cost += f * d[T];
    }
    best = cost;
    return true;
}

// Feasibility + complementary slackness => optimal (LP duality).
static bool certificate(const Inst &in, const NetworkSimplex &ns, i64 &cost) {
    std::vector<i64> bal(in.n, 0);
    cost = 0;
    for (int e = 0; e < (int)in.arcs.size(); ++e) {
        const Arc &a = in.arcs[e];
        i64 f = ns.flow(e);
        if (f < 0 || f > a.cap) return false;
        bal[a.u] += f; bal[a.v] -= f;
        cost += f * a.cost;
        i64 rc = a.cost + ns.potential(a.u) - ns.potential(a.v);
        if (f < a.cap && rc < 0) return false;
        if (f > 0 && rc > 0) return false;
    }
    for (int u = 0; u < in.n; ++u) if (bal[u] != in.supply[u]) return false;
    return cost == ns.total_cost();
}

static void load(NetworkSimplex &ns, const Inst &in) {
    ns.reset(in.n);
    for (const Arc &a : in.arcs) ns.add_arc(a.u, a.v, a.cap, a.cost);
    for (int u = 0; u < in.n; ++u) ns.set_supply(u, in.supply[u]);
}

// Random instance; INF arcs get non-negative cost so the problem is bounded.
static Inst random_inst(std::mt19937_64 &rng, bool circulation) {
    Inst in;
    in.n = 2 + (int)(rng() % 11);
    int m = (int)(rng() % (4 * in.n + 1));
    for (int k = 0; k < m; ++k) {
        Arc a;
        a.u = (int)(rng() % in.n);
        a.v = (int)(rng() % in.n);
        if (a.u == a.v) a.v = (a.u + 1) % in.n;
        bool inf = rng() % 5 == 0;
        a.cap = inf ? NetworkSimplex::INF : (i64)(rng() % 8);
        a.cost = (i64)(rng() % 41) - 20;
        if (inf && a.cost < 0) a.cost = -a.cost;
        in.arcs.push_back(a);
        if (rng() % 6 == 0) in.arcs.push_back(a);  // parallel arc
    }
    in.supply.assign(in.n, 0);
    if (!circulation) {
        for (int k = 0; k < in.n; ++k) {
            int u = (int)(rng() % in.n), v = (int)(rng() % in.n);
            i64 s = (i64)(rng() % 6);
            in.supply[u] += s; in.supply[v] -= s;
        }
    }
    return in;
}

int main() {
    std::mt19937_64 rng(12345);
    int solved = 0, infeasible = 0, warm_checks = 0;
    for (int it = 0; it < 20000; ++it) {
        Inst in = random_inst(rng, it % 3 == 0);
        NetworkSimplex ns;
        load(ns, in);
        NetworkSimplex::Status st = ns.solve();
        i64 best = 0;
        bool feas = oracle(in, best);
        CHECK(st != NetworkSimplex::UNBOUNDED, "iter %d: unexpected UNBOUNDED", it);
        if (!feas) {
            CHECK(st == NetworkSimplex::INFEASIBLE, "iter %d: oracle infeasible, ns status %d", it, (int)st);
            ++infeasible;
            continue;
        }
        CHECK(st == NetworkSimplex::OPTIMAL, "iter %d: oracle feasible, ns status %d", it, (int)st);
        if (st != NetworkSimplex::OPTIMAL) continue;
        i64 c;
        CHECK(certificate(in, ns, c), "iter %d: certificate failed", it);
        CHECK(c == best, "iter %d: cost %lld vs oracle %lld", it, (long long)c, (long long)best);
        ++solved;
        // Warm re-solve after cost changes, appended arcs and retargeted arcs.
        for (int rep = 0; rep < 3; ++rep) {
            for (int e = 0; e < (int)in.arcs.size(); ++e) {
                if (rng() % 3) continue;
                Arc &a = in.arcs[e];
                a.cost = (i64)(rng() % 41) - 20;
                if (a.cap >= NetworkSimplex::INF && a.cost < 0) a.cost = -a.cost;
                ns.set_cost(e, a.cost);
            }
            if (rng() % 2) {
                Arc a{(int)(rng() % in.n), (int)(rng() % in.n), (i64)(rng() % 8), (i64)(rng() % 41) - 20};
                if (a.u == a.v) a.v = (a.u + 1) % in.n;
                in.arcs.push_back(a);
                ns.add_arc(a.u, a.v, a.cap, a.cost);
            }
            if (!in.arcs.empty() && rng() % 2) {
                int e = (int)(rng() % in.arcs.size());
                int u = (int)(rng() % in.n), v = (int)(rng() % in.n);
                if (u != v && !ns.is_basic(e) && ns.flow(e) == 0) {
                    CHECK(ns.set_endpoints(e, u, v), "set_endpoints refused");
                    in.arcs[e].u = u; in.arcs[e].v = v;
                }
            }
            st = ns.solve();
            feas = oracle(in, best);
            CHECK(feas == (st == NetworkSimplex::OPTIMAL), "warm iter %d: status %d vs oracle %d", it, (int)st, (int)feas);
            if (feas && st == NetworkSimplex::OPTIMAL) {
                CHECK(certificate(in, ns, c), "warm iter %d: certificate failed", it);
                CHECK(c == best, "warm iter %d: cost %lld vs oracle %lld", it, (long long)c, (long long)best);
                ++warm_checks;
            }
            if (st != NetworkSimplex::OPTIMAL) break;
        }
    }
    // Unbounded: a negative-cost cycle of INF arcs.
    {
        NetworkSimplex ns(3);
        ns.add_arc(0, 1, NetworkSimplex::INF, -2);
        ns.add_arc(1, 2, NetworkSimplex::INF, 1);
        ns.add_arc(2, 0, NetworkSimplex::INF, 0);
        ns.add_arc(0, 2, 5, -7);
        CHECK(ns.solve() == NetworkSimplex::UNBOUNDED, "expected UNBOUNDED");
        ns.set_cost(0, 2);   // now bounded; cold restart after UNBOUNDED
        CHECK(ns.solve() == NetworkSimplex::OPTIMAL, "expected OPTIMAL after fix");
        CHECK(ns.total_cost() == 5 * (-7), "cost %lld", (long long)ns.total_cost());
    }
    std::printf("test_ns: %d optimal (oracle-equal), %d infeasible agreed, %d warm re-solves checked, %d failures\n",
                solved, infeasible, warm_checks, g_fail);
    return g_fail ? 1 : 0;
}
