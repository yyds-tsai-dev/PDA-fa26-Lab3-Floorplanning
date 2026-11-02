// test_timing.cpp -- per-axis solve times at n = 200 (cases/dev/hidden4_*),
// cold vs warm re-solve after a small change, and polish() end to end.
// Warm results are checked against cold solves (exact objective).
// Usage: test_timing <repo_root>
#include "wl_lp.h"
#include "test_util.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace mcf;
using i64 = int64_t;
using Clock = std::chrono::steady_clock;

static int g_fail = 0;
#define CHECK(c, ...) do { if (!(c)) { ++g_fail; std::printf("FAIL %s:%d: ", __FILE__, __LINE__); \
    std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

static double us_since(Clock::time_point t) {
    return std::chrono::duration<double, std::micro>(Clock::now() - t).count();
}

static void report(const char *what, std::vector<double> v) {
    std::sort(v.begin(), v.end());
    auto q = [&](double p) { return v[std::min(v.size() - 1, (size_t)(p * (double)v.size()))]; };
    std::printf("  %-44s n=%5zu  median %8.1f us  p90 %8.1f us  max %8.1f us\n", what, v.size(), q(0.5), q(0.9), v.back());
}

static AxisProblem axis_problem(const lab3::Instance &inst, const std::vector<lab3::Rect> &r, int axis,
                                const std::vector<Edge> &edges) {
    AxisProblem p;
    int n = (int)r.size();
    p.size.resize(n);
    for (int i = 0; i < n; ++i) p.size[i] = axis == 0 ? r[i].w : r[i].h;
    instance_bounds(inst, r, axis, p.lo, p.hi);
    p.edges = edges;
    p.nets = instance_nets(inst, axis);
    return p;
}

int main(int argc, char **argv) {
    if (argc != 2) { std::fprintf(stderr, "usage: test_timing <repo_root>\n"); return 2; }
    std::string dev = std::string(argv[1]) + "/cases/dev";
    std::vector<double> t_cold, t_cold_solve, t_build, t_polish, t_shape_cold, t_shape_warm, t_topo_cold, t_topo_warm;
    std::vector<double> piv_cold, piv_shape, piv_topo;
    int arcs_max = 0;
    std::mt19937_64 rng(2026);
    for (int s = 1001; s <= 1004; ++s) {
        std::string base = dev + "/hidden4_s" + std::to_string(s);
        lab3::Instance inst;
        std::string err;
        std::vector<lab3::Rect> cons;
        if (!lab3::parse_input_file(base + ".in", inst, err) || !tu::read_output(base + ".construction.out", inst, cons)) {
            CHECK(false, "cannot read %s", base.c_str());
            continue;
        }
        int n = (int)inst.blocks.size();
        for (int rep = 0; rep < 10; ++rep) {
            auto t0 = Clock::now();
            PolishInfo info;
            polish(inst, cons, &info);
            t_polish.push_back(us_since(t0) / 1000.0);
            CHECK(info.ok, "polish failed");
        }
        for (int v = 0; v < 8; ++v) {
            ConstraintGraphs g;
            auto t0 = Clock::now();
            if (v == 0) build_constraint_graphs(cons, g);
            else build_constraint_graphs(cons, g, [&](int, int) { return rng() & 1; });
            std::vector<Edge> E[2] = {transitive_reduction(n, g.h), transitive_reduction(n, g.v)};
            t_build.push_back(us_since(t0));
            for (int axis = 0; axis < 2; ++axis) {
                AxisProblem p = axis_problem(inst, cons, axis, E[axis]);
                AxisResult ref;
                for (int rep = 0; rep < 25; ++rep) {
                    auto t1 = Clock::now();
                    ref = solve_axis(p);
                    t_cold.push_back(us_since(t1));
                    AxisSolver as;
                    as.set_problem(p);
                    auto t2 = Clock::now();
                    AxisResult r2 = as.solve(false);
                    t_cold_solve.push_back(us_since(t2));
                    arcs_max = std::max(arcs_max, as.simplex().num_arcs());
                    if (rep == 0) piv_cold.push_back((double)r2.pivots);
                }
                CHECK(ref.ok, "cold infeasible");
                // Warm: small shape change (3 non-preplaced blocks shrunk by up to 10%).
                AxisSolver ws;
                ws.set_problem(p);
                ws.solve(false);
                for (int trial = 0; trial < 40; ++trial) {
                    AxisProblem q = p;
                    for (int k = 0; k < 3; ++k) {
                        int i = (int)(rng() % n);
                        if (inst.blocks[i].kind == lab3::Kind::Preplaced) continue;
                        i64 d = (i64)(rng() % (uint64_t)std::max<i64>(1, q.size[i] / 10)) + 1;   // up to 10%
                        q.size[i] = std::max<i64>(1, q.size[i] - d);   // shrink keeps feasibility
                    }
                    // recompute bounds for the new sizes, keep pins of boundary blocks
                    for (int i = 0; i < n; ++i) {
                        if (q.lo[i] == q.hi[i] && inst.blocks[i].kind != lab3::Kind::Preplaced) {
                            if (p.lo[i] != 0) q.lo[i] = q.hi[i] = p.hi[i] + p.size[i] - q.size[i];
                        } else if (inst.blocks[i].kind != lab3::Kind::Preplaced) {
                            q.hi[i] = p.hi[i] + p.size[i] - q.size[i];
                        }
                    }
                    auto t3 = Clock::now();
                    AxisResult c = solve_axis(q);
                    t_shape_cold.push_back(us_since(t3));
                    ws.set_problem(q);
                    auto t4 = Clock::now();
                    AxisResult w = ws.solve(true);
                    t_shape_warm.push_back(us_since(t4));
                    piv_shape.push_back((double)w.pivots);
                    CHECK(c.ok == w.ok && (!c.ok || c.obj2 == w.obj2), "shape warm mismatch");
                    CHECK(w.warm, "shape re-solve was not warm");
                    // restore p for the next trial (also a warm solve)
                    ws.set_problem(p);
                    ws.solve(true);
                }
            }
            // Warm: topology change -- flip one diagonal pair whose edge is in the
            // transitive reduction between H and V (infeasible flips are compared too).
            AxisSolver ws[2];
            for (int axis = 0; axis < 2; ++axis) {
                ws[axis].set_problem(axis_problem(inst, cons, axis, E[axis]));
                ws[axis].solve(false);
            }
            for (int trial = 0; trial < 40; ++trial) {
                int i, j;
                const lab3::Rect *a, *b;
                for (;;) {
                    i = (int)(rng() % n); j = (int)(rng() % n);
                    if (i == j) continue;
                    a = &cons[i]; b = &cons[j];
                    i64 ox = std::min(a->x + a->w, b->x + b->w) - std::max(a->x, b->x);
                    i64 oy = std::min(a->y + a->h, b->y + b->h) - std::max(a->y, b->y);
                    if (ox > 0 || oy > 0) continue;
                    // only pairs whose relation is a reduced (direct) edge, so the flip matters
                    Edge eh0 = a->x + a->w <= b->x ? Edge(i, j) : Edge(j, i);
                    Edge ev0 = a->y + a->h <= b->y ? Edge(i, j) : Edge(j, i);
                    if (std::find(E[0].begin(), E[0].end(), eh0) != E[0].end() ||
                        std::find(E[1].begin(), E[1].end(), ev0) != E[1].end()) break;
                }
                ConstraintGraphs g2 = g;
                Edge eh = a->x + a->w <= b->x ? Edge(i, j) : Edge(j, i);
                Edge ev = a->y + a->h <= b->y ? Edge(i, j) : Edge(j, i);
                auto ith = std::find(g2.h.begin(), g2.h.end(), eh);
                if (ith != g2.h.end()) { g2.h.erase(ith); g2.v.push_back(ev); }
                else { auto itv = std::find(g2.v.begin(), g2.v.end(), ev); if (itv != g2.v.end()) g2.v.erase(itv); g2.h.push_back(eh); }
                std::vector<Edge> E2[2] = {transitive_reduction(n, g2.h), transitive_reduction(n, g2.v)};
                for (int axis = 0; axis < 2; ++axis) {
                    AxisProblem q = axis_problem(inst, cons, axis, E2[axis]);
                    auto t3 = Clock::now();
                    AxisResult c = solve_axis(q);
                    t_topo_cold.push_back(us_since(t3));
                    ws[axis].set_problem(q);
                    auto t4 = Clock::now();
                    AxisResult w = ws[axis].solve(true);
                    t_topo_warm.push_back(us_since(t4));
                    piv_topo.push_back((double)w.pivots);
                    CHECK(c.ok == w.ok && (!c.ok || c.obj2 == w.obj2), "topology warm mismatch");
                    CHECK(w.warm, "topology re-solve was not warm");
                    ws[axis].set_problem(axis_problem(inst, cons, axis, E[axis]));
                    ws[axis].solve(true);
                }
            }
        }
    }
    auto med = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return v.empty() ? 0.0 : v[v.size() / 2]; };
    std::printf("timing at n=200 (hidden4 x 4 cases x 8 topologies x 2 axes), up to %d arcs per axis:\n", arcs_max);
    report("cold solve_axis (build + solve)", t_cold);
    report("cold solve only (AxisSolver::solve(false))", t_cold_solve);
    report("H/V graphs + 2x transitive reduction", t_build);
    report("warm re-solve, 3 sizes shrunk by <=10%", t_shape_warm);
    report("  cold solve_axis of the same problems", t_shape_cold);
    report("warm re-solve, 1 direct diagonal pair flipped", t_topo_warm);
    report("  cold solve_axis of the same problems", t_topo_cold);
    auto mean = [](const std::vector<double> &v) { double t = 0; for (double x : v) t += x; return v.empty() ? 0.0 : t / (double)v.size(); };
    auto nz = [](const std::vector<double> &v) { int c = 0; for (double x : v) c += x > 0; return c; };
    std::printf("  pivots median/mean: cold %.0f/%.0f, warm(shape) %.0f/%.1f (%d of %zu re-solves pivot), "
                "warm(flip) %.0f/%.1f (%d of %zu pivot)\n", med(piv_cold), mean(piv_cold), med(piv_shape),
                mean(piv_shape), nz(piv_shape), piv_shape.size(), med(piv_topo), mean(piv_topo), nz(piv_topo), piv_topo.size());
    std::vector<double> pms = t_polish;
    std::sort(pms.begin(), pms.end());
    std::printf("  polish() end to end (n=200, both axes, incl. gate check): median %.3f ms\n", pms[pms.size() / 2]);
    std::printf("test_timing: %d failures\n", g_fail);
    return g_fail ? 1 : 0;
}
