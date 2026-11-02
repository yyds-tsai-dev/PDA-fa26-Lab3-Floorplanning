// test_wl.cpp -- wirelength LP tests.
//   1. random tiny axis problems vs exhaustive enumeration (cold, full vs
//      reduced edges, and an AxisSolver warm-start sequence);
//   2. every cases/dev construction: the 8 topology variants of
//      tools/lp_polish_proto.py (Python RNG reproduced) solved and compared
//      with the scipy LP optimum in results.json (exact, 2*WL integers);
//      polish(construction) == variant 0; outputs checked by the verifier;
//   3. polish() of every E3 packed/LP/MCF placement: ok, WL <= input,
//      verifier FEASIBLE with full boundary and MIB.
// Usage: test_wl <repo_root> <e3_dir|-> <out_dir>
#include "wl_lp.h"
#include "test_util.h"

#include <cmath>
#include <cstdio>
#include <dirent.h>
#include <random>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace mcf;
using i64 = int64_t;

static int g_fail = 0;
#define CHECK(c, ...) do { if (!(c)) { ++g_fail; std::printf("FAIL %s:%d: ", __FILE__, __LINE__); \
    std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

// ---------------------------------------------------------------- brute force

static bool brute(const AxisProblem &p, i64 &best) {
    int n = (int)p.size.size();
    std::vector<i64> x(n);
    for (int i = 0; i < n; ++i) { if (p.lo[i] > p.hi[i]) return false; x[i] = p.lo[i]; }
    bool found = false;
    for (;;) {
        bool ok = true;
        for (const Edge &e : p.edges) if (x[e.second] - x[e.first] < p.size[e.first]) { ok = false; break; }
        if (ok) {
            i64 s = 0;
            for (const AxisNet &nt : p.nets) {
                i64 a2 = 2 * x[nt.a] + p.size[nt.a];
                i64 b2 = nt.b >= 0 ? 2 * x[nt.b] + p.size[nt.b] : nt.p2;
                s += nt.weight * std::llabs(a2 - b2);
            }
            if (!found || s < best) best = s;
            found = true;
        }
        int k = 0;
        while (k < n && x[k] == p.hi[k]) { x[k] = p.lo[k]; ++k; }
        if (k == n) break;
        ++x[k];
    }
    return found;
}

static AxisProblem random_axis(std::mt19937_64 &rng) {
    AxisProblem p;
    int n = 1 + (int)(rng() % 4);
    i64 L = 4 + (i64)(rng() % 9);
    p.size.resize(n); p.lo.resize(n); p.hi.resize(n);
    for (int i = 0; i < n; ++i) {
        p.size[i] = 1 + (i64)(rng() % 4);
        p.lo[i] = 0; p.hi[i] = L - p.size[i];
        if (rng() % 6 == 0) p.lo[i] = p.hi[i] = (i64)(rng() % (L - p.size[i] + 1));
        else if (rng() % 6 == 0) p.lo[i] = p.hi[i];
    }
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && rng() % 4 == 0 && i < j) p.edges.emplace_back(rng() % 2 ? Edge(i, j) : Edge(j, i));
    int m = (int)(rng() % 6);
    for (int k = 0; k < m; ++k) {
        AxisNet nt;
        nt.a = (int)(rng() % n);
        nt.weight = 1 + (i64)(rng() % 10);
        if (n > 1 && rng() % 3) { nt.b = (int)(rng() % n); if (nt.b == nt.a) nt.b = (nt.a + 1) % n; }
        else { nt.b = -1; nt.p2 = 2 * (i64)(rng() % (L + 1)); }
        p.nets.push_back(nt);
    }
    return p;
}

static void test_brute() {
    std::mt19937_64 rng(777);
    int opt = 0, inf = 0, warm = 0;
    for (int it = 0; it < 20000; ++it) {
        AxisProblem p = random_axis(rng);
        i64 best = 0;
        bool feas = brute(p, best);
        AxisResult r = solve_axis(p);
        CHECK(r.ok == feas, "brute %d: ok %d vs brute %d", it, (int)r.ok, (int)feas);
        if (feas && r.ok) {
            CHECK(r.obj2 == best, "brute %d: obj2 %lld vs %lld", it, (long long)r.obj2, (long long)best);
            ++opt;
        } else if (!feas) ++inf;
        // transitive reduction gives the same optimum (when acyclic)
        bool cyc = false;
        AxisProblem q = p;
        q.edges = transitive_reduction((int)p.size.size(), p.edges, &cyc);
        if (!cyc) {
            AxisResult r2 = solve_axis(q);
            CHECK(r2.ok == r.ok && (!r.ok || r2.obj2 == r.obj2), "brute %d: TR mismatch", it);
        }
    }
    // Warm sequences: same nets, changing sizes / bounds / edges.
    for (int it = 0; it < 3000; ++it) {
        AxisProblem p = random_axis(rng);
        int n = (int)p.size.size();
        AxisSolver s;
        for (int step = 0; step < 6; ++step) {
            s.set_problem(p);
            AxisResult r = s.solve();
            i64 best = 0;
            bool feas = brute(p, best);
            CHECK(r.ok == feas && (!feas || r.obj2 == best), "warm %d/%d: ok %d/%d obj %lld vs %lld",
                  it, step, (int)r.ok, (int)feas, (long long)r.obj2, (long long)best);
            warm += r.warm;
            // mutate
            i64 L = p.hi[0] + p.size[0];
            if (p.lo[0] == p.hi[0]) L = 12;
            int i = (int)(rng() % n);
            if (rng() % 2) {
                p.size[i] = 1 + (i64)(rng() % 4);
                p.lo[i] = 0; p.hi[i] = std::max<i64>(0, L - p.size[i]);
            }
            if (rng() % 2 && n > 1) {
                int a = (int)(rng() % n), b = (int)(rng() % n);
                if (a != b) {
                    bool found = false;
                    for (size_t k = 0; k < p.edges.size(); ++k)
                        if ((p.edges[k].first == a && p.edges[k].second == b) ||
                            (p.edges[k].first == b && p.edges[k].second == a)) {
                            p.edges.erase(p.edges.begin() + k); found = true; break;
                        }
                    if (!found) p.edges.emplace_back(a, b);
                }
            }
        }
    }
    {   // contradictory bounds (lo > hi) and a self-edge are infeasible
        AxisProblem p;
        p.size = {2, 3}; p.lo = {5, 0}; p.hi = {4, 7};
        p.nets.push_back(AxisNet{0, 1, 3, 0});
        CHECK(!solve_axis(p).ok, "lo > hi accepted");
        p.lo = {0, 0}; p.edges = {Edge(1, 1)};
        CHECK(!solve_axis(p).ok, "self-edge accepted");
        p.edges = {Edge(0, 1), Edge(1, 0)};
        CHECK(!solve_axis(p).ok, "cycle accepted");
    }
    std::printf("brute force: %d optimal equal, %d infeasible agreed, %d warm solves\n", opt, inf, warm);
}

// ---------------------------------------------------------------- dev cases

static std::vector<std::string> list_dir(const std::string &d, const std::string &suffix) {
    std::vector<std::string> r;
    DIR *dp = opendir(d.c_str());
    if (!dp) return r;
    while (dirent *e = readdir(dp)) {
        std::string s = e->d_name;
        if (s.size() > suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0) r.push_back(s);
    }
    closedir(dp);
    std::sort(r.begin(), r.end());
    return r;
}

static bool file_exists(const std::string &p) { struct stat st; return stat(p.c_str(), &st) == 0; }

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

static i64 two(double wl) { return (i64)std::llround(2 * wl); }

int main(int argc, char **argv) {
    if (argc != 4) { std::fprintf(stderr, "usage: test_wl <repo_root> <e3_dir|-> <out_dir>\n"); return 2; }
    std::string root = argv[1], e3 = argv[2], outd = argv[3];
    std::string verifier = root + "/verifier/verifier", dev = root + "/cases/dev";
    mkdir(outd.c_str(), 0755);

    // Python RNG port sanity (values from CPython 3).
    {
        uint32_t c = tu::crc32("public_s1001/1");
        tu::PyRandom r(c);
        double a = r.random(), b = r.random();
        CHECK(c == 3816122789u, "crc32 %u", c);
        CHECK(a == 0.6245961245097115 && b == 0.5449374463212875, "PyRandom %.17g %.17g", a, b);
    }
    test_brute();

    tu::Json results;
    bool have_e3 = e3 != "-" && file_exists(e3 + "/results.json");
    if (have_e3) {
        std::ifstream f(e3 + "/results.json");
        std::stringstream ss; ss << f.rdbuf();
        std::string txt = ss.str();
        tu::JsonParser jp(txt);
        if (!jp.parse(results)) { CHECK(false, "cannot parse results.json"); have_e3 = false; }
    } else {
        std::printf("SKIP: E3 dir '%s' has no results.json; scipy comparisons and packing polish skipped\n", e3.c_str());
    }

    int exact = 0, compared = 0, ver_ok = 0, ver_runs = 0, cons_ok = 0, pol_ok = 0, pol_runs = 0, edge_cnt_ok = 0;
    int mcf_axis_ok = 0, full_edges_ok = 0;
    double gain_sum = 0;
    for (const std::string &fn : list_dir(dev, ".in")) {
        std::string base = fn.substr(0, fn.size() - 3), inp = dev + "/" + fn;
        lab3::Instance inst;
        std::string err;
        if (!lab3::parse_input_file(inp, inst, err)) { CHECK(false, "%s: %s", inp.c_str(), err.c_str()); continue; }
        std::vector<lab3::Rect> cons;
        if (!tu::read_output(dev + "/" + base + ".construction.out", inst, cons)) { CHECK(false, "%s: construction", base.c_str()); continue; }
        int n = (int)inst.blocks.size();
        i64 cons_wl2 = lab3::wirelength2(inst, cons);

        // polish(construction)
        PolishInfo info;
        std::vector<lab3::Rect> pc = polish(inst, cons, &info);
        CHECK(info.ok, "%s: polish(construction) not ok", base.c_str());
        CHECK(info.wl2_out <= cons_wl2, "%s: polish WL up", base.c_str());
        {
            std::string of = outd + "/" + base + ".polish.out";
            tu::write_output(of, inst, pc);
            tu::VerifierResult v = tu::run_verifier(verifier, inp, of);
            ++ver_runs;
            bool good = v.feasible && v.bnd_full && v.mib_full && v.computed_wl == lab3::format_half(info.wl2_out);
            CHECK(good, "%s: verifier on polish(construction) feas=%d bnd=%d mib=%d wl=%s", base.c_str(),
                  (int)v.feasible, (int)v.bnd_full, (int)v.mib_full, v.computed_wl.c_str());
            ver_ok += good;
        }
        gain_sum += (double)(cons_wl2 - info.wl2_out) / (double)cons_wl2;

        if (!have_e3) continue;
        const tu::Json *rec = nullptr;
        for (const tu::Json &r : results.arr) if (r["case"].str == base) rec = &r;
        if (!rec) { CHECK(false, "%s: not in results.json", base.c_str()); continue; }
        const std::vector<tu::Json> &vars = (*rec)["variants"].arr;

        // Sorted-name order, as in the prototype.
        std::vector<int> perm(n);
        for (int i = 0; i < n; ++i) perm[i] = i;
        std::sort(perm.begin(), perm.end(), [&](int a, int b) { return inst.blocks[a].name < inst.blocks[b].name; });
        std::vector<lab3::Rect> sr(n);
        for (int k = 0; k < n; ++k) sr[k] = cons[perm[k]];

        for (int v = 0; v < (int)vars.size(); ++v) {
            const tu::Json &vr = vars[v];
            ConstraintGraphs g;
            bool built;
            if (v == 0) built = build_constraint_graphs(sr, g);
            else {
                tu::PyRandom rng(tu::crc32(base + "/" + std::to_string(v)));
                built = build_constraint_graphs(sr, g, [&](int, int) { return rng.random() < 0.5; });
            }
            CHECK(built, "%s v%d: overlap", base.c_str(), v);
            for (Edge &e : g.h) e = Edge(perm[e.first], perm[e.second]);
            for (Edge &e : g.v) e = Edge(perm[e.first], perm[e.second]);
            std::vector<Edge> hr = transitive_reduction(n, g.h), vrd = transitive_reduction(n, g.v);
            bool cnt = (i64)g.h.size() == (i64)vr["edges_h"].num && (i64)g.v.size() == (i64)vr["edges_v"].num &&
                       (i64)hr.size() == (i64)vr["edges_h_tr"].num && (i64)vrd.size() == (i64)vr["edges_v_tr"].num;
            CHECK(cnt, "%s v%d: edge counts %zu/%zu/%zu/%zu vs %g/%g/%g/%g", base.c_str(), v, g.h.size(), g.v.size(),
                  hr.size(), vrd.size(), vr["edges_h"].num, vr["edges_v"].num, vr["edges_h_tr"].num, vr["edges_v_tr"].num);
            edge_cnt_ok += cnt;
            AxisResult rx = solve_axis(axis_problem(inst, cons, 0, hr));
            AxisResult ry = solve_axis(axis_problem(inst, cons, 1, vrd));
            CHECK(rx.ok && ry.ok, "%s v%d: LP infeasible", base.c_str(), v);
            if (!(rx.ok && ry.ok)) continue;
            ++compared;
            i64 want = two(vr["wl_lp"].num);
            bool eq = rx.obj2 + ry.obj2 == want;
            CHECK(eq, "%s v%d: 2WL %lld vs scipy %lld", base.c_str(), v, (long long)(rx.obj2 + ry.obj2), (long long)want);
            exact += eq;
            if (vr.has("mcf")) {
                bool ax = rx.obj2 == (i64)vr["mcf"]["opt_x"].num && ry.obj2 == (i64)vr["mcf"]["opt_y"].num;
                CHECK(ax, "%s v%d: per-axis %lld/%lld vs %g/%g", base.c_str(), v, (long long)rx.obj2,
                      (long long)ry.obj2, vr["mcf"]["opt_x"].num, vr["mcf"]["opt_y"].num);
                mcf_axis_ok += ax;
            }
            if (v == 0) {
                // full (unreduced) edge lists give the same optimum
                AxisResult fx = solve_axis(axis_problem(inst, cons, 0, g.h));
                AxisResult fy = solve_axis(axis_problem(inst, cons, 1, g.v));
                bool fe = fx.ok && fy.ok && fx.obj2 == rx.obj2 && fy.obj2 == ry.obj2;
                CHECK(fe, "%s: full-edge LP differs", base.c_str());
                full_edges_ok += fe;
                bool ce = info.wl2_out == want;
                CHECK(ce, "%s: polish(construction) %lld vs v0 %lld", base.c_str(), (long long)info.wl2_out, (long long)want);
                cons_ok += ce;
            }
            std::vector<lab3::Rect> out = cons;
            for (int i = 0; i < n; ++i) { out[i].x = rx.x[i]; out[i].y = ry.x[i]; }
            i64 wl2 = lab3::wirelength2(inst, out);
            CHECK(wl2 == want, "%s v%d: placement 2WL %lld vs %lld", base.c_str(), v, (long long)wl2, (long long)want);
            std::string of = outd + "/" + base + "_v" + std::to_string(v) + ".mcf.out";
            tu::write_output(of, inst, out);
            tu::VerifierResult vv = tu::run_verifier(verifier, inp, of);
            ++ver_runs;
            bool good = vv.feasible && vv.bnd_full && vv.mib_full && vv.computed_wl == lab3::format_half(want);
            CHECK(good, "%s v%d: verifier feas=%d bnd=%d mib=%d wl=%s", base.c_str(), v, (int)vv.feasible,
                  (int)vv.bnd_full, (int)vv.mib_full, vv.computed_wl.c_str());
            ver_ok += good;

            // polish() of the E3 placements of this variant
            for (const char *kind : {"packed", "lp", "mcf"}) {
                std::string pf = e3 + "/" + base + "_v" + std::to_string(v) + "." + kind + ".out";
                if (!file_exists(pf)) continue;
                std::vector<lab3::Rect> pr;
                if (!tu::read_output(pf, inst, pr)) { CHECK(false, "%s: read", pf.c_str()); continue; }
                PolishInfo pi;
                std::vector<lab3::Rect> po = polish(inst, pr, &pi);
                CHECK(pi.ok && pi.wl2_out <= pi.wl2_in, "%s: polish ok=%d %lld -> %lld", pf.c_str(), (int)pi.ok,
                      (long long)pi.wl2_in, (long long)pi.wl2_out);
                std::string of2 = outd + "/" + base + "_v" + std::to_string(v) + "." + kind + ".polish.out";
                tu::write_output(of2, inst, po);
                tu::VerifierResult pv = tu::run_verifier(verifier, inp, of2);
                ++pol_runs;
                bool pg = pi.ok && pv.feasible && pv.bnd_full && pv.mib_full &&
                          pv.computed_wl == lab3::format_half(pi.wl2_out) && pi.wl2_out <= pi.wl2_in;
                CHECK(pg, "%s: polished verifier feas=%d bnd=%d mib=%d", pf.c_str(), (int)pv.feasible,
                      (int)pv.bnd_full, (int)pv.mib_full);
                pol_ok += pg;
            }
        }
    }
    std::printf("construction polish: %d verifier-clean of %d runs incl. variants; mean gain %.2f%% over 20 cases\n",
                ver_ok, ver_runs, 100.0 * gain_sum / 20);
    if (have_e3) {
        std::printf("scipy exact match: %d / %d topology variants (2*WL integers); per-axis vs networkx MCF: %d; "
                    "edge counts equal: %d; full-edge == TR: %d; polish(construction) == v0 LP: %d\n",
                    exact, compared, mcf_axis_ok, edge_cnt_ok, full_edges_ok, cons_ok);
        std::printf("polish of E3 packed/lp/mcf placements: %d / %d ok + verifier FEASIBLE, full boundary/MIB, WL <= input\n",
                    pol_ok, pol_runs);
    }
    std::printf("test_wl: %d failures\n", g_fail);
    return g_fail ? 1 : 0;
}
