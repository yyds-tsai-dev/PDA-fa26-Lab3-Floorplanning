// test_construct.cpp -- harness for construct.h (legal-first constructive greedy).
//
//   test_construct <case.in> [--seed N=1] [--blind] [--score 0|1] [--repair] [--repair-ms MS=1500]
//                  [--out DIR] [--construction FILE] [--gp-iters N=25] [--gp-anchor X=0.05]
//
// Without --blind, targets are the GlobalPlacer leaf centres (G.tx, G.ty) of
// one GP run with Rng(seed) (as the first seed of fp's make_seeds). Prints one
// line per case:
//   CONS <tag> placed=<k>/<n> unpl=<u> qarea%=<..> maxmer=<w>x<h>(<area%>) nmer=<..>/<max>
//        gate=<0|1> wl=<WL/cons> bsat= msat= msgp= msgreedy= msrep= rep=<acc>/<tries>
// and, with --out, writes <DIR>/<tag>.out when every block is placed.
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "construct.h"

using namespace fp;

int main(int argc, char **argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s case.in [options]\n", argv[0]); return 2; }
    std::string in = argv[1], outdir, cons;
    cx::Opts o;
    int gp_iters = 25;
    double gp_anchor = 0.05;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto val = [&]() -> const char * { return i + 1 < argc ? argv[++i] : ""; };
        if (a == "--seed") o.seed = std::strtoull(val(), nullptr, 10);
        else if (a == "--blind") o.target = false;
        else if (a == "--score") o.score = std::atoi(val());
        else if (a == "--repair") o.repair = true;
        else if (a == "--order") o.order = std::atoi(val());
        else if (a == "--bestfit") o.bestfit = std::atoi(val());
        else if (a == "--evict-max") o.evict_max = std::atoi(val());
        else if (a == "--evict-fixed") o.evict_fixed = true;
        else if (a == "--plateau") o.plateau = true;
        else if (a == "--max-fails") o.max_fails = std::atoi(val());
        else if (a == "--repair-ms") o.repair_ms = std::atof(val());
        else if (a == "--out") outdir = val();
        else if (a == "--construction") cons = val();
        else if (a == "--gp-iters") gp_iters = std::atoi(val());
        else if (a == "--gp-anchor") gp_anchor = std::atof(val());
        else { std::fprintf(stderr, "bad option %s\n", a.c_str()); return 2; }
    }
    if (cons.empty()) {
        cons = in;
        if (cons.size() > 3 && cons.compare(cons.size() - 3, 3, ".in") == 0) cons.resize(cons.size() - 3);
        cons += ".construction.out";
    }
    std::string tag = in.substr(in.find_last_of('/') + 1);
    if (tag.size() > 3) tag.resize(tag.size() - 3);

    Problem P;
    std::string err;
    if (!build_problem(in, P, err)) { std::fprintf(stderr, "input: %s\n", err.c_str()); return 2; }
    std::vector<lab3::Rect> cr;
    double wcons = 0;
    if (read_placement(P, cons, cr, err)) wcons = (double)lab3::wirelength2(P.inst, cr);
    else std::fprintf(stderr, "construction: %s\n", err.c_str());

    Rng rng(o.seed);
    double ms_gp = 0;
    std::vector<double> tx, ty;
    if (o.target) {
        auto t0 = std::chrono::steady_clock::now();
        GlobalPlacer G(P, rng);
        G.run(gp_iters, gp_anchor);
        ms_gp = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        tx = G.tx;
        ty = G.ty;
    }
    cx::Builder B(P, o);
    if (o.target) { B.tx = &tx; B.ty = &ty; }
    cx::Stats st = B.run(rng);

    if (std::getenv("CX_CHECK")) {   // MER invariants + partial-placement legality
        int bad = 0;
        std::vector<cx::Rc> pl;
        for (int i = 0; i < P.n; ++i) if (B.placed[i]) pl.push_back(B.pos[i]);
        for (size_t i = 0; i < pl.size(); ++i)
            for (size_t j = i + 1; j < pl.size(); ++j) if (cx::ovl(pl[i], pl[j])) ++bad;
        auto freeR = [&](const cx::Rc &r) {
            if (r.x < 0 || r.y < 0 || r.x1() > P.W || r.y1() > P.H) return false;
            for (auto &q : pl) if (cx::ovl(q, r)) return false;
            return true;
        };
        int nonmax = 0, notfree = 0;
        for (const cx::Rc &F : B.M.f) {
            if (!freeR(F)) ++notfree;
            cx::Rc a = F; a.x -= 1; a.w += 1;
            cx::Rc b = F; b.w += 1;
            cx::Rc c = F; c.y -= 1; c.h += 1;
            cx::Rc d = F; d.h += 1;
            if (freeR(a) || freeR(b) || freeR(c) || freeR(d)) ++nonmax;
        }
        std::printf("CHECK overlaps=%d mer_notfree=%d mer_nonmax=%d\n", bad, notfree, nonmax);
    }
    if (std::getenv("CX_VERB")) {
        for (int i = 0; i < P.n; ++i) if (!B.placed[i])
            std::printf("UNPL %s cls=%d kind=%d area=%lld minside=%d bmask=%d\n", P.inst.blocks[i].name.c_str(), B.cls(i), P.kind[i],
                        (long long)B.barea[i], B.minside[i], P.bmask[i]);
        std::vector<cx::Rc> fr = B.M.f;
        std::sort(fr.begin(), fr.end(), [](const cx::Rc &a, const cx::Rc &b) { return a.area() > b.area(); });
        for (size_t k = 0; k < std::min<size_t>(8, fr.size()); ++k) std::printf("MER %d %d %d %d\n", fr[k].x, fr[k].y, fr[k].w, fr[k].h);
        int64_t freea = (int64_t)P.W * P.H; for (int i = 0; i < P.n; ++i) if (B.placed[i]) freea -= B.pos[i].area();
        std::printf("FREE %.3f%% of outline\n", 100.0 * freea / ((double)P.W * P.H));
        // union area of MERs with both sides >= the smallest min side of the unplaced blocks
        int ms = INT32_MAX; int64_t maxq = 0;
        for (int i = 0; i < P.n; ++i) if (!B.placed[i]) { ms = std::min(ms, B.minside[i]); maxq = std::max(maxq, B.barea[i]); }
        std::vector<cx::Rc> fat;
        for (auto &F : B.M.f) if (F.w >= ms && F.h >= ms) fat.push_back(F);
        std::vector<int> xs, ys;
        for (auto &F : fat) { xs.push_back(F.x); xs.push_back(F.x1()); ys.push_back(F.y); ys.push_back(F.y1()); }
        std::sort(xs.begin(), xs.end()); xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
        std::sort(ys.begin(), ys.end()); ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
        double ua = 0;
        for (size_t a = 0; a + 1 < xs.size(); ++a) for (size_t b = 0; b + 1 < ys.size(); ++b) {
            cx::Rc c{xs[a], ys[b], xs[a + 1] - xs[a], ys[b + 1] - ys[b]};
            for (auto &F : fat) if (cx::contains(F, c)) { ua += (double)c.area(); break; }
        }
        std::printf("FAT %.3f%% of outline (minside %d), %zu fat MERs, max queued area %lld\n", 100.0 * ua / ((double)P.W * P.H), ms, fat.size(), (long long)maxq);
    }
    int gate = 0, bsat = 0, msat = 0;
    double wl = 0;
    if (st.unplaced == 0) {
        std::vector<lab3::Rect> r = B.rects();
        GateInfo gi = gate_check(P, r);
        gate = gi.ok;
        bsat = gi.bsat;
        msat = gi.msat;
        if (!gi.ok) std::fprintf(stderr, "gate: %s\n", gi.why.c_str());
        if (wcons > 0) wl = (double)lab3::wirelength2(P.inst, r) / wcons;
        if (!outdir.empty()) write_output(P, r, outdir + "/" + tag + ".out");
    }
    std::printf("CONS %s placed=%d/%d unpl=%d qarea%%=%.3f maxmer=%dx%d(%.3f%%) nmer=%d/%d gate=%d wl=%.4f bsat=%d/%zu msat=%d/%zu "
                "msgp=%.1f msgreedy=%.1f msrep=%.1f rep=%d/%d\n",
                tag.c_str(), st.placed, P.n, st.unplaced, 100.0 * st.qarea / st.total_area, st.maxmer.w, st.maxmer.h,
                100.0 * st.maxmer.area() / st.total_area, st.nmer, st.max_nmer, gate, wl, bsat, P.inst.boundaries.size(), msat,
                P.inst.mibs.size(), ms_gp, st.ms_greedy, st.ms_repair, st.rep_acc, st.rep_tries);
    return 0;
}
