// test_gls.cpp -- harness for gls.h (coordinate-space GLS overlap removal).
//
//   test_gls <case.in> [--seed N=1] [--start seed|ctr] [--time T=60] [--seeds K=8]
//            [--out DIR] [--construction FILE] [--verbose] [--wl SECONDS] [--wlw W] [--wlw-half S]
//   test_gls <case.in> --polish IN.out --polish-out OUT.out [--time T] [--seed N] [--wl-T X]
//            (post-legal polish: only the WL phase of gls_wl_polish on a legal
//            placement for T seconds; prints POL case= ok= wlin= wlout= gain= steps=)
//
// Start points (both from fp's make_seeds: K GlobalPlacer runs with Rng(seed),
// ranked as fp ranks them, best one taken):
//   seed : the seed SP after shape descent, packed by the Packer (may overflow
//          the outline; gls clamps it in);
//   ctr  : the GP bisection-leaf centres with the sp_from_gp shapes.
// Then gls_legalize with the remaining budget (T minus the seed time). On a
// legal result: extract an SP (extract.h), pack it (+ flush) and check gate,
// Packer feasibility and the SP invariant (inv=1: invariant_ok and no
// forbidden block-preplaced relation, forb_count == 0, fp's strong_inv).
// Prints one line:
//   GLS case= start= seed= legal= ttl= ttl_tot= ms_seed= steps= bumps= minov=
//       wl0= wl= sp= spwl= bsat=a/b msat=c/d inv=
// (wl0: at the first legal state, wl: after the optional --wl phase, spwl: the
// packed extracted SP; all relative to the construction WL). With --out writes
// <DIR>/<tag>_<start>_s<seed>_gls.out and _sp.out.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "gls.h"

using namespace fp;

int main(int argc, char **argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s case.in [options]\n", argv[0]); return 2; }
    std::string in = argv[1], outdir, cons, start = "seed", pol_in, pol_out;
    uint64_t seed = 1;
    int nseeds = 8;
    double T = 60;
    GlsOpts go;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto val = [&]() -> const char * { return i + 1 < argc ? argv[++i] : ""; };
        if (a == "--seed") seed = std::strtoull(val(), nullptr, 10);
        else if (a == "--seeds") nseeds = std::atoi(val());
        else if (a == "--start") start = val();
        else if (a == "--time") T = std::atof(val());
        else if (a == "--out") outdir = val();
        else if (a == "--construction") cons = val();
        else if (a == "--verbose") go.verbose = true;
        else if (a == "--gls-every") go.gls_every = std::atoll(val());
        else if (a == "--gls-inc") go.gls_inc = std::atof(val());
        else if (a == "--period") go.t_period = std::atof(val());
        else if (a == "--wl") go.wl_time = std::atof(val());
        else if (a == "--wlw") go.wlw = std::atof(val());
        else if (a == "--wlw-half") go.wlw_half = std::atof(val());
        else if (a == "--wl-T") go.wl_T = std::atof(val());
        else if (a == "--polish") pol_in = val();
        else if (a == "--polish-out") pol_out = val();
        else { std::fprintf(stderr, "bad option %s\n", a.c_str()); return 2; }
    }
    if (start != "seed" && start != "ctr") { std::fprintf(stderr, "bad --start\n"); return 2; }
    auto t0 = std::chrono::steady_clock::now();
    auto el = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
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
    build_forbidden(P);
    if (!pol_in.empty()) {   // post-legal WL polish only
        std::vector<lab3::Rect> r;
        if (!read_placement(P, pol_in, r, err)) { std::fprintf(stderr, "polish input: %s\n", err.c_str()); return 2; }
        std::vector<int> g(P.gmem.size(), 0);
        bool exact = true;
        for (int k = 0; k < (int)P.gmem.size(); ++k) {
            const auto &ws = P.gshape[k].widths;
            int w = (int)r[P.gmem[k][0]].w;
            auto it = std::lower_bound(ws.begin(), ws.end(), w);
            if (it == ws.end() || *it != w) { exact = false; break; }
            g[k] = (int)(it - ws.begin());
            for (int j : P.gmem[k]) if (r[j].w != w) exact = false;
        }
        double wlin = (double)lab3::wirelength2(P.inst, r);
        GlsStats st;
        Rng prng(seed * 0x9E3779B97F4A7C15ULL + 777);
        bool ok = exact && gls_wl_polish(P, r, g, T, prng, st, go);
        GateInfo gi = gate_check(P, r);
        ok = ok && gi.ok;
        if (ok && !pol_out.empty()) write_output(P, r, pol_out);
        double wlout = ok ? (double)st.wl2_out : wlin;
        std::printf("POL case=%s ok=%d exact=%d wlin=%.1f wlstart=%.1f wlout=%.1f gain=%.5f steps=%lld bsat=%d/%zu msat=%d/%zu\n", tag.c_str(),
                    ok ? 1 : 0, exact ? 1 : 0, wlin / 2, st.wl2_legal / 2.0, wlout / 2, 1.0 - wlout / wlin, (long long)st.steps, gi.bsat,
                    P.inst.boundaries.size(), gi.msat, P.inst.mibs.size());
        return ok ? 0 : 1;
    }
    std::vector<lab3::Rect> cr;
    double wcons = 0;
    if (read_placement(P, cons, cr, err)) wcons = (double)lab3::wirelength2(P.inst, cr);
    auto rel = [&](double wl2) { return wcons > 0 ? wl2 / wcons : 0.0; };

    // make_seeds (fp.cpp) with the best-ranked seed kept, plus its GP centres.
    Rng rng(seed);
    Packer K;
    K.init(P.n);
    double bestc = 1e300;
    std::vector<lab3::Rect> sr;
    std::vector<int> sg;
    std::vector<std::pair<int64_t, int64_t>> sc;
    std::vector<State> seeds;
    std::vector<std::vector<double>> ctx, cty;
    std::vector<std::vector<int>> cg;
    for (int v = 0; v < nseeds; ++v) {
        GlobalPlacer G(P, rng);
        G.rand_axis = v > 0;
        G.run(25, 0.05);
        State s;
        if (!sp_from_gp(P, s, G, true, true, err)) continue;
        if (!invariant_ok(P, s)) continue;
        std::vector<int> g0 = s.gidx;
        shape_descent(P, s, K);
        K.pack(P, s);
        sc.push_back({K.ovfx + K.ovfy + K.viol, K.wl2});
        seeds.push_back(s);
        ctx.push_back(G.tx); cty.push_back(G.ty); cg.push_back(g0);
    }
    if (seeds.empty()) { std::printf("GLS case=%s no seeds\n", tag.c_str()); return 1; }
    int64_t wmin = INT64_MAX;
    for (auto &q : sc) wmin = std::min(wmin, q.second);
    size_t bi = 0;
    for (size_t i = 0; i < seeds.size(); ++i) {
        double c = (double)sc[i].second / std::max<int64_t>(1, wmin) + 10.0 * sc[i].first / (P.W + P.H);
        if (c < bestc) { bestc = c; bi = i; }
    }
    if (start == "seed") {
        K.pack(P, seeds[bi]);
        sr = to_rects(K, seeds[bi]);
        sg = seeds[bi].gidx;
    } else {
        sg = cg[bi];
        State t;
        t.gidx = sg;
        init_dims(P, t);
        sr.assign(P.n, lab3::Rect());
        for (int i = 0; i < P.n; ++i)
            sr[i] = lab3::Rect{(int64_t)std::lround(ctx[bi][i] - t.w[i] / 2.0), (int64_t)std::lround(cty[bi][i] - t.h[i] / 2.0), t.w[i], t.h[i]};
    }
    double ms_seed = el() * 1000;

    Rng grng(seed * 0x9E3779B97F4A7C15ULL + 12345);
    GlsStats st;
    std::vector<lab3::Rect> r = sr;
    std::vector<int> gidx = sg;
    bool legal = gls_legalize(P, r, gidx, std::max(0.0, T - el()), grng, st, go);
    double ttl_tot = legal ? ms_seed / 1000 + st.ttl : -1;
    GateInfo gi = gate_check(P, r);
    legal = legal && gi.ok;
    int spok = 0, inv = 0;
    double spwl = 0;
    if (legal) {
        if (!outdir.empty()) write_output(P, r, outdir + "/" + tag + "_" + start + "_s" + std::to_string(seed) + "_gls.out");
        std::vector<XRect> xr(P.n);
        std::vector<int64_t> kp(P.n), km(P.n);
        for (int i = 0; i < P.n; ++i) {
            xr[i] = {r[i].x, r[i].y, r[i].w, r[i].h};
            int64_t X2 = 2 * r[i].x + r[i].w, Y2 = 2 * r[i].y + r[i].h;
            kp[i] = X2 - Y2; km[i] = X2 + Y2;
        }
        State t;
        t.gidx = gidx;
        init_dims(P, t);
        if (extract_sp(xr, kp, km, t.gp, t.gm, err)) {
            t.pp.assign(P.n, 0); t.pm.assign(P.n, 0); t.sync_pos();
            inv = invariant_ok(P, t) && forb_count(P, t) == 0;
            K.pack(P, t);
            std::vector<lab3::Rect> pr = to_rects(K, t);
            GateInfo g2 = gate_check(P, pr);
            spok = g2.ok && K.feasible() && inv && g2.bsat == (int)P.inst.boundaries.size() && g2.msat == (int)P.inst.mibs.size();
            spwl = rel((double)K.wl2);
            if (!outdir.empty()) write_output(P, pr, outdir + "/" + tag + "_" + start + "_s" + std::to_string(seed) + "_sp.out");
        } else {
            std::fprintf(stderr, "extract: %s\n", err.c_str());
        }
    }
    std::printf("GLS case=%s start=%s seed=%llu legal=%d ttl=%.3f ttl_tot=%.3f ms_seed=%.1f steps=%lld acc=%lld bumps=%lld minov=%.5f "
                "wl0=%.4f wl=%.4f sp=%d spwl=%.4f bsat=%d/%zu msat=%d/%zu inv=%d why=%s\n",
                tag.c_str(), start.c_str(), (unsigned long long)seed, legal ? 1 : 0, legal ? st.ttl : -1.0, ttl_tot, ms_seed,
                (long long)st.steps, (long long)st.accepted, (long long)st.bumps, st.min_overlap, legal ? rel((double)st.wl2_legal) : 0.0, legal ? rel((double)st.wl2_out) : 0.0,
                spok, spwl, gi.bsat, P.inst.boundaries.size(), gi.msat, P.inst.mibs.size(), inv, gi.ok ? "-" : gi.why.c_str());
    return legal ? 0 : 1;
}
