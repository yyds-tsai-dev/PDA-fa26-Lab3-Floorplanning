// test_legalize.cpp -- harness for legalize.h on fp's analytical seeds.
//
//   test_legalize <case.in> [--seed N=1] [--seeds K=8] [--out DIR] [--construction FILE]
//                 [--gp2] [--noshape] [--noflip] [--iters N] [--gp-iters N=25] [--gp-anchor X=0.05]
//
// Reproduces make_seeds() of fp.cpp (Rng(seed), K GlobalPlacer runs, the
// v-th with rand_axis = v > 0, shapes from sp_from_gp), legalizes each seed
// from its bisection-leaf centres and sp_from_gp shapes, then extracts an SP
// from the legal placement and packs it. Per seed it prints one line:
//   SEED v ctr=<WL ctr/cons> legal=<0|1> lp=<WL/cons> asap=<WL/cons> sp=<0|1> spwl=<WL/cons>
//        iters= flips= shapes= exc0= ms= mslp= bsat= msat=
// and, with --out, writes <DIR>/<tag>_s<v>_lp.out and _sp.out for the verifier.
// It also prints the fp-parity line "PARITY best_seed pen=.. wl2=.." (what fp
// prints as STATS best_seed for the same seed).
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "legalize.h"

using namespace fp;

int main(int argc, char **argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s case.in [options]\n", argv[0]); return 2; }
    std::string in = argv[1], outdir, cons;
    uint64_t seed = 1;
    int nseeds = 8, gp_iters = 25;
    double gp_anchor = 0.05;
    LegalizeOpts lo;
    bool leaf = false, qp = false, gp2 = false;   // gp2: obstacle-aware bisection + fit-first shapes   // qp: legalize the last quadratic positions G.x/G.y instead of leaf centres
    std::string centres;
    double jitter = 0;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto val = [&]() -> const char * { return i + 1 < argc ? argv[++i] : ""; };
        if (a == "--seed") seed = std::strtoull(val(), nullptr, 10);
        else if (a == "--seeds") nseeds = std::atoi(val());
        else if (a == "--out") outdir = val();
        else if (a == "--construction") cons = val();
        else if (a == "--noshape") lo.shape = false;
        else if (a == "--noflip") lo.flip = false;
        else if (a == "--noreloc") lo.reloc = false;
        else if (a == "--nogflip") lo.gflip = false;
        else if (a == "--notele") lo.teleport = false;
        else if (a == "--spread") lo.spread_iters = std::atoi(val());
        else if (a == "--gk") lo.gflip_k = std::atoi(val());
        else if (a == "--leaf") leaf = true;
        else if (a == "--qp") qp = true;
        else if (a == "--gp2") gp2 = true;
        else if (a == "--centres") centres = val();
        else if (a == "--jitter") jitter = std::atof(val());
        else if (a == "--iters") lo.max_iter = std::atoi(val());
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
    build_forbidden(P);   // fp strong mode builds it before make_seeds (does not touch rng)
    std::vector<lab3::Rect> cr;
    double wcons = 0;
    if (read_placement(P, cons, cr, err)) wcons = (double)lab3::wirelength2(P.inst, cr);
    else { std::fprintf(stderr, "construction: %s\n", err.c_str()); }
    auto rel = [&](double wl2) { return wcons > 0 ? wl2 / wcons : 0.0; };

    if (!centres.empty()) {   // oracle: legalize a given placement's centres and shapes
        std::vector<lab3::Rect> r;
        if (!read_placement(P, centres, r, err)) { std::fprintf(stderr, "centres: %s\n", err.c_str()); return 2; }
        std::vector<double> cx(P.n), cy(P.n);
        Rng jr(seed);
        for (int i = 0; i < P.n; ++i) {
            cx[i] = r[i].x + r[i].w / 2.0 + jitter * r[i].w * (jr.uni() - 0.5);
            cy[i] = r[i].y + r[i].h / 2.0 + jitter * r[i].h * (jr.uni() - 0.5);
        }
        std::vector<int> gidx(P.gmem.size(), 0);
        for (int g = 0; g < (int)P.gmem.size(); ++g) {
            const auto &ws = P.gshape[g].widths;
            int w = (int)r[P.gmem[g][0]].w, best = 0;
            for (int i = 0; i < (int)ws.size(); ++i) if (std::abs(ws[i] - w) < std::abs(ws[best] - w)) best = i;
            gidx[g] = best;
        }
        std::vector<lab3::Rect> rects;
        LegalizeStats st;
        bool ok = legalize(P, cx, cy, gidx, rects, nullptr, st, lo);
        GateInfo gi = ok ? gate_check(P, rects) : GateInfo();
        if (ok && !outdir.empty()) write_output(P, rects, outdir + "/" + tag + "_oracle_lp.out");
        std::printf("ORACLE in=%.4f legal=%d lp=%.4f asap=%.4f lpok=%d iters=%d flips=%d shapes=%d exc0=%lld ms=%.1f bsat=%d msat=%d\n",
                    rel((double)lab3::wirelength2(P.inst, r)), ok && gi.ok, ok ? rel((double)st.wl2) : 0.0,
                    ok ? rel((double)st.wl2_asap) : 0.0, st.lp_ok, st.iters, st.flips, st.shapes,
                    (long long)st.excess0, st.ms, gi.bsat, gi.msat);
        return 0;
    }
    Rng rng(seed);
    Packer K;
    K.init(P.n);
    std::vector<std::pair<int64_t, int64_t>> par;   // fp ranking parity
    for (int v = 0; v < nseeds; ++v) {
        auto t0 = std::chrono::steady_clock::now();
        GlobalPlacer G(P, rng);
        G.rand_axis = v > 0;
        G.obstacle = gp2;
        G.run(gp_iters, gp_anchor);
        double ms_gp = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        State s;
        if (!sp_from_gp(P, s, G, true, true, err)) { std::printf("SEED %d sp_from_gp failed: %s\n", v, err.c_str()); continue; }
        // fp parity: shape descent + pack of the seed SP (as make_seeds)
        {
            State t = s;
            if (invariant_ok(P, t)) {
                shape_descent(P, t, K);
                K.pack(P, t);
                par.push_back({K.ovfx + K.ovfy + K.viol, K.wl2});
            }
        }
        if (gp2) {
            gp_fit_shapes(P, G, s.gidx);
            for (int g = 0; g < (int)s.gidx.size(); ++g) if (G.sgidx[g] >= 0) s.gidx[g] = G.sgidx[g];
            init_dims(P, s);
            if (std::getenv("LG_LEAFDIAG")) std::printf("SIZEDIAG %d free_regions=%d unfit=%d\n", v, G.ob_nfree, G.ob_nmiss);
        }
        // GP centres WL (leaf centres, boundary snapped, sp_from_gp shapes).
        std::vector<double> cx(qp ? G.x : gp2 ? G.sx : G.tx), cy(qp ? G.y : gp2 ? G.sy : G.ty);
        double wctr;
        {
            double s2 = 0;
            std::vector<double> X(P.n), Y(P.n);
            for (int i = 0; i < P.n; ++i) {
                X[i] = cx[i]; Y[i] = cy[i];
                if (P.kind[i] == K_PRE) { X[i] = P.px[i] + P.fw[i] / 2.0; Y[i] = P.py[i] + P.fh[i] / 2.0; }
                int k = P.bmask[i];
                if (k & lab3::EDGE_LEFT) X[i] = s.w[i] / 2.0;
                if (k & lab3::EDGE_RIGHT) X[i] = P.W - s.w[i] / 2.0;
                if (k & lab3::EDGE_BOTTOM) Y[i] = s.h[i] / 2.0;
                if (k & lab3::EDGE_TOP) Y[i] = P.H - s.h[i] / 2.0;
            }
            for (const BBConn &c : P.bb) s2 += c.w * (std::fabs(X[c.a] - X[c.b]) + std::fabs(Y[c.a] - Y[c.b]));
            for (const BTConn &c : P.bt) s2 += c.w * (std::fabs(X[c.a] - c.tx2 / 2.0) + std::fabs(Y[c.a] - c.ty2 / 2.0));
            wctr = 2 * s2;
        }
        if (std::getenv("LG_LEAFDIAG")) {
            int nover = 0, npre = 0, nmov = 0, nfix = 0, nsmall = 0;
            double oh = 0, preov = 0;
            for (int i = 0; i < P.n; ++i) {
                if (P.kind[i] == K_PRE) continue;
                ++nmov;
                double ow = std::max(0.0, s.w[i] - G.leafw[i]), ohh = std::max(0.0, s.h[i] - G.leafh[i]);
                if (ow > 0 || ohh > 0) { ++nover; oh += ow + ohh; nfix += P.kind[i] == K_FIXED; nsmall += G.leafw[i] * G.leafh[i] < (double)s.w[i] * s.h[i]; }
                double x0 = G.tx[i] - G.leafw[i] / 2, x1 = G.tx[i] + G.leafw[i] / 2;
                double y0 = G.ty[i] - G.leafh[i] / 2, y1 = G.ty[i] + G.leafh[i] / 2;
                double a = 0;
                for (int p : P.pre) {
                    double ox = std::min(x1, (double)P.px[p] + P.fw[p]) - std::max(x0, (double)P.px[p]);
                    double oy = std::min(y1, (double)P.py[p] + P.fh[p]) - std::max(y0, (double)P.py[p]);
                    if (ox > 0 && oy > 0) a += ox * oy;
                }
                if (a > 0) { ++npre; preov += a / (G.leafw[i] * G.leafh[i]); }
            }
            std::printf("LEAFDIAG %d movable=%d overhang_blocks=%d overhang_sum=%.0f leaves_on_pre=%d mean_pre_frac=%.2f over_fixed=%d over_smallleaf=%d\n", v, nmov,
                        nover, oh, npre, npre ? preov / npre : 0.0, nfix, nsmall);
        }
        std::vector<int> gidx = s.gidx;
        std::vector<lab3::Rect> rects;
        mcf::ConstraintGraphs g;
        LegalizeStats st;
        bool ok = legalize(P, cx, cy, gidx, rects, &g, st, lo, leaf ? &G.leafw : nullptr, leaf ? &G.leafh : nullptr);
        int spok = 0, bsat = 0, msat = 0, lpgate = 0;
        double spwl = 0, sppol = 0;
        double ms_sp = 0;
        if (ok) {
            GateInfo gi = gate_check(P, rects);
            lpgate = gi.ok;
            bsat = gi.bsat; msat = gi.msat;
            if (!outdir.empty()) write_output(P, rects, outdir + "/" + tag + "_s" + std::to_string(v) + "_lp.out");
            // SP from the legal placement, pack + flush.
            auto t1 = std::chrono::steady_clock::now();
            std::vector<XRect> xr(P.n);
            std::vector<int64_t> kp(P.n), km(P.n);
            for (int i = 0; i < P.n; ++i) {
                xr[i] = {rects[i].x, rects[i].y, rects[i].w, rects[i].h};
                int64_t X2 = 2 * rects[i].x + rects[i].w, Y2 = 2 * rects[i].y + rects[i].h;
                kp[i] = X2 - Y2; km[i] = X2 + Y2;
            }
            State t;
            t.gidx = gidx;
            init_dims(P, t);
            if (extract_sp(xr, kp, km, t.gp, t.gm, err)) {
                t.pp.assign(P.n, 0); t.pm.assign(P.n, 0); t.sync_pos();
                K.pack(P, t);
                std::vector<lab3::Rect> pr = to_rects(K, t);
                GateInfo g2 = gate_check(P, pr);
                spok = g2.ok && K.feasible() && invariant_ok(P, t);
                spwl = rel((double)K.wl2);
                // LP polish of the packed SP on its own graphs (as fp's on_incumbent)
                mcf::ConstraintGraphs sg;
                sp_graphs(t.gp, t.pm, sg);
                mcf::PolishInfo pi;
                std::vector<lab3::Rect> pp2 = mcf::polish_topology(P.inst, pr, sg, &pi);
                sppol = rel((double)pi.wl2_out);
                if (!outdir.empty()) write_output(P, pr, outdir + "/" + tag + "_s" + std::to_string(v) + "_sp.out");
            }
            ms_sp = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count();
        }
        std::printf("SEED %d ctr=%.4f legal=%d lp=%.4f asap=%.4f lpok=%d sp=%d spwl=%.4f sppol=%.4f iters=%d flips=%d relocs=%d tele=%d "
                    "shapes=%d exc0=%lld exc=%lld ms=%.1f mslp=%.1f mssp=%.1f msgp=%.1f bsat=%d/%zu msat=%d/%zu why=%s\n",
                    v, rel(wctr), ok && lpgate, ok ? rel((double)st.wl2) : 0.0, ok ? rel((double)st.wl2_asap) : 0.0,
                    st.lp_ok, spok, spwl, sppol, st.iters, st.flips, st.relocs, st.teleports, st.shapes, (long long)st.excess0,
                    (long long)st.excess, st.ms, st.ms_lp, ms_sp, ms_gp, bsat, P.inst.boundaries.size(), msat,
                    P.inst.mibs.size(), st.why.empty() ? "-" : st.why.c_str());
        std::fflush(stdout);
    }
    if (!par.empty()) {
        int64_t wmin = INT64_MAX;
        for (auto &q : par) wmin = std::min(wmin, q.second);
        double bc = 1e300;
        size_t bi = 0;
        for (size_t i = 0; i < par.size(); ++i) {
            double c = (double)par[i].second / std::max<int64_t>(1, wmin) + 10.0 * par[i].first / (P.W + P.H);
            if (c < bc) { bc = c; bi = i; }
        }
        std::printf("PARITY best_seed pen=%lld wl2=%lld\n", (long long)par[bi].first, (long long)par[bi].second);
    }
    return 0;
}
