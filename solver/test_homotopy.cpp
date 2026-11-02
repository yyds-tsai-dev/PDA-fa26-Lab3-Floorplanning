// test_homotopy.cpp -- feasibility-homotopy compressibility diagnostic
// (docs/research/adhd-solver-round3.md, "Deepened 2. Feasibility homotopy").
//
//   test_homotopy <true.in> <relaxed.out> [--s0 0.92] [--budget 40] [--corr 1.0]
//                 [--seed N=1] [--out FILE] [--construction FILE] [--verbose]
//
// Loads the TRUE instance and a placement that is legal for the relaxed
// instance (movable areas scaled by s0), extracts its SP (keys X2-Y2 / X2+Y2
// on doubled centroids) and per soft group the TRUE shape index nearest to the
// relaxed aspect ratio. scale_dims(P, S, s): fixed = max(1, floor(true*sqrt s)),
// soft = true shape at gidx * sqrt s (rounded), preplaced untouched; the Packer
// reads the dims from State.
//
// feasible(s): pack inside the outline, preplaced pins satisfied, invariant_ok,
// forb_count == 0 (true forbidden table); at s = 1 also gate_check.
// s is searched on the grid 0.001 (rounding makes feasibility only roughly
// monotone in s), lower bound 0.80.
//
//  (i)  pack alone: largest feasible grid s in [0.80, 1] (also without the
//       forb_count term: s_geom).
//  (ii) predictor-corrector: s = bisect; while s < 1 and t < budget: corrector
//       SA (<= corr seconds) at target min(1, s + 0.005) with cost
//       wl2/wl_cons + 30*(ovf+viol)/(W+H) + forb, low T (1% increase accepted
//       with p = 0.1), stopping at the first state feasible at the target;
//       then bisect again. At s = 1: LP polish on the SP graphs, write --out.
//       On a stall: predecessor-chain diagnosis at min(1, s + 0.01).
//
// Prints one line:
//   HOMO case= s_start= forb0= sstar= sgeom= reached= s_final= steps= t1= wl_rel=
//        wl_s0pack= wl1pack= wl1pol= bsat= msat= trace= stall=
// --mode ramp (default) | absorb | direct, --gls T=20:
//   absorb: ramp, but stop when a step makes no progress at s >= 0.95 (or the
//           budget ends); pack the SP at s = 1, clamp into the outline,
//           gls_legalize (<= T s), LP polish (mcf::polish), write.
//   direct: no ramp; the relaxed rects re-centred with the TRUE nearest
//           shapes, clamped, gls_legalize (<= T s), LP polish, write.
//   Extra fields: path=ramp|gls|none gls_legal= gls_t= (GLS wall seconds).
// (trace: every 4th step s@t:packed WL ratio at s). --minwl: the corrector
// spends its whole budget and keeps the lowest-WL state feasible at the target.
// (WL values relative to the construction WL; wl_rel is the relaxed input).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "gls.h"

using namespace fp;

static void scale_dims(const Problem &P, State &S, double s) {
    S.w.assign(P.n, 0);
    S.h.assign(P.n, 0);
    const double r = std::sqrt(s);
    const bool one = s >= 1.0;
    for (int i = 0; i < P.n; ++i) {
        if (P.kind[i] == K_PRE) { S.w[i] = P.fw[i]; S.h[i] = P.fh[i]; }
        else if (P.kind[i] == K_FIXED) {
            S.w[i] = one ? P.fw[i] : std::max(1, (int)std::floor(P.fw[i] * r));
            S.h[i] = one ? P.fh[i] : std::max(1, (int)std::floor(P.fh[i] * r));
        }
    }
    for (int g = 0; g < (int)P.gmem.size(); ++g) {
        int w = P.gshape[g].widths[S.gidx[g]], h = P.gshape[g].heights[S.gidx[g]];
        if (!one) { w = std::max(1, (int)std::lround(w * r)); h = std::max(1, (int)std::lround(h * r)); }
        for (int b : P.gmem[g]) { S.w[b] = w; S.h[b] = h; }
    }
}

struct Eval {
    bool feas = false;
    int64_t ovf = 0, viol = 0, forb = 0, wl2 = 0;
    bool inv = false;
};

static Eval eval_at(const Problem &P, State &S, Packer &K, double s) {
    scale_dims(P, S, s);
    Eval e;
    e.inv = invariant_ok(P, S);
    e.forb = forb_count(P, S);
    K.pack(P, S);
    e.ovf = K.ovfx + K.ovfy;
    e.viol = K.viol;
    e.wl2 = K.wl2;
    e.feas = e.inv && e.forb == 0 && K.feasible();
    if (e.feas && s >= 1.0) e.feas = gate_check(P, to_rects(K, S)).ok;
    return e;
}

// Largest grid s in [lo, 1] (grid 0.001) with feasible(s); -1 if none.
// Bisection assumes monotonicity; the result is re-verified by a pack.
static double bisect(const Problem &P, State &S, Packer &K, double lo, bool use_forb) {
    auto ok = [&](int k) {
        Eval e = eval_at(P, S, K, k / 1000.0);
        return use_forb ? e.feas : (e.inv && K.feasible());
    };
    int a = (int)std::lround(lo * 1000), b = 1000;
    if (ok(b)) return 1.0;
    if (!ok(a)) {   // search downward to 0.80
        int lo2 = 800;
        if (a <= lo2 || !ok(lo2)) return -1;
        b = a; a = lo2;
    }
    while (b - a > 1) {   // ok(a), !ok(b)
        int m = (a + b) / 2;
        if (ok(m)) a = m; else b = m;
    }
    return a / 1000.0;
}

struct Corr {
    bool success = false;
    int64_t moves = 0;
};

// Corrector SA at target st; floor sf (state must stay feasible there for a
// fallback). Returns the first state feasible at st, or the lowest-cost
// visited state that is feasible at sf (else the input).
static Corr corrector(const Problem &P, State &S, Packer &K, double st, double sf, double budget, double wl_norm, Rng &rng,
                      bool minwl) {
    const int n = P.n;
    const double lam = 30.0, lamf = 1.0;
    auto cost_of = [&](const Eval &e) {
        return (double)e.wl2 / wl_norm + lam * (double)(e.ovf + e.viol) / (P.W + P.H) + lamf * (double)e.forb;
    };
    auto t0 = std::chrono::steady_clock::now();
    Corr R;
    State cur = S;
    Eval ec = eval_at(P, cur, K, st);
    if (ec.feas && !minwl) { S = cur; R.success = true; return R; }
    double cc = cost_of(ec);
    State bfeas;
    double bfc = 1e300;
    if (ec.feas) { bfeas = cur; bfc = cc; }
    const double T = 0.01 * cc / std::log(10.0);
    State best = S;
    double bc = 1e300;
    int ng = (int)P.gmem.size();
    State nx;
    while (true) {
        if ((R.moves & 255) == 0) {
            double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (el > budget) break;
        }
        ++R.moves;
        nx = cur;
        int mv = rng.below(ng > 0 ? 6 : 5);
        if (mv == 0) {          // swap in Γ+ / Γ- / both
            int a = rng.below(n), b = rng.below(n);
            if (a == b) continue;
            int which = rng.below(3);
            if (which != 1) std::swap(nx.gp[nx.pp[a]], nx.gp[nx.pp[b]]);
            if (which != 0) std::swap(nx.gm[nx.pm[a]], nx.gm[nx.pm[b]]);
        } else if (mv == 1 || mv == 2) {   // insert in one sequence (mv 2: both)
            int a = rng.below(n);
            auto ins = [&](std::vector<int> &q, std::vector<int> &pos) {
                int i = pos[a], j = rng.below(n);
                if (i == j) return;
                q.erase(q.begin() + i);
                q.insert(q.begin() + j, a);
            };
            if (mv == 2) { ins(nx.gp, nx.pp); ins(nx.gm, nx.pm); }
            else if (rng.below(2)) ins(nx.gp, nx.pp);
            else ins(nx.gm, nx.pm);
        } else if (mv == 3 || mv == 4) {   // adjacent swap in Γ+ or Γ-
            int i = rng.below(n - 1);
            if (mv == 3) std::swap(nx.gp[i], nx.gp[i + 1]);
            else std::swap(nx.gm[i], nx.gm[i + 1]);
        } else {                // shape index +-k for one group (MIB members move together)
            int g = rng.below(ng), sz = (int)P.gshape[g].widths.size();
            if (sz < 2) continue;
            int k = 1 + rng.below(std::max(1, std::min(5, sz / 10)));
            int v = nx.gidx[g] + (rng.below(2) ? k : -k);
            if (v < 0 || v >= sz) continue;
            nx.gidx[g] = v;
        }
        nx.sync_pos();
        if (!invariant_ok(P, nx)) continue;
        Eval e = eval_at(P, nx, K, st);
        if (e.feas && !minwl) { S = nx; R.success = true; return R; }
        double c = cost_of(e);
        if (e.feas && c < bfc) { bfeas = nx; bfc = c; }
        if (c <= cc || rng.uni() < std::exp(-(c - cc) / T)) {
            cur = nx; cc = c;
            if (c < bc) {
                State t = nx;
                if (eval_at(P, t, K, sf).feas) { best = nx; bc = c; }
            }
        }
    }
    if (bfc < 1e300) { S = bfeas; R.success = true; return R; }
    if (bc < 1e300) S = best;
    return R;
}

// O(n^2) longest-path pack at s with argmax predecessors; describes the
// blocking chains (overflow / preplaced-pin violation).
static std::string diagnose(const Problem &P, State S, double s) {
    scale_dims(P, S, s);
    const int n = P.n;
    std::vector<int64_t> x(n, 0), y(n, 0);
    std::vector<int> px(n, -1), py(n, -1);
    std::string out;
    int nviolx = 0, nvioly = 0;
    for (int i = 0; i < n; ++i) {   // x in Γ+ order
        int b = S.gp[i];
        int64_t q = 0; int arg = -1;
        for (int j = 0; j < i; ++j) {
            int a = S.gp[j];
            if (S.pm[a] < S.pm[b] && x[a] + S.w[a] > q) { q = x[a] + S.w[a]; arg = a; }
        }
        if (P.kind[b] == K_PRE) { x[b] = P.px[b]; if (q > P.px[b]) { ++nviolx; px[b] = arg; } }
        else { x[b] = q; px[b] = arg; }
    }
    for (int i = n - 1; i >= 0; --i) {
        int b = S.gp[i];
        int64_t q = 0; int arg = -1;
        for (int j = n - 1; j > i; --j) {
            int a = S.gp[j];
            if (S.pm[a] < S.pm[b] && y[a] + S.h[a] > q) { q = y[a] + S.h[a]; arg = a; }
        }
        if (P.kind[b] == K_PRE) { y[b] = P.py[b]; if (q > P.py[b]) { ++nvioly; py[b] = arg; } }
        else { y[b] = q; py[b] = arg; }
    }
    auto chain = [&](int b, bool ax) {
        // walk predecessors back to a preplaced pin or the origin
        std::string c;
        int cur = b, len = 0;
        bool pin = false;
        while (cur >= 0 && len < 60) {
            const char *nm = P.inst.blocks[cur].name.c_str();
            const char *k = P.kind[cur] == K_PRE ? "P" : P.kind[cur] == K_FIXED ? "F" : "S";
            c = std::string(nm) + k + (c.empty() ? "" : "<" + c);
            ++len;
            if (P.kind[cur] == K_PRE && cur != b) { pin = true; break; }
            cur = ax ? py[cur] : px[cur];
        }
        return std::string(pin ? "pin:" : "org:") + std::to_string(len) + ":" + c;
    };
    // preplaced violations
    int shown = 0;
    for (int p : P.pre) {
        int64_t qx = 0, qy = 0;
        if (px[p] >= 0) qx = x[px[p]] + S.w[px[p]] - P.px[p];
        if (py[p] >= 0) qy = y[py[p]] + S.h[py[p]] - P.py[p];
        if (qx > 0 && shown < 2) { out += " PRE-X(" + P.inst.blocks[p].name + ",+" + std::to_string(qx) + ")=" + chain(px[p], false) + ">" + P.inst.blocks[p].name; ++shown; }
        if (qy > 0 && shown < 2) { out += " PRE-Y(" + P.inst.blocks[p].name + ",+" + std::to_string(qy) + ")=" + chain(py[p], true) + ">" + P.inst.blocks[p].name; ++shown; }
    }
    int bx = -1, by = -1;
    int64_t mx = P.W, my = P.H;
    for (int b = 0; b < n; ++b) {
        if (x[b] + S.w[b] > mx) { mx = x[b] + S.w[b]; bx = b; }
        if (y[b] + S.h[b] > my) { my = y[b] + S.h[b]; by = b; }
    }
    if (bx >= 0) out += " OVF-X(+" + std::to_string(mx - P.W) + ")=" + chain(bx, false);
    if (by >= 0) out += " OVF-Y(+" + std::to_string(my - P.H) + ")=" + chain(by, true);
    int64_t f = forb_count(P, S);
    if (f) {
        out += " FORB=" + std::to_string(f);
        int k = 0;
        for (const auto &fb : P.forb) {
            int r = sp_rel(S.pp, S.pm, fb[0], fb[1]);
            if ((fb[2] >> r) & 1 && k++ < 3) {
                static const char *rn[4] = {"left-of", "right-of", "below", "above"};
                out += "(" + P.inst.blocks[fb[0]].name + " " + rn[r] + " " + P.inst.blocks[fb[1]].name + ")";
            }
        }
    }
    out += " npreviol=" + std::to_string(nviolx + nvioly);
    return out;
}

int main(int argc, char **argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: %s true.in relaxed.out [options]\n", argv[0]); return 2; }
    std::string in = argv[1], rel_out = argv[2], out, cons;
    double s0 = 0.92, budget = 40, corr_t = 1.0;
    uint64_t seed = 1;
    bool verbose = false, minwl = false;
    std::string mode = "ramp";
    double gls_t = 20;
    for (int i = 3; i < argc; ++i) {
        std::string a = argv[i];
        auto val = [&]() -> const char * { return i + 1 < argc ? argv[++i] : ""; };
        if (a == "--s0") s0 = std::atof(val());
        else if (a == "--budget") budget = std::atof(val());
        else if (a == "--corr") corr_t = std::atof(val());
        else if (a == "--seed") seed = std::strtoull(val(), nullptr, 10);
        else if (a == "--out") out = val();
        else if (a == "--construction") cons = val();
        else if (a == "--verbose") verbose = true;
        else if (a == "--minwl") minwl = true;
        else if (a == "--mode") mode = val();
        else if (a == "--gls") gls_t = std::atof(val());   // corrector: full budget, keep the lowest-WL state feasible at the target
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
    build_forbidden(P);
    std::vector<lab3::Rect> cr, rr;
    if (!read_placement(P, cons, cr, err)) { std::fprintf(stderr, "construction: %s\n", err.c_str()); return 2; }
    const double wcons = (double)lab3::wirelength2(P.inst, cr);
    if (!read_placement(P, rel_out, rr, err)) { std::fprintf(stderr, "relaxed: %s\n", err.c_str()); return 2; }
    const double wl_rel = (double)lab3::wirelength2(P.inst, rr) / wcons;

    // SP + nearest true shape (aspect ratio) from the relaxed placement.
    State S;
    S.gidx.assign(P.gmem.size(), 0);
    for (int g = 0; g < (int)P.gmem.size(); ++g) {
        const lab3::Rect &q = rr[P.gmem[g][0]];
        double ar = std::log((double)q.w / q.h), bd = 1e300;
        const ShapeRange &sr = P.gshape[g];
        for (int i = 0; i < (int)sr.widths.size(); ++i) {
            double d = std::fabs(std::log((double)sr.widths[i] / sr.heights[i]) - ar);
            if (d < bd) { bd = d; S.gidx[g] = i; }
        }
    }
    {
        std::vector<XRect> xr(P.n);
        std::vector<int64_t> kp(P.n), km(P.n);
        for (int i = 0; i < P.n; ++i) {
            xr[i] = {rr[i].x, rr[i].y, rr[i].w, rr[i].h};
            int64_t X2 = 2 * rr[i].x + rr[i].w, Y2 = 2 * rr[i].y + rr[i].h;
            kp[i] = X2 - Y2; km[i] = X2 + Y2;
        }
        if (!extract_sp(xr, kp, km, S.gp, S.gm, err)) { std::fprintf(stderr, "extract: %s\n", err.c_str()); return 2; }
        S.pp.assign(P.n, 0); S.pm.assign(P.n, 0); S.sync_pos();
    }
    Packer K;
    K.init(P.n);
    const int64_t forb0 = forb_count(P, S);
    const bool inv0 = invariant_ok(P, S);
    Eval e0 = eval_at(P, S, K, s0);
    const double wl_s0pack = (double)e0.wl2 / wcons;

    // (i) pack alone
    State T = S;
    const double sstar = bisect(P, T, K, s0, true);
    const double sgeom = bisect(P, T, K, s0, false);

    // (ii) predictor-corrector
    auto t0 = std::chrono::steady_clock::now();
    auto el = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    Rng rng(seed);
    double s = sstar;
    const double s_start = s;
    int steps = 0, nsucc = 0;
    int64_t moves = 0;
    std::string trace;
    if (mode != "ramp" && mode != "absorb" && mode != "direct") { std::fprintf(stderr, "bad mode\n"); return 2; }
    while (mode != "direct" && (s < 0 || s < 1.0) && el() < budget) {
        double st = s < 0 ? s0 : std::min(1.0, s + 0.005);
        double sf = s < 0 ? s0 : s;
        Corr c = corrector(P, S, K, st, sf, std::min(corr_t, budget - el()), wcons, rng, minwl);
        ++steps;
        nsucc += c.success;
        moves += c.moves;
        double sn = bisect(P, S, K, s < 0 ? s0 : s, true);
        if (verbose) std::fprintf(stderr, "step %d t=%.2f target=%.3f succ=%d moves=%lld s=%.3f -> %.3f\n", steps, el(), st, c.success,
                                  (long long)c.moves, s, sn);
        const bool stalled = !(sn > s);
        if (sn > s) s = sn;
        if (mode == "absorb" && stalled && s >= 0.95) break;
        if (steps <= 400 && (steps % 4 == 0 || s >= 1.0)) {
            Eval ew = eval_at(P, S, K, s < 0 ? s0 : s);
            char b[48]; std::snprintf(b, sizeof b, "%.3f@%.1f:%.3f,", s, el(), (double)ew.wl2 / wcons); trace += b;
        }
    }
    const bool reached = s >= 1.0;
    const double t1 = el();
    double wl1pack = 0, wl1pol = 0;
    int bsat = 0, msat = 0;
    std::string stall = "-", path = reached ? "ramp" : "none";
    int gls_legal = 0;
    double gls_secs = 0;
    if (!reached && mode != "ramp") {
        std::vector<lab3::Rect> r(P.n);
        std::vector<int> gidx = S.gidx;
        if (mode == "absorb") {
            Eval e1 = eval_at(P, S, K, 1.0);
            r = to_rects(K, S);
            wl1pack = (double)e1.wl2 / wcons;
        } else {   // direct: relaxed centres, true dims
            State t = S;
            scale_dims(P, t, 1.0);
            for (int i = 0; i < P.n; ++i) {
                if (P.kind[i] == K_PRE) { r[i] = lab3::Rect{P.px[i], P.py[i], P.fw[i], P.fh[i]}; continue; }
                int64_t cx2 = 2 * rr[i].x + rr[i].w, cy2 = 2 * rr[i].y + rr[i].h;
                r[i] = lab3::Rect{(cx2 - t.w[i]) / 2, (cy2 - t.h[i]) / 2, t.w[i], t.h[i]};
            }
        }
        for (int i = 0; i < P.n; ++i) {   // clamp into the outline
            r[i].x = std::max<int64_t>(0, std::min<int64_t>(r[i].x, P.W - r[i].w));
            r[i].y = std::max<int64_t>(0, std::min<int64_t>(r[i].y, P.H - r[i].h));
        }
        Rng grng(seed * 0x9E3779B97F4A7C15ULL + 12345);
        GlsStats gs;
        auto tg = std::chrono::steady_clock::now();
        bool ok = gls_legalize(P, r, gidx, gls_t, grng, gs);
        gls_secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - tg).count();
        GateInfo gi = gate_check(P, r);
        if (ok && gi.ok) {
            gls_legal = 1;
            path = "gls";
            wl1pack = (double)lab3::wirelength2(P.inst, r) / wcons;
            mcf::PolishInfo pi;
            std::vector<lab3::Rect> pol = mcf::polish(P.inst, r, &pi);
            GateInfo g2 = gate_check(P, pol);
            if (!g2.ok) { pol = r; g2 = gi; }
            wl1pol = (double)lab3::wirelength2(P.inst, pol) / wcons;
            bsat = g2.bsat; msat = g2.msat;
            if (!out.empty()) write_output(P, pol, out);
        }
    }
    if (reached) {
        Eval e1 = eval_at(P, S, K, 1.0);
        std::vector<lab3::Rect> pr = to_rects(K, S);
        wl1pack = (double)e1.wl2 / wcons;
        mcf::ConstraintGraphs sg;
        sp_graphs(S.gp, S.pm, sg);
        mcf::PolishInfo pi;
        std::vector<lab3::Rect> pol = mcf::polish_topology(P.inst, pr, sg, &pi);
        wl1pol = (double)pi.wl2_out / wcons;
        GateInfo gi = gate_check(P, pol);
        if (!gi.ok) { pol = pr; gi = gate_check(P, pr); wl1pol = wl1pack; }
        bsat = gi.bsat; msat = gi.msat;
        if (!out.empty()) write_output(P, pol, out);
    } else if (mode != "direct") {
        double sd = s < 0 ? s0 : std::min(1.0, s + 0.01);
        char b[32]; std::snprintf(b, sizeof b, "@%.3f", sd);
        stall = b + diagnose(P, S, sd);
    }
    std::printf("HOMO case=%s inv0=%d forb0=%lld s_start=%.3f sstar=%.3f sgeom=%.3f reached=%d s_final=%.3f steps=%d succ=%d moves=%lld "
                "t1=%.2f wl_rel=%.4f wl_s0pack=%.4f wl1pack=%.4f wl1pol=%.4f bsat=%d/%zu msat=%d/%zu path=%s gls_legal=%d gls_t=%.2f trace=%s stall=%s\n",
                tag.c_str(), inv0, (long long)forb0, s_start, sstar, sgeom, reached, s, steps, nsucc, (long long)moves, t1, wl_rel,
                wl_s0pack, wl1pack, wl1pol, bsat, P.inst.boundaries.size(), msat, P.inst.mibs.size(), path.c_str(), gls_legal, gls_secs, trace.empty() ? "-" : trace.c_str(),
                stall.c_str());
    return 0;
}
