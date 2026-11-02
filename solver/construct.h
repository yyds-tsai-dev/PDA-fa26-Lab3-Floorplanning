// construct.h -- legal-first constructive greedy (go/no-go prototype, not used
// by solver/fp). See docs/research/adhd-solver-round2.md, shortlist item 2.
//
// State = legal partial placement + queue of unplaced blocks. Free space is the
// set of maximal empty rectangles (MERs) inside the outline, kept by the
// MaxRects split rule (split every MER that intersects a placed rect into up
// to four remnants, drop remnants contained in another MER).
//
// Order: preplaced (rocks), MIB groups (one shared shape, members placed one
// after another), corner/edge boundary blocks, other fixed blocks, soft blocks;
// within a class, largest area first. Candidates: a block anchored at one of
// the four corners of a MER (legal by construction). Soft shapes snap one side
// to the MER (MER side, or MER side minus the min side of a following block)
// and take the other side from the area. Score (lexicographic): dead area
// (strips thinner than the smallest min side of the still-unplaced blocks),
// then walls touched (more is better), then distance to the analytical target
// (or bottom-left in blind mode). Blocks with no candidate go to the queue.
//
// Repair (optional): evict 2-8 soft blocks in a window around a queued block's
// target / the largest residual pocket, re-pour queued + evicted blocks, and
// keep the result only if the queued area drops.
#ifndef FP_CONSTRUCT_H
#define FP_CONSTRUCT_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <vector>

#include "gp.h"

namespace fp {
namespace cx {

struct Rc {
    int x = 0, y = 0, w = 0, h = 0;
    int x1() const { return x + w; }
    int y1() const { return y + h; }
    int64_t area() const { return (int64_t)w * h; }
};
inline bool ovl(const Rc &a, const Rc &b) {
    return std::max(a.x, b.x) < std::min(a.x1(), b.x1()) && std::max(a.y, b.y) < std::min(a.y1(), b.y1());
}
inline bool contains(const Rc &a, const Rc &b) {
    return a.x <= b.x && a.y <= b.y && a.x1() >= b.x1() && a.y1() >= b.y1();
}

struct MerSet {
    std::vector<Rc> f;
    void init(int W, int H) { f.assign(1, Rc{0, 0, W, H}); }
    void place(const Rc &r) {
        std::vector<Rc> keep, fresh;
        keep.reserve(f.size() + 8);
        for (const Rc &q : f) {
            if (!ovl(q, r)) { keep.push_back(q); continue; }
            if (r.x > q.x) fresh.push_back({q.x, q.y, r.x - q.x, q.h});
            if (r.x1() < q.x1()) fresh.push_back({r.x1(), q.y, q.x1() - r.x1(), q.h});
            if (r.y > q.y) fresh.push_back({q.x, q.y, q.w, r.y - q.y});
            if (r.y1() < q.y1()) fresh.push_back({q.x, r.y1(), q.w, q.y1() - r.y1()});
        }
        size_t nk = keep.size();
        for (size_t i = 0; i < fresh.size(); ++i) {
            bool dom = false;
            for (size_t k = 0; k < nk && !dom; ++k) dom = contains(keep[k], fresh[i]);
            for (size_t j = 0; j < fresh.size() && !dom; ++j) {
                if (j == i || !contains(fresh[j], fresh[i])) continue;
                dom = !contains(fresh[i], fresh[j]) || j < i;   // strict container, or earlier duplicate
            }
            if (!dom) keep.push_back(fresh[i]);
        }
        f.swap(keep);
    }
};

struct Opts {
    bool target = true;      // distance to analytical target as last key (else bottom-left)
    int score = 0;           // 0: dead, walls, dist   1: dead, dist (walls ignored)
    int bestfit = 0;         // >0: best-fit choice among the next K plain (non-MIB, non-boundary) blocks
    int evict_max = 8;       // repair evicts 2..evict_max blocks
    bool evict_fixed = false;
    bool plateau = false;    // repair also accepts equal queued area (sideways moves)
    int max_fails = 400;
    int order = 0;           // 0: corners, MIB, edges, fixed, soft   1: all boundary before MIB
    bool repair = false;
    double repair_ms = 1500; // wall budget for repair
    uint64_t seed = 1;
};

struct Stats {
    int placed = 0, unplaced = 0;
    int64_t qarea = 0, total_area = 0;
    Rc maxmer;
    int nmer = 0, max_nmer = 0;
    int rep_tries = 0, rep_acc = 0;
    double ms_greedy = 0, ms_repair = 0;
};

struct Builder {
    const Problem &P;
    const std::vector<double> *tx = nullptr, *ty = nullptr;
    Opts o;
    int n;
    std::vector<Rc> pos;
    std::vector<char> placed;
    MerSet M;
    std::vector<int> gsw, gsh;      // per shape group: chosen shape (-1 = free)
    std::vector<int> minside;       // per block
    std::vector<int64_t> barea;     // per block (soft: A)
    std::vector<int> rest;          // blocks not yet processed (for dead threshold)
    int max_nmer = 0;

    Builder(const Problem &p, const Opts &op) : P(p), o(op), n(p.n) {
        pos.assign(n, Rc{});
        placed.assign(n, 0);
        gsw.assign(P.gmem.size(), -1);
        gsh.assign(P.gmem.size(), -1);
        minside.assign(n, 0);
        barea.assign(n, 0);
        for (int i = 0; i < n; ++i) {
            if (P.kind[i] == K_SOFT) {
                const ShapeRange &s = P.gshape[P.group[i]];
                minside[i] = s.widths.front();
                barea[i] = s.A;
            } else {
                minside[i] = std::min(P.fw[i], P.fh[i]);
                barea[i] = (int64_t)P.fw[i] * P.fh[i];
            }
        }
    }

    void reset_mers() {
        M.init(P.W, P.H);
        for (int i = 0; i < n; ++i) if (placed[i]) M.place(pos[i]);
    }
    void commit(int b, const Rc &r) {
        pos[b] = r;
        placed[b] = 1;
        M.place(r);
        max_nmer = std::max(max_nmer, (int)M.f.size());
        if (P.kind[b] == K_SOFT && gsw[P.group[b]] < 0) { gsw[P.group[b]] = r.w; gsh[P.group[b]] = r.h; }
    }

    // ---- scoring
    struct Sc {
        int64_t dead = INT64_MAX;
        int walls = 0;
        double dist = 1e300;
    };
    bool better(const Sc &a, const Sc &b) const {   // a strictly better than b
        if (a.dead != b.dead) return a.dead < b.dead;
        if (o.score == 0 && a.walls != b.walls) return a.walls > b.walls;
        return a.dist < b.dist;
    }
    int64_t dead_of(const Rc &c, int thr) const {
        int64_t d = 0;
        for (const Rc &q : M.f) {
            if (!ovl(q, c)) continue;
            int t;
            if ((t = c.x - q.x) > 0 && t < thr) d += (int64_t)t * q.h;
            if ((t = q.x1() - c.x1()) > 0 && t < thr) d += (int64_t)t * q.h;
            if ((t = c.y - q.y) > 0 && t < thr) d += (int64_t)t * q.w;
            if ((t = q.y1() - c.y1()) > 0 && t < thr) d += (int64_t)t * q.w;
        }
        return d;
    }
    int walls_of(const Rc &c) const {
        bool L = c.x == 0, R = c.x1() == P.W, B = c.y == 0, T = c.y1() == P.H;
        for (int i = 0; i < n && !(L && R && B && T); ++i) {
            if (!placed[i]) continue;
            const Rc &q = pos[i];
            bool yo = std::max(q.y, c.y) < std::min(q.y1(), c.y1());
            bool xo = std::max(q.x, c.x) < std::min(q.x1(), c.x1());
            if (yo && q.x1() == c.x) L = true;
            if (yo && q.x == c.x1()) R = true;
            if (xo && q.y1() == c.y) B = true;
            if (xo && q.y == c.y1()) T = true;
        }
        return (int)L + R + B + T;
    }
    double dist_of(int b, const Rc &c) const {
        double cxx = c.x + c.w / 2.0, cyy = c.y + c.h / 2.0;
        if (o.target && tx) return std::fabs(cxx - (*tx)[b]) + std::fabs(cyy - (*ty)[b]);
        return (double)c.y * (P.W + 1.0) + c.x;
    }

    // ---- shapes
    // Soft shapes for block b in a pocket F (all with w <= F.w, h <= F.h).
    void soft_shapes(int b, const Rc &F, const std::vector<int> &ms, std::vector<std::pair<int, int>> &out) const {
        const ShapeRange &s = P.gshape[P.group[b]];
        const std::vector<int> &W = s.widths;
        auto add = [&](int w, int h) {
            if (w <= F.w && h <= F.h) out.push_back({w, h});
        };
        // widest legal w <= tw with hmin(w) <= hcap; push (w,hmin) and (w, min(hmax, hcap)).
        auto snap = [&](int tw, int hcap, bool swap) {
            auto it = std::upper_bound(W.begin(), W.end(), tw);
            int k = (int)(it - W.begin()) - 1;
            for (int steps = 0; k >= 0 && steps < 64; --k, ++steps) {
                int64_t w = W[k], lo = shape_hmin(s.whmin, s.R, w), hi = shape_hmax(s.whmax, s.R, w);
                if (lo > hcap) continue;
                int h1 = (int)lo, h2 = (int)std::min<int64_t>(hi, hcap);
                if (swap) { add(h1, (int)w); if (h2 != h1) add(h2, (int)w); }
                else { add((int)w, h1); if (h2 != h1) add((int)w, h2); }
                return;
            }
        };
        std::vector<int> tws = {F.w}, ths = {F.h};
        for (int m : ms) { if (F.w - m > 0) tws.push_back(F.w - m); if (F.h - m > 0) ths.push_back(F.h - m); }
        for (int tw : tws) snap(tw, F.h, false);
        for (int th : ths) snap(th, F.w, true);   // symmetric shape set: snap the height
        // near-square
        int sq = (int)std::sqrt((double)s.A);
        auto it = std::lower_bound(W.begin(), W.end(), sq);
        if (it != W.end()) {
            int w = *it;
            add(w, (int)shape_hmin(s.whmin, s.R, w));
        }
    }

    struct Cand {
        Rc r;
        Sc s;
        bool ok = false;
    };
    // Best candidate for block b. fw/fh > 0: forced shape.
    Cand best_for(int b, int fwid, int fhei, int thr, const std::vector<int> &ms) const {
        Cand best;
        int bm = P.bmask[b];
        std::vector<std::pair<int, int>> shp;
        for (const Rc &F : M.f) {
            if (F.w < minside[b] || F.h < minside[b]) continue;
            if ((bm & lab3::EDGE_LEFT) && F.x != 0) continue;
            if ((bm & lab3::EDGE_RIGHT) && F.x1() != P.W) continue;
            if ((bm & lab3::EDGE_BOTTOM) && F.y != 0) continue;
            if ((bm & lab3::EDGE_TOP) && F.y1() != P.H) continue;
            shp.clear();
            if (fwid > 0) { if (fwid <= F.w && fhei <= F.h) shp.push_back({fwid, fhei}); }
            else soft_shapes(b, F, ms, shp);
            for (auto [w, h] : shp) {
                for (int corner = 0; corner < 4; ++corner) {
                    Rc c{(corner & 1) ? F.x1() - w : F.x, (corner & 2) ? F.y1() - h : F.y, w, h};
                    if ((bm & lab3::EDGE_LEFT) && c.x != 0) continue;
                    if ((bm & lab3::EDGE_RIGHT) && c.x1() != P.W) continue;
                    if ((bm & lab3::EDGE_BOTTOM) && c.y != 0) continue;
                    if ((bm & lab3::EDGE_TOP) && c.y1() != P.H) continue;
                    Sc s;
                    s.dead = dead_of(c, thr);
                    if (best.ok && s.dead > best.s.dead) continue;
                    if (o.score == 0) s.walls = walls_of(c);
                    s.dist = dist_of(b, c);
                    if (!best.ok || better(s, best.s)) { best.ok = true; best.r = c; best.s = s; }
                }
            }
        }
        return best;
    }

    // Dead threshold and snap offsets from the blocks still to place after `cur`.
    void thr_ms(const std::vector<int> &todo, size_t cur, int &thr, std::vector<int> &ms) const {
        thr = INT32_MAX;
        ms.clear();
        for (size_t k = cur + 1; k < todo.size(); ++k) {
            int b = todo[k];
            thr = std::min(thr, minside[b]);
            if (ms.size() < 3 && std::find(ms.begin(), ms.end(), minside[b]) == ms.end()) ms.push_back(minside[b]);
        }
        if (thr == INT32_MAX) thr = 1;
    }

    // Forced shape for b (fixed block, or soft group already sized).
    void forced(int b, int &w, int &h) const {
        w = h = -1;
        if (P.kind[b] != K_SOFT) { w = P.fw[b]; h = P.fh[b]; return; }
        int g = P.group[b];
        if (gsw[g] > 0) { w = gsw[g]; h = gsh[g]; }
    }

    // Place an MIB group (members in `mem`, todo index of its first member `cur`).
    // Tries up to 8 shapes ranked by the first member's best candidates.
    void place_mib(const std::vector<int> &mem, const std::vector<int> &todo, size_t cur_last, std::vector<int> &queue) {
        int g = P.group[mem[0]];
        int thr;
        std::vector<int> ms;
        thr_ms(todo, cur_last, thr, ms);
        if (gsw[g] > 0) {   // shape already fixed (repair re-pour): sequential
            for (int b : mem) {
                Cand c = best_for(b, gsw[g], gsh[g], thr, ms);
                if (c.ok) commit(b, c.r); else queue.push_back(b);
            }
            return;
        }
        // candidate shapes: best placement of the first member per distinct shape
        std::vector<std::pair<Sc, std::pair<int, int>>> shapes;
        {
            int b = mem[0];
            std::vector<std::pair<int, int>> shp;
            std::vector<std::pair<int, int>> seen;
            for (const Rc &F : M.f) {
                shp.clear();
                soft_shapes(b, F, ms, shp);
                for (auto wh : shp)
                    if (std::find(seen.begin(), seen.end(), wh) == seen.end()) seen.push_back(wh);
            }
            for (auto wh : seen) {
                Cand c = best_for(b, wh.first, wh.second, thr, ms);
                if (c.ok) shapes.push_back({c.s, wh});
            }
            std::sort(shapes.begin(), shapes.end(), [&](const auto &a, const auto &b2) { return better(a.first, b2.first); });
            if (shapes.size() > 8) shapes.resize(8);
        }
        int bestcnt = -1;
        int64_t bestdead = 0;
        std::pair<int, int> bestwh{-1, -1};
        std::vector<Rc> bestpos;
        MerSet M0 = M;
        for (auto &sh : shapes) {
            int cnt = 0;
            int64_t dead = 0;
            std::vector<Rc> rs;
            for (int b : mem) {
                Cand c = best_for(b, sh.second.first, sh.second.second, thr, ms);
                if (!c.ok) { rs.push_back(Rc{-1, -1, 0, 0}); continue; }
                ++cnt;
                dead += c.s.dead;
                rs.push_back(c.r);
                pos[b] = c.r; placed[b] = 1; M.place(c.r);
            }
            for (int b : mem) placed[b] = 0;
            M = M0;
            if (cnt > bestcnt || (cnt == bestcnt && dead < bestdead)) { bestcnt = cnt; bestdead = dead; bestwh = sh.second; bestpos = rs; }
            if (cnt == (int)mem.size() && dead == 0) break;
        }
        if (bestcnt <= 0) {
            for (int b : mem) queue.push_back(b);
            return;
        }
        gsw[g] = bestwh.first; gsh[g] = bestwh.second;
        for (size_t k = 0; k < mem.size(); ++k) {
            if (bestpos[k].w > 0) commit(mem[k], bestpos[k]);
            else queue.push_back(mem[k]);
        }
    }

    // Pour the given blocks in the given order (MIB members must be adjacent).
    void pour(const std::vector<int> &todo0, std::vector<int> &queue) {
        int thr;
        std::vector<int> ms;
        std::vector<int> todo = todo0;
        for (size_t k = 0; k < todo.size();) {
            if (o.bestfit > 0 && !is_mib(todo[k]) && P.bmask[todo[k]] == 0) {
                // best fit: among the next K plain blocks pick the (block, spot) with the least
                // dead area, then most walls, then the larger block, then distance.
                size_t jb = k;
                Cand cb;
                bool qd = false;
                for (size_t j = k; j < todo.size() && j < k + (size_t)o.bestfit; ++j) {
                    int b = todo[j];
                    if (is_mib(b) || P.bmask[b]) break;
                    std::swap(todo[k], todo[j]);   // thr_ms excludes position k
                    thr_ms(todo, k, thr, ms);
                    std::swap(todo[k], todo[j]);
                    int w, h;
                    forced(b, w, h);
                    Cand c = best_for(b, w, h, thr, ms);
                    if (!c.ok) { jb = j; qd = true; break; }
                    bool bt = !cb.ok || c.s.dead < cb.s.dead ||
                              (c.s.dead == cb.s.dead && (o.score == 0 && c.s.walls != cb.s.walls ? c.s.walls > cb.s.walls
                                                         : barea[b] != barea[todo[jb]] ? barea[b] > barea[todo[jb]] : c.s.dist < cb.s.dist));
                    if (bt) { cb = c; jb = j; }
                }
                std::swap(todo[k], todo[jb]);
                if (qd) queue.push_back(todo[k]);
                else commit(todo[k], cb.r);
                ++k;
                continue;
            }
            int b = todo[k];
            int g = P.kind[b] == K_SOFT ? P.group[b] : -1;
            if (g >= 0 && P.gmem[g].size() > 1) {
                std::vector<int> mem;
                size_t e = k;
                while (e < todo.size() && P.kind[todo[e]] == K_SOFT && P.group[todo[e]] == g) mem.push_back(todo[e++]);
                place_mib(mem, todo, e - 1, queue);
                k = e;
                continue;
            }
            thr_ms(todo, k, thr, ms);
            int w, h;
            forced(b, w, h);
            Cand c = best_for(b, w, h, thr, ms);
            if (c.ok) commit(b, c.r); else queue.push_back(b);
            ++k;
        }
    }

    // Order key: class then area desc. Classes: 0 corner blocks and MIB groups
    // with a corner member, 1 other MIB groups, 2 edge blocks, 3 fixed, 4 soft.
    bool is_mib(int b) const { return P.kind[b] == K_SOFT && P.gmem[P.group[b]].size() > 1; }
    int cls(int b) const {
        if (is_mib(b)) {
            for (int m : P.gmem[P.group[b]]) if (__builtin_popcount(P.bmask[m]) == 2) return 0;
            return 1;
        }
        int bm = P.bmask[b];
        if (bm && __builtin_popcount(bm) == 2) return 0;
        if (o.order == 1 && bm) return 0;   // all boundary blocks before MIB groups
        if (bm) return 2;
        if (P.kind[b] == K_FIXED) return 3;
        return 4;
    }
    std::vector<int> order(const std::vector<int> &ids) const {
        std::vector<int> v = ids;
        // MIB members stay adjacent (unit key = group), boundary members first.
        auto unit = [&](int b) { return is_mib(b) ? n + P.group[b] : b; };
        std::stable_sort(v.begin(), v.end(), [&](int a, int b) {
            int ca = cls(a), cb = cls(b);
            if (ca != cb) return ca < cb;
            int ua = unit(a), ub = unit(b);
            if (ua != ub) return barea[a] != barea[b] ? barea[a] > barea[b] : ua < ub;
            return (P.bmask[a] != 0) > (P.bmask[b] != 0);
        });
        return v;
    }

    int64_t qarea(const std::vector<int> &q) const {
        int64_t s = 0;
        for (int b : q) s += barea[b];
        return s;
    }

    // ---- repair
    void repair(std::vector<int> &queue, Rng &rng, Stats &st) {
        auto t0 = std::chrono::steady_clock::now();
        auto ms = [&]() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(); };
        int fails = 0;
        while (!queue.empty() && ms() < o.repair_ms && fails < o.max_fails) {
            ++st.rep_tries;
            // pick the window centre: a queued block's target (largest first, then random) or a big MER
            std::vector<int> qs = queue;
            std::sort(qs.begin(), qs.end(), [&](int a, int b) { return barea[a] > barea[b]; });
            int qb = qs[std::min<int>((int)qs.size() - 1, fails % 3 == 0 ? 0 : rng.below((int)qs.size()))];
            double cxw, cyw;
            int mode = rng.below(2);
            if (mode == 0 && tx) { cxw = (*tx)[qb]; cyw = (*ty)[qb]; }
            else {
                // a random one of the largest MERs
                std::vector<int> idx(M.f.size());
                for (size_t i = 0; i < idx.size(); ++i) idx[i] = (int)i;
                std::sort(idx.begin(), idx.end(), [&](int a, int b) { return M.f[a].area() > M.f[b].area(); });
                const Rc &F = M.f[idx[rng.below(std::min<int>(4, (int)idx.size()))]];
                cxw = F.x + F.w / 2.0; cyw = F.y + F.h / 2.0;
            }
            int k = 2 + rng.below(o.evict_max - 1);   // evict 2..evict_max nearest soft (non-MIB) neighbours
            std::vector<std::pair<double, int>> nb;
            for (int i = 0; i < n; ++i) {
                if (!placed[i] || P.kind[i] == K_PRE || is_mib(i)) continue;
                if (P.kind[i] == K_FIXED && !o.evict_fixed) continue;
                double d = std::fabs(pos[i].x + pos[i].w / 2.0 - cxw) + std::fabs(pos[i].y + pos[i].h / 2.0 - cyw);
                nb.push_back({d, i});
            }
            if (nb.empty()) break;
            std::sort(nb.begin(), nb.end());
            k = std::min<int>(k, (int)nb.size());
            // save
            std::vector<Rc> spos = pos;
            std::vector<char> spl = placed;
            std::vector<int> sgw = gsw, sgh = gsh;
            MerSet sM = M;
            int64_t before = qarea(queue);
            std::vector<int> ids = queue;
            for (int j = 0; j < k; ++j) {
                int b = nb[j].second;
                placed[b] = 0;
                if (P.kind[b] == K_SOFT) gsw[P.group[b]] = gsh[P.group[b]] = -1;
                ids.push_back(b);
            }
            reset_mers();
            std::vector<int> todo = order(ids);
            std::vector<int> nq;
            pour(todo, nq);
            int64_t after = qarea(nq);
            if (after < before || (o.plateau && after == before)) { queue = nq; ++st.rep_acc; if (after < before) fails = 0; else ++fails; }
            else { pos = spos; placed = spl; gsw = sgw; gsh = sgh; M = sM; ++fails; }
        }
        st.ms_repair = ms();
    }

    // ---- driver
    Stats run(Rng &rng) {
        Stats st;
        auto t0 = std::chrono::steady_clock::now();
        M.init(P.W, P.H);
        std::vector<int> ids;
        for (int i = 0; i < n; ++i) {
            st.total_area += barea[i];
            if (P.kind[i] == K_PRE) commit(i, Rc{P.px[i], P.py[i], P.fw[i], P.fh[i]});
            else ids.push_back(i);
        }
        std::vector<int> todo = order(ids), queue;
        pour(todo, queue);
        st.ms_greedy = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        if (o.repair && !queue.empty()) repair(queue, rng, st);
        st.unplaced = (int)queue.size();
        st.placed = n - st.unplaced;
        st.qarea = qarea(queue);
        for (const Rc &F : M.f) if (F.area() > st.maxmer.area()) st.maxmer = F;
        st.nmer = (int)M.f.size();
        st.max_nmer = max_nmer;
        return st;
    }

    std::vector<lab3::Rect> rects() const {
        std::vector<lab3::Rect> r(n);
        for (int i = 0; i < n; ++i) r[i] = lab3::Rect{pos[i].x, pos[i].y, pos[i].w, pos[i].h};
        return r;
    }
};

}  // namespace cx
}  // namespace fp

#endif  // FP_CONSTRUCT_H
