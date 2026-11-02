#!/usr/bin/env python3
"""E3 gating experiment: exact LP "WL polish" with the topology fixed.

Research prototype (not part of the release). For every dev case it:
  1. builds H/V constraint graphs from the construction geometry (diagonal
     pairs resolved by the larger-gap rule or by coin flips -> topology variants),
  2. computes the packed placement (ASAP longest path, preplaced pinned,
     RIGHT/TOP sinks flushed),
  3. solves the per-axis WL LP on the same graphs (scipy HiGHS dual simplex),
  4. cross-checks a few LPs with a min-cost-circulation dual (networkx
     network_simplex) and recovers the placement from its potentials,
  5. verifies every written placement with ./verifier/verifier.

LP formulation (integer lower-left coordinates, 2*WL units).  For a
connection of weight c between a and b, with d = x_a - x_b and
delta = w_a - w_b (terminal: x_b = 0 is the root and delta = w_a - 2*x_t):
    2 * c * |cx_a - cx_b| = c * |2d + delta|
  even delta: one hinge  2c * |d + delta/2|
  odd  delta: two hinges  c * |d + (delta-1)/2| + c * |d + (delta+1)/2|
(exact for integer d because |k| + |k+1| = |2k+1|).  All constraints are
differences x_j - x_i >= l, so the LP is the dual of a network flow and has
an integral optimal vertex.

Usage:
  python tools/lp_polish_proto.py [--cases 'cases/dev/*.in'] [--variants 8]
        [--mcf 1] [--out-dir DIR] [--json results.json]
"""
import argparse
import glob
import json
import os
import random
import re
import statistics
import subprocess
import sys
import time
import zlib

import networkx as nx
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VERIFIER = os.path.join(ROOT, "verifier", "verifier")


# ----------------------------------------------------------------- parsing

def parse_input(path):
    tok = open(path).read().split()
    p = 0

    def take(k=1):
        nonlocal p
        r = tok[p:p + k]
        p += k
        return r

    case = {"blocks": {}, "terms": {}, "conns": [], "mibs": [], "bnd": {}}
    _, W, H = take(3)
    case["W"], case["H"] = int(W), int(H)
    take(2)  # WLTHRESHOLD T
    take(2)  # ASPECTRATIO R
    _, n = take(2)
    for _ in range(int(n)):
        kind = take(1)[0]
        if kind == "SOFT":
            name, area = take(2)
            case["blocks"][name] = {"kind": "SOFT", "area": int(area)}
        elif kind == "FIXED":
            name, w, h = take(3)
            case["blocks"][name] = {"kind": "FIXED", "w": int(w), "h": int(h)}
        else:
            name, w, h, x, y = take(5)
            case["blocks"][name] = {"kind": "PREPLACED", "w": int(w), "h": int(h),
                                    "x": int(x), "y": int(y)}
    _, r = take(2)
    for _ in range(int(r)):
        _, name, x, y = take(4)
        case["terms"][name] = (int(x), int(y))
    _, m = take(2)
    for _ in range(int(m)):
        _, a, b, c = take(4)
        case["conns"].append((a, b, int(c)))
    _, q = take(2)
    for _ in range(int(q)):
        _, k = take(2)
        case["mibs"].append(take(int(k)))
    _, pb = take(2)
    for _ in range(int(pb)):
        take(1)
        name = take(1)[0]
        edges = []
        while p < len(tok) and tok[p] in ("LEFT", "RIGHT", "TOP", "BOTTOM"):
            edges.append(take(1)[0])
        case["bnd"][name] = edges
    return case


def parse_output(path):
    lines = open(path).read().split("\n")
    rects = {}
    for ln in lines[2:]:
        t = ln.split()
        if len(t) == 5:
            rects[t[0]] = tuple(int(v) for v in t[1:])
    return rects  # name -> (x, y, w, h)


def write_output(path, names, X, Y, Wd, Hd, wl):
    with open(path, "w") as f:
        f.write("WIRELENGTH %.1f\nBLOCKS %d\n" % (wl, len(names)))
        for i, nm in enumerate(names):
            f.write("%s %d %d %d %d\n" % (nm, X[i], Y[i], Wd[i], Hd[i]))


# ------------------------------------------------------------- wirelength

def wirelength2(case, idx, X, Y, Wd, Hd):
    """Exact 2*WL (integer) in doubled coordinates."""
    tot = 0
    for a, b, c in case["conns"]:
        pa = (2 * X[idx[a]] + Wd[idx[a]], 2 * Y[idx[a]] + Hd[idx[a]]) if a in idx \
            else (2 * case["terms"][a][0], 2 * case["terms"][a][1])
        pb = (2 * X[idx[b]] + Wd[idx[b]], 2 * Y[idx[b]] + Hd[idx[b]]) if b in idx \
            else (2 * case["terms"][b][0], 2 * case["terms"][b][1])
        tot += c * (abs(pa[0] - pb[0]) + abs(pa[1] - pb[1]))
    return tot  # doubled-coordinate distances, so tot = 2*WL


# --------------------------------------------------------- constraint graphs

def build_graphs(X, Y, Wd, Hd, rule, rng):
    """Return (Hedges, Vedges) as lists of (i, j): i left of / below j."""
    n = len(X)
    He, Ve = [], []
    for i in range(n):
        for j in range(i + 1, n):
            ox = min(X[i] + Wd[i], X[j] + Wd[j]) - max(X[i], X[j])
            oy = min(Y[i] + Hd[i], Y[j] + Hd[j]) - max(Y[i], Y[j])
            assert not (ox > 0 and oy > 0), "construction overlaps"
            if oy > 0:
                horiz = True
            elif ox > 0:
                horiz = False
            elif rule == "gap":
                horiz = -ox >= -oy  # larger separation axis
            else:
                horiz = rng.random() < 0.5
            if horiz:
                He.append((i, j) if X[i] + Wd[i] <= X[j] else (j, i))
            else:
                Ve.append((i, j) if Y[i] + Hd[i] <= Y[j] else (j, i))
    return He, Ve


def transitive_reduction(n, edges):
    g = nx.DiGraph()
    g.add_nodes_from(range(n))
    g.add_edges_from(edges)
    return list(nx.transitive_reduction(g).edges())


# ----------------------------------------------------------------- packing

def pack_axis(n, edges, size, L, pinned, low_pin, high_pin):
    """ASAP longest path. pinned: i -> coord. Returns coords or None."""
    g = nx.DiGraph()
    g.add_nodes_from(range(n))
    g.add_edges_from(edges)
    pos = [0] * n
    for v in nx.topological_sort(g):
        lo = max([pos[u] + size[u] for u in g.predecessors(v)], default=0)
        if v in pinned:
            if lo > pinned[v]:
                return None
            pos[v] = pinned[v]
        else:
            pos[v] = lo
    for v in low_pin:
        pos[v] = 0
    for v in high_pin:
        assert g.out_degree(v) == 0, "boundary-high block is not a sink"
        pos[v] = L - size[v]
    return pos


# ------------------------------------------------------------------ LP

def hinges_axis(case, idx, coord_size, term_axis):
    """List of (a, b, k, weight) for weight*|x_a - x_b + k|, b = -1 is root."""
    out = []
    for a, b, c in case["conns"]:
        if a not in idx:
            a, b = b, a
        ia = idx[a]
        if b in idx:
            ib = idx[b]
            delta = coord_size[ia] - coord_size[ib]
        else:
            ib = -1
            delta = coord_size[ia] - 2 * term_axis(case["terms"][b])
        if delta % 2 == 0:
            out.append((ia, ib, delta // 2, 2 * c))
        else:
            out.append((ia, ib, (delta - 1) // 2, c))
            out.append((ia, ib, (delta + 1) // 2, c))
    return out


def solve_lp_axis(n, edges, size, L, fixed, hinges):
    """min sum wt*t_h, t_h >= +-(x_a - x_b + k).  Returns (x, obj2, info)."""
    m = len(hinges)
    nv = n + m
    rows, cols, vals, rhs = [], [], [], []
    r = 0
    for h, (a, b, k, wt) in enumerate(hinges):
        for s in (1, -1):  # s*(x_a - x_b + k) - t_h <= 0
            rows += [r, r]
            cols += [a, n + h]
            vals += [s, -1]
            if b >= 0:
                rows.append(r)
                cols.append(b)
                vals.append(-s)
            rhs.append(-s * k)
            r += 1
    for i, j in edges:  # x_i - x_j <= -w_i
        rows += [r, r]
        cols += [i, j]
        vals += [1, -1]
        rhs.append(-size[i])
        r += 1
    A = coo_matrix((vals, (rows, cols)), shape=(r, nv)).tocsr()
    cost = np.zeros(nv)
    cost[n:] = [h[3] for h in hinges]
    bounds = []
    for i in range(n):
        if i in fixed:
            bounds.append((fixed[i], fixed[i]))
        else:
            bounds.append((0, L - size[i]))
    bounds += [(0, None)] * m
    t0 = time.perf_counter()
    res = linprog(cost, A_ub=A, b_ub=np.array(rhs, float), bounds=bounds,
                  method="highs-ds")
    dt = time.perf_counter() - t0
    assert res.status == 0, res.message
    xs = res.x[:n]
    frac = int(np.sum(np.abs(xs - np.round(xs)) > 1e-6))
    x = [int(round(v)) for v in xs]
    return x, res.fun, {"time": dt, "vars": nv, "rows": r, "frac": frac,
                        "hinges": m}


def solve_mcf_axis(n, edges, size, L, fixed, hinges):
    """Dual min-cost circulation; returns (x from potentials, LP optimum)."""
    R = n  # root node id
    G = nx.MultiDiGraph()
    G.add_nodes_from(range(n + 1))

    def node(v):
        return R if v < 0 else v

    for a, b, k, wt in hinges:  # profit k on a->b, -k on b->a, cap wt
        G.add_edge(node(a), node(b), weight=-k, capacity=wt)
        G.add_edge(node(b), node(a), weight=k, capacity=wt)
    for i, j in edges:  # x_j - x_i >= w_i: arc i->j, profit w_i, cap inf
        G.add_edge(i, j, weight=-size[i])
    for i in range(n):
        lo, hi = (fixed[i], fixed[i]) if i in fixed else (0, L - size[i])
        G.add_edge(R, i, weight=-lo)  # x_i - x_R >= lo
        G.add_edge(i, R, weight=hi)   # x_R - x_i >= -hi
    t0 = time.perf_counter()
    cost, flow = nx.network_simplex(G)
    dt = time.perf_counter() - t0
    # residual graph -> shortest-path potentials -> x = -pi (x_R = 0)
    Rg = nx.DiGraph()
    Rg.add_nodes_from(range(n + 1))

    def add(u, v, w):
        if not Rg.has_edge(u, v) or Rg[u][v]["weight"] > w:
            Rg.add_edge(u, v, weight=w)

    for u, v, key, d in G.edges(keys=True, data=True):
        f = flow[u][v][key]
        if "capacity" not in d or f < d["capacity"]:
            add(u, v, d["weight"])
        if f > 0:
            add(v, u, -d["weight"])
    pi = nx.single_source_bellman_ford_path_length(Rg, R)
    x = [-(pi[i] - pi[R]) for i in range(n)]
    return x, -cost, dt


def lp_objective(x, hinges):
    return sum(wt * abs(x[a] - (x[b] if b >= 0 else 0) + k)
               for a, b, k, wt in hinges)


# ---------------------------------------------------------------- verifier

def verify(inp, out):
    p = subprocess.run([VERIFIER, inp, out], capture_output=True, text=True)
    txt = p.stdout + p.stderr
    m = re.search(r"Computed WL\s+([\d.]+)", txt)
    bnd = re.search(r"Boundary\s+[\d.]+ / 3\s+\((\d+)/(\d+)", txt)
    mib = re.search(r"MIB\s+[\d.]+ / 2\s+\((\d+)/(\d+)", txt)
    return {
        "feasible": p.returncode == 0 and "Result: FEASIBLE" in txt,
        "wl": float(m.group(1)) if m else None,
        "bnd_full": bool(bnd) and bnd.group(1) == bnd.group(2),
        "mib_full": bool(mib) and mib.group(1) == mib.group(2),
    }


# ------------------------------------------------------------------ driver

def run_case(inp, args):
    case = parse_input(inp)
    cons = parse_output(inp[:-3] + ".construction.out")
    names = sorted(case["blocks"])
    idx = {nm: i for i, nm in enumerate(names)}
    n = len(names)
    X = [cons[nm][0] for nm in names]
    Y = [cons[nm][1] for nm in names]
    Wd = [cons[nm][2] for nm in names]
    Hd = [cons[nm][3] for nm in names]
    W, H = case["W"], case["H"]
    pinx = {idx[nm]: b["x"] for nm, b in case["blocks"].items() if b["kind"] == "PREPLACED"}
    piny = {idx[nm]: b["y"] for nm, b in case["blocks"].items() if b["kind"] == "PREPLACED"}
    side = {e: [idx[nm] for nm, es in case["bnd"].items() if e in es]
            for e in ("LEFT", "RIGHT", "BOTTOM", "TOP")}
    fx = dict(pinx)
    fx.update({i: 0 for i in side["LEFT"]})
    fx.update({i: W - Wd[i] for i in side["RIGHT"]})
    fy = dict(piny)
    fy.update({i: 0 for i in side["BOTTOM"]})
    fy.update({i: H - Hd[i] for i in side["TOP"]})
    hx = hinges_axis(case, idx, Wd, lambda t: t[0])
    hy = hinges_axis(case, idx, Hd, lambda t: t[1])

    base = os.path.splitext(os.path.basename(inp))[0]
    cons_wl2 = wirelength2(case, idx, X, Y, Wd, Hd)
    rec = {"case": base, "n": n, "construction_wl": cons_wl2 / 2, "variants": [],
           "infeasible_packings": 0}
    for v in range(args.variants):
        rule = "gap" if v == 0 else "rand"
        rng = random.Random(zlib.crc32(("%s/%d" % (base, v)).encode()))
        He, Ve = build_graphs(X, Y, Wd, Hd, rule, rng)
        HeR, VeR = transitive_reduction(n, He), transitive_reduction(n, Ve)
        px = pack_axis(n, HeR, Wd, W, pinx, side["LEFT"], side["RIGHT"])
        py = pack_axis(n, VeR, Hd, H, piny, side["BOTTOM"], side["TOP"])
        if px is None or py is None or any(px[i] + Wd[i] > W or py[i] + Hd[i] > H
                                           for i in range(n)):
            rec["infeasible_packings"] += 1
            continue
        pk_wl = wirelength2(case, idx, px, py, Wd, Hd) / 2
        f_pk = os.path.join(args.out_dir, "%s_v%d.packed.out" % (base, v))
        write_output(f_pk, names, px, py, Wd, Hd, pk_wl)
        vp = verify(inp, f_pk)
        if not (vp["feasible"] and vp["bnd_full"] and vp["mib_full"]):
            rec["infeasible_packings"] += 1
            continue
        assert abs(vp["wl"] - pk_wl) < 0.01

        lx, ox, ix = solve_lp_axis(n, HeR, Wd, W, fx, hx)
        ly, oy, iy = solve_lp_axis(n, VeR, Hd, H, fy, hy)
        assert lp_objective(lx, hx) == round(ox) and lp_objective(ly, hy) == round(oy)
        lp_wl = wirelength2(case, idx, lx, ly, Wd, Hd) / 2
        assert abs(lp_wl - (ox + oy) / 2) < 0.01, (lp_wl, ox, oy)
        f_lp = os.path.join(args.out_dir, "%s_v%d.lp.out" % (base, v))
        write_output(f_lp, names, lx, ly, Wd, Hd, lp_wl)
        vl = verify(inp, f_lp)
        vr = {"rule": rule, "wl_packed": pk_wl, "wl_lp": lp_wl,
              "gain": (pk_wl - lp_wl) / pk_wl,
              "lp_feasible": vl["feasible"], "lp_bnd_full": vl["bnd_full"],
              "lp_mib_full": vl["mib_full"],
              "lp_wl_match": vl["wl"] is not None and abs(vl["wl"] - lp_wl) < 0.01,
              "edges_h": len(He), "edges_v": len(Ve),
              "edges_h_tr": len(HeR), "edges_v_tr": len(VeR),
              "lp_x": ix, "lp_y": iy}
        if v == 0:
            # Extra: construction as-is -> LP.  The LP depends only on the
            # graphs, so this is the deterministic variant's LP optimum.
            vr["construction_gain"] = (cons_wl2 / 2 - lp_wl) / (cons_wl2 / 2)
        if v < args.mcf:
            mx, mox, tx = solve_mcf_axis(n, HeR, Wd, W, fx, hx)
            my, moy, ty = solve_mcf_axis(n, VeR, Hd, H, fy, hy)
            mwl = wirelength2(case, idx, mx, my, Wd, Hd) / 2
            f_mc = os.path.join(args.out_dir, "%s_v%d.mcf.out" % (base, v))
            write_output(f_mc, names, mx, my, Wd, Hd, mwl)
            vm = verify(inp, f_mc)
            vr["mcf"] = {"opt_x": mox, "opt_y": moy, "lp_x": ox, "lp_y": oy,
                         "wl_from_potentials": mwl, "time_x": tx, "time_y": ty,
                         "match": round(mox) == round(ox) and round(moy) == round(oy)
                         and abs(mwl - lp_wl) < 0.01,
                         "verifier_feasible": vm["feasible"] and vm["bnd_full"]
                         and vm["mib_full"]}
        rec["variants"].append(vr)
    wp = [v["wl_packed"] for v in rec["variants"]]
    wlp = [v["wl_lp"] for v in rec["variants"]]
    if len(wp) >= 3 and len(set(wp)) > 1 and len(set(wlp)) > 1:
        rec["spearman"] = float(spearmanr(wp, wlp).statistic)
    else:
        rec["spearman"] = None
    return rec


def pct(xs, q):
    return float(np.percentile(xs, q))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cases", default=os.path.join(ROOT, "cases/dev/*.in"))
    ap.add_argument("--variants", type=int, default=8)
    ap.add_argument("--mcf", type=int, default=1,
                    help="per case, cross-check the first K variants with MCF")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    cases = sorted(glob.glob(args.cases))
    recs = []
    for inp in cases:
        rec = run_case(inp, args)
        recs.append(rec)
        g = [v["gain"] * 100 for v in rec["variants"]]
        print("%-16s n=%3d var=%d inf=%d gain med=%.2f%% [%.2f..%.2f] rho=%s cons->LP=%.2f%%"
              % (rec["case"], rec["n"], len(g), rec["infeasible_packings"],
                 statistics.median(g), min(g), max(g),
                 "%.2f" % rec["spearman"] if rec["spearman"] is not None else "n/a",
                 rec["variants"][0].get("construction_gain", float("nan")) * 100),
              flush=True)
    allv = [v for r in recs for v in r["variants"]]
    gains = [v["gain"] * 100 for v in allv]
    print("packings=%d  gain median=%.2f%% p10=%.2f%% p90=%.2f%%"
          % (len(gains), statistics.median(gains), pct(gains, 10), pct(gains, 90)))
    bad = [v for v in allv if not (v["lp_feasible"] and v["lp_bnd_full"]
                                   and v["lp_mib_full"] and v["lp_wl_match"])]
    print("LP outputs failing verifier/boundary/MIB/WL-match: %d" % len(bad))
    print("LP vertices with fractional x: %d" % sum(v["lp_x"]["frac"] + v["lp_y"]["frac"]
                                                    for v in allv))
    mc = [v["mcf"] for v in allv if "mcf" in v]
    print("MCF cross-checks: %d, all match=%s, all verifier-feasible=%s"
          % (len(mc), all(m["match"] for m in mc), all(m["verifier_feasible"] for m in mc)))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(recs, f, indent=1)


if __name__ == "__main__":
    sys.exit(main())
