#!/usr/bin/env python3
"""Lab3 fixed-outline floorplanning test-case generator (docs/design.md section 8).

Builds a known zero-whitespace dissection of the outline (slicing cuts mixed
with pinwheel nodes, MIB groups as equal strips), shrinks every cell to a
block to inject whitespace, assigns block kinds and constraints, and emits
`<out>.in` plus `<out>.construction.out` (a feasible solution with full MIB
and boundary points whose WIRELENGTH is the exact construction WL).

Python 3.10, stdlib + numpy only. Deterministic for a given seed.
"""

import argparse
import math
import os
import sys

# single-threaded numpy (avoid BLAS thread spin on shared hosts)
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402

ASPECT_R = 3            # ASPECTRATIO written to the input
CELL_ASPECT = 2.5       # max aspect ratio of a dissection cell
BLOCK_ASPECT = 2.8      # max aspect ratio of a shrunk block (< ASPECT_R)
AVG_CELL_SIDE = 480     # target mean cell side before shrinking
MIN_SOFT_AREA = 10000
COORD_MAX = 100000      # docs/design.md section 3 bounds
CONN_MAX = 1000

PRESETS = {
    #           soft fixed pre  term mib bnd corners ws
    "public":  (30, 14, 6, 24, 4, 10, 2, 0.15),
    "hidden1": (48, 22, 10, 36, 6, 16, 3, 0.12),
    "hidden2": (72, 34, 14, 48, 8, 24, 4, 0.10),
    "hidden3": (96, 45, 19, 60, 10, 32, 4, 0.09),
    "hidden4": (120, 56, 24, 72, 12, 40, 4, 0.08),
}

EDGE_NAMES = ("LEFT", "RIGHT", "BOTTOM", "TOP")


class GenError(Exception):
    """Recoverable failure; the generator retries with the next attempt."""


# ---------------------------------------------------------------------------
# Dissection
# ---------------------------------------------------------------------------

def aspect(w, h):
    return max(w, h) / min(w, h)


class Dissector:
    """Recursive zero-whitespace dissection of a rectangle into cells.

    A unit (mib_id or -1, count, weight, uid) is either a normal cell (count 1)
    or an MIB unit (count k) that is split into k equal strips. Cuts divide the
    area in proportion to the weights (= counts unless --size-tail/--big-cells
    is used); aspect plausibility uses the counts. Output cells:
    (x0, y0, x1, y1, mib_id or -1, uid).
    """

    def __init__(self, rng, budget=20000, p_pinwheel=0.5):
        self.rng = rng
        self.budget = budget
        self.p_pinwheel = p_pinwheel
        self.pinwheels = 0

    def plausible(self, rect, units):
        x0, y0, x1, y1 = rect
        w, h = x1 - x0, y1 - y0
        if w < 2 or h < 2:
            return False
        if len(units) == 1:
            return self.leaf_ok(rect, units[0]) is not None
        total = sum(u[1] for u in units)
        return aspect(w, h) <= CELL_ASPECT * total

    def leaf_ok(self, rect, unit):
        """Return orientation for a single unit ('n'/'v'/'h') or None."""
        x0, y0, x1, y1 = rect
        w, h = x1 - x0, y1 - y0
        mid, k = unit[0], unit[1]
        if mid < 0:
            return "n" if aspect(w, h) <= CELL_ASPECT else None
        opts = []
        if w // k >= 2 and aspect(w // k, h) <= CELL_ASPECT:
            opts.append("v")   # vertical strips side by side along x
        if h // k >= 2 and aspect(w, h // k) <= CELL_ASPECT:
            opts.append("h")
        if not opts:
            return None
        return opts[int(self.rng.integers(len(opts)))]

    def leaf(self, rect, unit):
        o = self.leaf_ok(rect, unit)
        if o is None:
            return None
        x0, y0, x1, y1 = rect
        mid, k, _, uid = unit
        if o == "n":
            return [(x0, y0, x1, y1, -1, uid)]
        cells = []
        if o == "v":
            s = (x1 - x0) // k
            for i in range(k):
                cx1 = x0 + s * (i + 1) if i < k - 1 else x1
                cells.append((x0 + s * i, y0, cx1, y1, mid, uid))
        else:
            s = (y1 - y0) // k
            for i in range(k):
                cy1 = y0 + s * (i + 1) if i < k - 1 else y1
                cells.append((x0, y0 + s * i, x1, cy1, mid, uid))
        return cells

    def run(self, rect, units):
        self.budget -= 1
        if self.budget < 0:
            raise GenError("dissection budget exhausted")
        if len(units) == 1:
            return self.leaf(rect, units[0])
        rng = self.rng
        for _ in range(10):
            order = [units[i] for i in rng.permutation(len(units))]
            if len(units) >= 5 and rng.random() < self.p_pinwheel:
                res = self.try_pinwheel(rect, order)
            else:
                res = self.try_slice(rect, order)
            if res is not None:
                return res
        return None

    def try_slice(self, rect, units):
        rng = self.rng
        x0, y0, x1, y1 = rect
        w, h = x1 - x0, y1 - y0
        vertical_cut = w >= h  # cut perpendicular to the longer side
        if rng.random() < 0.15:
            vertical_cut = not vertical_cut
        # area shares follow the unit weights (u[2]); they equal the strip
        # counts unless --size-tail / --big-cells is used
        total = sum(u[2] for u in units)
        frac = rng.uniform(0.3, 0.7)
        cum, best, besti = 0, None, 1
        for i in range(1, len(units)):
            cum += units[i - 1][2]
            d = abs(cum / total - frac)
            if best is None or d < best:
                best, besti = d, i
        a, b = units[:besti], units[besti:]
        ca = sum(u[2] for u in a)
        f = ca / total * rng.uniform(0.88, 1.14)
        L = w if vertical_cut else h
        cut = int(round(L * f))
        cut = min(max(cut, 2), L - 2)
        if vertical_cut:
            ra, rb = (x0, y0, x0 + cut, y1), (x0 + cut, y0, x1, y1)
        else:
            ra, rb = (x0, y0, x1, y0 + cut), (x0, y0 + cut, x1, y1)
        if not (self.plausible(ra, a) and self.plausible(rb, b)):
            return None
        la = self.run(ra, a)
        if la is None:
            return None
        lb = self.run(rb, b)
        if lb is None:
            return None
        return la + lb

    def try_pinwheel(self, rect, units):
        """Five-region pinwheel (non-slicing) node."""
        rng = self.rng
        x0, y0, x1, y1 = rect
        w, h = x1 - x0, y1 - y0
        a = x0 + int(round(w * rng.uniform(0.22, 0.45)))
        b = x0 + int(round(w * rng.uniform(0.55, 0.78)))
        c = y0 + int(round(h * rng.uniform(0.22, 0.45)))
        d = y0 + int(round(h * rng.uniform(0.55, 0.78)))
        if not (x0 < a < b < x1 and y0 < c < d < y1):
            return None
        regions = [
            (x0, y0, a, d),   # left, tall
            (a, y0, x1, c),   # bottom, wide
            (b, c, x1, y1),   # right, tall
            (x0, d, b, y1),   # top, wide
            (a, c, b, d),     # centre
        ]
        if rng.random() < 0.5:  # opposite chirality: mirror in x
            regions = [(x0 + x1 - rx1, ry0, x0 + x1 - rx0, ry1)
                       for (rx0, ry0, rx1, ry1) in regions]
        total = sum(u[2] for u in units)
        areas = [(r[2] - r[0]) * (r[3] - r[1]) for r in regions]
        tot_area = float(sum(areas))
        targets = [ar / tot_area * total for ar in areas]
        groups = [[] for _ in range(5)]
        counts = [0] * 5
        # every region gets one unit first (largest target first), then the
        # rest goes to the region with the largest remaining deficit
        order5 = sorted(range(5), key=lambda g: -targets[g])
        idx = 0
        for g in order5:
            groups[g].append(units[idx])
            counts[g] += units[idx][2]
            idx += 1
        for u in units[idx:]:
            g = max(range(5), key=lambda g: targets[g] - counts[g])
            groups[g].append(u)
            counts[g] += u[2]
        for g in range(5):
            if not self.plausible(regions[g], groups[g]):
                return None
        out = []
        for g in range(5):
            r = self.run(regions[g], groups[g])
            if r is None:
                return None
            out += r
        self.pinwheels += 1
        return out


# ---------------------------------------------------------------------------
# Generator
# ---------------------------------------------------------------------------

def touches(cell, W, H):
    x0, y0, x1, y1 = cell[:4]
    t = []
    if x0 == 0:
        t.append("LEFT")
    if x1 == W:
        t.append("RIGHT")
    if y0 == 0:
        t.append("BOTTOM")
    if y1 == H:
        t.append("TOP")
    return t


def wl2(blocks, terms, conns):
    """Exact 2*WL as an int. blocks: idx -> (x,y,w,h); terms: idx -> (x,y)."""
    def p2(e):
        kind, i = e
        if kind == "b":
            x, y, w, h = blocks[i]
            return 2 * x + w, 2 * y + h
        x, y = terms[i]
        return 2 * x, 2 * y
    tot = 0
    for a, b, c in conns:
        ax, ay = p2(a)
        bx, by = p2(b)
        tot += c * (abs(ax - bx) + abs(ay - by))
    return tot


def fmt_half(v2):
    """Format v2/2 with one decimal, exactly."""
    return f"{v2 // 2}.{5 if v2 % 2 else 0}"


def arranged_preplaced(args, rng, cells, interior, taken, want, W, H):
    """Interior non-MIB cells for --pre-arrange line/cross/ring/largest.

    Returns up to `want` cell indices not in `taken`; the caller fills the
    remainder with the default spread rule. Draws from rng only when a
    non-default arrangement is requested.
    """
    if args.pre_arrange == "spread" or want <= 0:
        return []
    avail = [i for i in interior if i not in taken]
    if args.pre_arrange == "largest":
        avail.sort(key=lambda i: (-(cells[i][2] - cells[i][0]) * (cells[i][3] - cells[i][1]), i))
        return avail[:want]

    def crossing(horizontal, pos):
        # cells whose span contains the line, ordered along the line
        if horizontal:
            cs = [i for i in avail if cells[i][1] <= pos < cells[i][3]]
            cs.sort(key=lambda i: cells[i][0])
        else:
            cs = [i for i in avail if cells[i][0] <= pos < cells[i][2]]
            cs.sort(key=lambda i: cells[i][1])
        return cs

    def contiguous(cs, k):
        if len(cs) <= k:
            return cs
        s = (len(cs) - k) // 2
        return cs[s:s + k]

    if args.pre_arrange in ("line", "cross"):
        horizontal = bool(rng.integers(2))
        u = rng.uniform(0.4, 0.6)
        out = contiguous(crossing(horizontal, int(round((H if horizontal else W) * u))), want)
        if args.pre_arrange == "cross" and len(out) < want:
            v = rng.uniform(0.4, 0.6)
            cs = [i for i in crossing(not horizontal, int(round((W if horizontal else H) * v)))
                  if i not in out]
            out += contiguous(cs, want - len(out))
        return out
    # ring: cells intersecting the boundary of the central rectangle
    # [0.25W, 0.75W] x [0.25H, 0.75H], ordered by angle around the centre
    rx0, rx1, ry0, ry1 = 0.25 * W, 0.75 * W, 0.25 * H, 0.75 * H
    cs = []
    for i in avail:
        x0, y0, x1, y1 = cells[i][:4]
        if x1 <= rx0 or x0 >= rx1 or y1 <= ry0 or y0 >= ry1:
            continue            # outside the ring rectangle
        if x0 >= rx0 and x1 <= rx1 and y0 >= ry0 and y1 <= ry1:
            continue            # strictly inside: not on the fence
        cs.append(i)
    cs.sort(key=lambda i: math.atan2((cells[i][1] + cells[i][3]) / 2 - H / 2,
                                     (cells[i][0] + cells[i][2]) / 2 - W / 2))
    if len(cs) > want:
        cs = [cs[int(j * len(cs) / want)] for j in range(want)]
    return cs


def generate(args, attempt):
    rng = np.random.default_rng([args.seed, attempt])
    S, F, P = args.soft, args.fixed, args.preplaced
    n = S + F + P
    Q, B, C, R = args.mib_groups, args.boundary, args.corners, args.terminals
    ws = args.whitespace

    # ---- MIB group sizes ------------------------------------------------
    if args.mib_k_min is None and args.mib_k_max is None:
        mib_k = [int(rng.integers(2, 5)) for _ in range(Q)]
    else:
        lo = args.mib_k_min if args.mib_k_min is not None else 2
        hi = args.mib_k_max if args.mib_k_max is not None else max(lo, 4)
        mib_k = [int(rng.integers(lo, hi + 1)) for _ in range(Q)]
    while sum(mib_k) > max(S - 2, 2 * Q):
        i = int(np.argmax(mib_k))
        if mib_k[i] <= 2:
            break
        mib_k[i] -= 1
    if sum(mib_k) > S:
        raise SystemExit("error: not enough soft blocks for MIB groups")

    # ---- outline and dissection ----------------------------------------
    ratio = rng.uniform(0.75, 1.33)
    if args.outline_aspect is not None:
        ratio = args.outline_aspect   # W/H
    area = n * AVG_CELL_SIDE ** 2
    W = int(round(math.sqrt(area * ratio)))
    H = int(round(area / W))
    nplain = n - sum(mib_k)
    # unit weights: integer counts by default (byte-identical output); with
    # --size-tail the per-strip weight is lognormal, with --big-cells K plain
    # units together take the area share --big-frac
    big_uids = set()
    if args.size_tail is None and args.big_cells is None:
        units = [(-1, 1, 1, u) for u in range(nplain)] + \
                [(g, k, k, nplain + g) for g, k in enumerate(mib_k)]
    else:
        wt = np.ones(nplain + Q)
        if args.size_tail is not None:
            wt = np.exp(rng.normal(0.0, args.size_tail, nplain + Q))
        for g, k in enumerate(mib_k):
            wt[nplain + g] *= k
        if args.big_cells:
            K = args.big_cells
            if K > nplain:
                raise SystemExit("error: --big-cells exceeds the non-MIB block count")
            big_uids = set(int(u) for u in rng.permutation(nplain)[:K])
            rest = sum(float(wt[u]) for u in range(nplain + Q) if u not in big_uids)
            frac = args.big_frac
            wbig = frac * rest / ((1.0 - frac) * K)
            for u in big_uids:
                wt[u] = wbig
        units = [(-1, 1, float(wt[u]), u) for u in range(nplain)] + \
                [(g, k, float(wt[nplain + g]), nplain + g) for g, k in enumerate(mib_k)]
    dis = Dissector(rng)
    cells = dis.run((0, 0, W, H), units)
    if cells is None or len(cells) != n:
        raise GenError("dissection failed")
    if dis.pinwheels < max(1, n // 40):
        raise GenError("too few non-slicing (pinwheel) nodes")

    cw = [c[2] - c[0] for c in cells]
    ch = [c[3] - c[1] for c in cells]
    cmib = [c[4] for c in cells]
    ctouch = [touches(c, W, H) for c in cells]
    assert sum(cw[i] * ch[i] for i in range(n)) == W * H

    # ---- kinds ------------------------------------------------------------
    kind = [None] * n
    for i in range(n):
        if cmib[i] >= 0:
            kind[i] = "S"
    interior = [i for i in range(n) if not ctouch[i] and cmib[i] < 0]
    if len(interior) < P:
        raise GenError("not enough interior non-MIB cells for preplaced")
    cx = [(c[0] + c[2]) / 2 for c in cells]
    cy = [(c[1] + c[3]) / 2 for c in cells]
    big_cells = [i for i in range(n) if cells[i][5] in big_uids]
    pre = []
    if args.big_kind == "pre":
        pre = [i for i in big_cells if i in interior][:P]
    pre += arranged_preplaced(args, rng, cells, interior, pre, P - len(pre), W, H)
    while len(pre) < P:
        avail = [i for i in interior if i not in pre]
        sample = [avail[j] for j in rng.permutation(len(avail))[:8]]
        if not pre:
            pick = sample[0]
        else:
            def mind(i):
                return min(abs(cx[i] - cx[j]) + abs(cy[i] - cy[j]) for j in pre)
            pick = max(sample, key=mind)
        pre.append(pick)
    for i in pre:
        kind[i] = "P"

    # boundary constraints
    corner_cells = []
    for want in (("LEFT", "BOTTOM"), ("RIGHT", "BOTTOM"), ("LEFT", "TOP"),
                 ("RIGHT", "TOP")):
        cs = [i for i in range(n) if want[0] in ctouch[i] and want[1] in ctouch[i]
              and kind[i] != "P"]
        assert len(cs) <= 1
        if cs:
            corner_cells.append((cs[0], list(want)))
    if len(corner_cells) < C:
        raise GenError("not enough corner cells")
    bnd = {}
    for j in rng.permutation(len(corner_cells))[:C]:
        i, e = corner_cells[j]
        bnd[i] = e
    edge_cands = [i for i in range(n) if ctouch[i] and kind[i] != "P"
                  and i not in bnd]
    nonmib = [i for i in edge_cands if cmib[i] < 0]
    mibc = [i for i in edge_cands if cmib[i] >= 0]
    order = [nonmib[j] for j in rng.permutation(len(nonmib))] + \
            [mibc[j] for j in rng.permutation(len(mibc))]
    if len(order) < B - C:
        raise GenError("not enough border cells for boundary constraints")
    for i in order[:B - C]:
        t = ctouch[i]
        bnd[i] = [t[int(rng.integers(len(t)))]]

    # soft / fixed among the rest (big cells first, when their kind is forced)
    if args.big_kind in ("pre", "fixed"):
        for i in big_cells:
            if kind[i] is None:
                kind[i] = "F"
        if sum(1 for k in kind if k == "F") > F:
            raise GenError("more forced-fixed big cells than --fixed")
    need_soft = S - sum(1 for k in kind if k == "S")
    rest = [i for i in range(n) if kind[i] is None]
    if need_soft < 0 or need_soft > len(rest):
        raise GenError("kind counts")
    rest = [rest[j] for j in rng.permutation(len(rest))]
    for i in rest[:need_soft]:
        kind[i] = "S"
    for i in rest[need_soft:]:
        kind[i] = "F"

    # ---- whitespace injection -----------------------------------------
    target = int(round((1.0 - ws) * W * H))
    # shrink "units": MIB group (shrunk identically, base = min member cell)
    groups = {}
    for i in range(n):
        if cmib[i] >= 0:
            groups.setdefault(cmib[i], []).append(i)
    shrink_units = []  # (members, base_w, base_h)
    for i in range(n):
        if cmib[i] < 0:
            shrink_units.append(([i], cw[i], ch[i]))
    for g in sorted(groups):
        mem = groups[g]
        shrink_units.append((mem, min(cw[i] for i in mem), min(ch[i] for i in mem)))
    U = len(shrink_units)
    base_area = np.array([len(m) * bw * bh for m, bw, bh in shrink_units], float)
    dmax = min(0.55, 3.0 * ws + 0.05)
    d = ws * np.exp(rng.normal(0.0, 0.6, U))
    filled = np.zeros(U, bool)
    if args.pre_fill:
        # preplaced blocks keep their whole cell: zero slack around obstacles
        filled = np.array([kind[m[0]] == "P" for m, _, _ in shrink_units])
        d[filled] = 0.0
    # include area lost to MIB remainder pixels
    lost = W * H - base_area.sum()
    need_cov = target
    scale = 1.0
    for _ in range(200):
        dd = np.clip(d * scale, 0.0, dmax)
        cov = (base_area * (1 - dd)).sum()
        if abs(cov - need_cov) < 1e-6 * need_cov:
            break
        # adjust scale
        miss = (cov - need_cov) / max((base_area * (dd < dmax)).sum(), 1.0)
        scale *= 1 + miss / max(dd.mean(), 1e-9)
        scale = max(scale, 0.0)
    dd = np.clip(d * scale, 0.0, dmax)
    del lost
    bw_ = [0] * n
    bh_ = [0] * n
    unit_dims = []
    for u, (mem, w0, h0) in enumerate(shrink_units):
        r = 1.0 - dd[u]
        best = None
        for _ in range(30):
            t = rng.uniform(0.0, 1.0)
            sx = r ** t
            sy = r ** (1.0 - t)
            w = max(1, min(w0, int(round(w0 * sx))))
            h = max(1, min(h0, int(round(h0 * sy))))
            if aspect(w, h) <= BLOCK_ASPECT:
                best = (w, h)
                break
        if best is None:
            sx = sy = math.sqrt(r)
            best = (max(1, min(w0, int(round(w0 * sx)))),
                    max(1, min(h0, int(round(h0 * sy)))))
            if aspect(*best) > BLOCK_ASPECT:
                raise GenError("aspect")
        unit_dims.append(list(best))

    def total_area():
        return sum(len(shrink_units[u][0]) * unit_dims[u][0] * unit_dims[u][1]
                   for u in range(U))

    # +-1 tweaks to hit the target closely
    tol = max(1, int(target * 2e-5))
    for _ in range(200000):
        diff = total_area() - target
        if abs(diff) <= tol:
            break
        u = int(rng.integers(U))
        if filled[u]:
            continue
        mem, w0, h0 = shrink_units[u]
        w, h = unit_dims[u]
        k = len(mem)
        if diff < 0:
            # grow by one in a dimension that fits
            opts = []
            if w + 1 <= w0 and aspect(w + 1, h) <= BLOCK_ASPECT and k * h <= -diff + tol:
                opts.append((w + 1, h))
            if h + 1 <= h0 and aspect(w, h + 1) <= BLOCK_ASPECT and k * w <= -diff + tol:
                opts.append((w, h + 1))
        else:
            opts = []
            if w > 1 and aspect(w - 1, h) <= BLOCK_ASPECT and k * h <= diff + tol:
                opts.append((w - 1, h))
            if h > 1 and aspect(w, h - 1) <= BLOCK_ASPECT and k * w <= diff + tol:
                opts.append((w, h - 1))
        if opts:
            unit_dims[u] = list(opts[int(rng.integers(len(opts)))])
    for u, (mem, _, _) in enumerate(shrink_units):
        for i in mem:
            bw_[i], bh_[i] = unit_dims[u]

    # ---- anchoring ----------------------------------------------------------
    bx = [0] * n
    by = [0] * n
    for i in range(n):
        x0, y0, x1, y1 = cells[i][:4]
        e = bnd.get(i, [])
        sx = x1 - x0 - bw_[i]
        sy = y1 - y0 - bh_[i]
        if "LEFT" in e:
            bx[i] = x0
        elif "RIGHT" in e:
            bx[i] = x1 - bw_[i]
        else:
            bx[i] = x0 + int(rng.integers(sx + 1))
        if "BOTTOM" in e:
            by[i] = y0
        elif "TOP" in e:
            by[i] = y1 - bh_[i]
        else:
            by[i] = y0 + int(rng.integers(sy + 1))

    # min soft area: integer upscaling preserves every property
    min_soft = min(bw_[i] * bh_[i] for i in range(n) if kind[i] == "S")
    if min_soft < MIN_SOFT_AREA:
        f = math.ceil(math.sqrt(MIN_SOFT_AREA / min_soft))
        W, H = W * f, H * f
        bx = [v * f for v in bx]
        by = [v * f for v in by]
        bw_ = [v * f for v in bw_]
        bh_ = [v * f for v in bh_]
        cells = [(c[0] * f, c[1] * f, c[2] * f, c[3] * f) + tuple(c[4:]) for c in cells]

    blocks = {i: (bx[i], by[i], bw_[i], bh_[i]) for i in range(n)}
    bcx = np.array([2 * bx[i] + bw_[i] for i in range(n)], float) / 2
    bcy = np.array([2 * by[i] + bh_[i] for i in range(n)], float) / 2

    # ---- terminals ------------------------------------------------------------
    border_cells = [i for i in range(n) if touches(cells[i], W, H)]
    terms = []
    used_pts = set()
    tconn = []  # (term idx, block idx)
    guard = 0
    while len(terms) < R:
        guard += 1
        if guard > 100000:
            raise GenError("terminals")
        i = border_cells[int(rng.integers(len(border_cells)))]
        t = touches(cells[i], W, H)
        e = t[int(rng.integers(len(t)))]
        x0, y0, x1, y1 = cells[i][:4]
        if e in ("LEFT", "RIGHT"):
            px = 0 if e == "LEFT" else W
            py = int(rng.integers(y0, y1 + 1))
        else:
            py = 0 if e == "BOTTOM" else H
            px = int(rng.integers(x0, x1 + 1))
        if (px, py) in used_pts:
            continue
        used_pts.add((px, py))
        ti = len(terms)
        terms.append((px, py))
        tconn.append((ti, i))
        if rng.random() < 0.4:
            dist = np.abs(bcx - px) + np.abs(bcy - py)
            dist[i] = np.inf
            j = int(np.argmin(dist))
            tconn.append((ti, j))
    tpos = {i: terms[i] for i in range(R)}

    # ---- block-block connections --------------------------------------------
    m_bb = int(round(1.75 * n))   # average block-block degree 3.5
    m_long = int(round(0.10 * m_bb))
    m_local = m_bb - m_long
    D = np.abs(bcx[:, None] - bcx[None, :]) + np.abs(bcy[:, None] - bcy[None, :])
    np.fill_diagonal(D, np.inf)
    K = 6
    knn = np.argsort(D, axis=1, kind="stable")[:, :K]
    cand = []
    seen = set()
    for i in range(n):
        for r_, j in enumerate(knn[i]):
            a, b = (i, int(j)) if i < j else (int(j), i)
            if (a, b) not in seen:
                seen.add((a, b))
                cand.append((a, b, r_))
    cand.sort()
    pw = np.array([1.0 / (1 + r_) for _, _, r_ in cand])
    pw /= pw.sum()
    m_local = min(m_local, len(cand))
    pick = rng.choice(len(cand), size=m_local, replace=False, p=pw)
    pairs = [(cand[k][0], cand[k][1]) for k in sorted(pick)]
    pairset = set(pairs)
    # every block gets at least one connection
    deg = [0] * n
    for a, b in pairs:
        deg[a] += 1
        deg[b] += 1
    for i in range(n):
        if deg[i] == 0:
            j = int(knn[i][0])
            p = (min(i, j), max(i, j))
            if p not in pairset:
                pairset.add(p)
                pairs.append(p)
                deg[i] += 1
                deg[j] += 1
    medD = float(np.median(D[np.isfinite(D)]))
    guard = 0
    added = 0
    while added < m_long:
        guard += 1
        if guard > 100000:
            break
        a, b = int(rng.integers(n)), int(rng.integers(n))
        if a == b:
            continue
        p = (min(a, b), max(a, b))
        if p in pairset or D[a, b] < medD:
            continue
        pairset.add(p)
        pairs.append(p)
        added += 1

    def weight():
        if rng.random() < 0.6:
            return int(rng.integers(1, 4))
        return int(rng.integers(1, 11))

    conns = [(("b", a), ("b", b), weight()) for a, b in pairs]
    conns += [(("t", t), ("b", b), weight()) for t, b in tconn]
    conns = [conns[j] for j in rng.permutation(len(conns))]

    # ---- names ------------------------------------------------------------------
    perm = rng.permutation(n)
    bname = {i: f"b{int(perm[i]) + 1}" for i in range(n)}
    tname = {t: f"t{t + 1}" for t in range(R)}

    v2 = wl2(blocks, tpos, conns)
    thr = fmt_half(v2)  # WL is already a multiple of 0.5

    case = dict(W=W, H=H, n=n, kind=kind, blocks=blocks, cells=cells,
                cmib=cmib, groups=groups, bnd=bnd, terms=terms, conns=conns,
                bname=bname, tname=tname, wl2=v2, thr=thr,
                pinwheels=dis.pinwheels, mib_k=mib_k, target=target,
                big_cells=big_cells)
    self_check(case, args)
    return case, rng


def self_check(c, args):
    W, H, n = c["W"], c["H"], c["n"]
    blocks, kind, cells = c["blocks"], c["kind"], c["cells"]
    if max(W, H) > COORD_MAX:
        raise GenError("outline exceeds the coordinate bound")
    if len(c["conns"]) > CONN_MAX:
        raise GenError("too many connections")
    assert kind.count("S") == args.soft
    assert kind.count("F") == args.fixed
    assert kind.count("P") == args.preplaced
    tot = 0
    for i in range(n):
        x, y, w, h = blocks[i]
        x0, y0, x1, y1 = cells[i][:4]
        assert w > 0 and h > 0
        assert x0 <= x and y0 <= y and x + w <= x1 and y + h <= y1, "outside cell"
        assert 0 <= x and 0 <= y and x + w <= W and y + h <= H, "outside outline"
        if kind[i] == "P":
            assert x > 0 and y > 0 and x + w < W and y + h < H, "preplaced touches outline edge"
        if kind[i] == "S":
            assert ASPECT_R * h >= w and ASPECT_R * w >= h, "aspect"
            assert w * h >= MIN_SOFT_AREA, "soft area"
        tot += w * h
    # overlap (blocks are inside disjoint cells, but check directly)
    order = sorted(range(n), key=lambda i: blocks[i][0])
    for a in range(n):
        i = order[a]
        xi, yi, wi, hi = blocks[i]
        for b in range(a + 1, n):
            j = order[b]
            xj, yj, wj, hj = blocks[j]
            if xj >= xi + wi:
                break
            assert not (max(yi, yj) < min(yi + hi, yj + hj)), "overlap"
    ws_act = 1 - tot / (W * H)
    assert abs(ws_act - args.whitespace) <= 0.005, f"whitespace {ws_act}"
    for g, mem in c["groups"].items():
        dims = {blocks[i][2:] for i in mem}
        assert len(dims) == 1, "MIB dims"
        assert all(kind[i] == "S" for i in mem)
    assert len(c["groups"]) == args.mib_groups
    nb = 0
    nc = 0
    for i, edges in c["bnd"].items():
        x, y, w, h = blocks[i]
        assert kind[i] in ("S", "F")
        for e in edges:
            assert {"LEFT": x == 0, "RIGHT": x + w == W,
                    "BOTTOM": y == 0, "TOP": y + h == H}[e], "boundary"
        nb += 1
        nc += len(edges) == 2
    assert nb == args.boundary and nc == args.corners
    for x, y in c["terms"]:
        assert (x in (0, W) and 0 <= y <= H) or (y in (0, H) and 0 <= x <= W)
    seen = set()
    for a, b, w in c["conns"]:
        assert 1 <= w <= 10
        assert a != b and not (a[0] == "t" and b[0] == "t")
        key = tuple(sorted([a, b]))
        assert key not in seen, "duplicate connection"
        seen.add(key)


def write_case(c, args, rng):
    n, W, H = c["n"], c["W"], c["H"]
    blocks, kind, bname, tname = c["blocks"], c["kind"], c["bname"], c["tname"]

    def nm(e):
        return bname[e[1]] if e[0] == "b" else tname[e[1]]

    lines = [f"OUTLINE {W} {H}", f"WLTHRESHOLD {c['thr']}",
             f"ASPECTRATIO {ASPECT_R}", f"BLOCKS {n}"]
    for i in [int(j) for j in rng.permutation(n)]:
        x, y, w, h = blocks[i]
        if kind[i] == "S":
            lines.append(f"SOFT {bname[i]} {w * h}")
        elif kind[i] == "F":
            lines.append(f"FIXED {bname[i]} {w} {h}")
        else:
            lines.append(f"PREPLACED {bname[i]} {w} {h} {x} {y}")
    lines.append(f"TERMINALS {len(c['terms'])}")
    for t, (x, y) in enumerate(c["terms"]):
        lines.append(f"TERMINAL {tname[t]} {x} {y}")
    lines.append(f"CONNECTIONS {len(c['conns'])}")
    for a, b, w in c["conns"]:
        if rng.random() < 0.5:
            a, b = b, a
        if a[0] == "t" and rng.random() < 0.7:
            a, b = b, a
        lines.append(f"CONN {nm(a)} {nm(b)} {w}")
    gl = [c["groups"][g] for g in sorted(c["groups"])]
    gl = [gl[j] for j in rng.permutation(len(gl))]
    lines.append(f"MIBS {len(gl)}")
    for mem in gl:
        mem = [mem[j] for j in rng.permutation(len(mem))]
        lines.append(f"MIB {len(mem)} " + " ".join(bname[i] for i in mem))
    bl = sorted(c["bnd"].items())
    bl = [bl[j] for j in rng.permutation(len(bl))]
    lines.append(f"BOUNDARIES {len(bl)}")
    for i, edges in bl:
        lines.append(f"BOUNDARY {bname[i]} " + " ".join(edges))

    out = [f"WIRELENGTH {fmt_half(c['wl2'])}", f"BLOCKS {n}"]
    for i in sorted(range(n), key=lambda i: int(bname[i][1:])):
        x, y, w, h = blocks[i]
        out.append(f"{bname[i]} {x} {y} {w} {h}")

    d = os.path.dirname(args.out)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(args.out + ".in", "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(args.out + ".construction.out", "w") as f:
        f.write("\n".join(out) + "\n")


def stats_line(c, args):
    n, W, H = c["n"], c["W"], c["H"]
    blocks, kind = c["blocks"], c["kind"]
    tot = sum(b[2] * b[3] for b in blocks.values())
    soft = [i for i in range(n) if kind[i] == "S"]
    nc = sum(1 for e in c["bnd"].values() if len(e) == 2)
    return (f"{os.path.basename(args.out)}: n={n} W={W} H={H} "
            f"ws={1 - tot / (W * H):.4f} conns={len(c['conns'])} "
            f"terms={len(c['terms'])} mib={len(c['groups'])}"
            f"(k={','.join(str(len(c['groups'][g])) for g in sorted(c['groups']))}) "
            f"bnd={len(c['bnd'])} corners={nc} pinwheels={c['pinwheels']} "
            f"minsoft={min(blocks[i][2] * blocks[i][3] for i in soft)} "
            f"maxsoftasp={max(aspect(blocks[i][2], blocks[i][3]) for i in soft):.3f} "
            f"maxasp={max(aspect(b[2], b[3]) for b in blocks.values()):.3f} "
            f"preA={sum(blocks[i][2] * blocks[i][3] for i in range(n) if kind[i] == 'P') / (W * H):.3f} "
            f"cellratio={max(b[2] * b[3] for b in blocks.values()) / min(b[2] * b[3] for b in blocks.values()):.1f} "
            f"bigpre={sum(1 for i in c['big_cells'] if kind[i] == 'P')}/{len(c['big_cells'])} "
            f"WL={fmt_half(c['wl2'])}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preset", choices=sorted(PRESETS))
    ap.add_argument("--blocks", type=int)
    ap.add_argument("--soft", type=int)
    ap.add_argument("--fixed", type=int)
    ap.add_argument("--preplaced", type=int)
    ap.add_argument("--terminals", type=int)
    ap.add_argument("--mib-groups", type=int)
    ap.add_argument("--boundary", type=int)
    ap.add_argument("--corners", type=int)
    ap.add_argument("--whitespace", type=float)
    st = ap.add_argument_group("structure knobs (optional; defaults reproduce the preset output exactly)")
    st.add_argument("--size-tail", type=float, metavar="SIGMA",
                    help="lognormal(0, SIGMA) area weights per cell: heavy-tailed block sizes")
    st.add_argument("--big-cells", type=int, metavar="K",
                    help="K non-MIB cells share the area fraction --big-frac")
    st.add_argument("--big-frac", type=float, default=0.3,
                    help="total area fraction of the --big-cells cells (default 0.3)")
    st.add_argument("--big-kind", choices=["any", "pre", "fixed"], default="any",
                    help="kind of the big cells: pre = preplaced when interior, else fixed; "
                         "fixed = fixed-shape; any = ordinary assignment")
    st.add_argument("--pre-arrange", choices=["spread", "line", "cross", "ring", "largest"],
                    default="spread",
                    help="interior preplaced cells: spread (default), along a mid line, a cross, "
                         "a ring around the centre, or the largest cells")
    st.add_argument("--pre-fill", action="store_true",
                    help="preplaced blocks fill their cells (no slack around obstacles)")
    st.add_argument("--outline-aspect", type=float, metavar="W/H",
                    help="fixed outline aspect instead of uniform(0.75, 1.33)")
    st.add_argument("--mib-k-min", type=int)
    st.add_argument("--mib-k-max", type=int, help="MIB group size range (default 2..4)")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fields = ["soft", "fixed", "preplaced", "terminals", "mib_groups",
              "boundary", "corners", "whitespace"]
    if args.preset:
        for f, v in zip(fields, PRESETS[args.preset]):
            if getattr(args, f) is None:
                setattr(args, f, v)
    missing = [f for f in fields if getattr(args, f) is None]
    if missing:
        ap.error("missing: " + ", ".join("--" + m.replace("_", "-") for m in missing))
    n = args.soft + args.fixed + args.preplaced
    if args.blocks is None:
        args.blocks = n
    if args.blocks != n:
        ap.error("--blocks must equal --soft + --fixed + --preplaced")
    if not 0.0 < args.whitespace < 0.5:
        ap.error("--whitespace must be in (0, 0.5)")
    if args.corners > 4:
        # A two-edge (corner) constraint forces the block to contain an
        # outline corner point; at most four blocks can do that.
        print(f"warning: --corners {args.corners} > 4 is infeasible; "
              f"clamped to 4 (the other {args.corners - 4} become single-edge "
              f"BOUNDARY records)", file=sys.stderr)
        args.corners = 4
    if args.corners > args.boundary:
        ap.error("--corners must be <= --boundary")
    if n > 200 or args.terminals > 100 or args.mib_groups > 20 or args.preplaced > 100 \
            or args.boundary > 100:
        ap.error("outside the docs/design.md section 3 bounds (n<=200, terminals<=100, "
                 "MIB groups<=20, preplaced<=100, boundaries<=100)")
    if args.big_cells is not None and not 0.0 < args.big_frac < 0.9:
        ap.error("--big-frac must be in (0, 0.9)")
    if args.size_tail is not None and args.size_tail < 0:
        ap.error("--size-tail must be >= 0")

    for attempt in range(200):
        try:
            case, rng = generate(args, attempt)
            break
        except GenError:
            continue
    else:
        sys.exit("error: generation failed after 200 attempts")
    write_case(case, args, rng)
    print(stats_line(case, args))


if __name__ == "__main__":
    main()
