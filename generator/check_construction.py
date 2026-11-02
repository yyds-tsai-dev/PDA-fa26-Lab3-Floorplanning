#!/usr/bin/env python3
"""Independent checker for a Lab3 input/output pair (docs/design.md sections 3-6).

Parses both files from disk, applies the gate in order, then scores WL,
boundary and MIB exactly like the verifier spec. Prints a short report and
exits 0 on a feasible output, 1 otherwise.

usage: check_construction.py <case.in> <case.out> [--quiet]
"""

import sys
from fractions import Fraction


def fail(msg):
    print(f"FAIL: {msg}")
    print("Total 0.00 / 20")
    sys.exit(1)


def parse_input(path):
    tok = open(path).read().split()
    p = 0

    def take(expect=None):
        nonlocal p
        if p >= len(tok):
            raise ValueError("unexpected EOF in input")
        t = tok[p]
        p += 1
        if expect is not None and t != expect:
            raise ValueError(f"input: expected {expect}, got {t}")
        return t

    inp = {}
    take("OUTLINE")
    inp["W"], inp["H"] = int(take()), int(take())
    take("WLTHRESHOLD")
    tstr = take()
    inp["T2"] = int(Fraction(tstr) * 2)
    assert Fraction(tstr) * 2 == inp["T2"], "T not a multiple of 0.5"
    take("ASPECTRATIO")
    inp["R"] = int(take())
    take("BLOCKS")
    n = int(take())
    blocks = {}
    for _ in range(n):
        k = take()
        name = take()
        if name in blocks:
            raise ValueError(f"duplicate block {name}")
        if k == "SOFT":
            blocks[name] = ("S", int(take()))
        elif k == "FIXED":
            blocks[name] = ("F", int(take()), int(take()))
        elif k == "PREPLACED":
            blocks[name] = ("P", int(take()), int(take()), int(take()), int(take()))
        else:
            raise ValueError(f"bad block kind {k}")
    inp["blocks"] = blocks
    take("TERMINALS")
    terms = {}
    for _ in range(int(take())):
        take("TERMINAL")
        name = take()
        if name in terms or name in blocks:
            raise ValueError(f"duplicate name {name}")
        terms[name] = (int(take()), int(take()))
    inp["terms"] = terms
    take("CONNECTIONS")
    conns = []
    for _ in range(int(take())):
        take("CONN")
        conns.append((take(), take(), int(take())))
    inp["conns"] = conns
    take("MIBS")
    mibs = []
    for _ in range(int(take())):
        take("MIB")
        k = int(take())
        mibs.append([take() for _ in range(k)])
    inp["mibs"] = mibs
    take("BOUNDARIES")
    bnds = []
    for _ in range(int(take())):
        take("BOUNDARY")
        name = take()
        edges = []
        while p < len(tok) and tok[p] in ("LEFT", "RIGHT", "TOP", "BOTTOM"):
            edges.append(take())
        bnds.append((name, edges))
    inp["bnds"] = bnds
    if p != len(tok):
        raise ValueError("trailing tokens in input")
    return inp


def validate_input(inp):
    """Input-side guarantees of section 3."""
    W, H, R = inp["W"], inp["H"], inp["R"]
    blocks, terms = inp["blocks"], inp["terms"]
    names = set(blocks) | set(terms)
    assert not names & {"LEFT", "RIGHT", "TOP", "BOTTOM"}, "reserved name"
    for name, b in blocks.items():
        if b[0] == "S":
            assert b[1] >= 10000, f"soft area < 10000: {name}"
        if b[0] == "P":
            _, w, h, x, y = b
            assert x > 0 and y > 0 and x + w < W and y + h < H, \
                f"preplaced {name} not strictly inside outline"
    pre = [(n_, b[1:]) for n_, b in blocks.items() if b[0] == "P"]
    for i in range(len(pre)):
        for j in range(i + 1, len(pre)):
            (w1, h1, x1, y1), (w2, h2, x2, y2) = pre[i][1], pre[j][1]
            assert not (max(x1, x2) < min(x1 + w1, x2 + w2) and
                        max(y1, y2) < min(y1 + h1, y2 + h2)), "preplaced overlap"
    for name, (x, y) in terms.items():
        assert (x in (0, W) and 0 <= y <= H) or (y in (0, H) and 0 <= x <= W), \
            f"terminal off border: {name}"
    seen = set()
    for a, b, w in inp["conns"]:
        assert a in names and b in names, f"unknown conn name {a} {b}"
        assert a != b, "self loop"
        assert not (a in terms and b in terms), "terminal-terminal"
        key = (min(a, b), max(a, b))
        assert key not in seen, f"duplicate conn {key}"
        seen.add(key)
        assert 1 <= w <= 10
    inmib = set()
    for mem in inp["mibs"]:
        assert len(mem) >= 2
        areas = set()
        for m in mem:
            assert blocks[m][0] == "S", f"MIB member not soft: {m}"
            assert m not in inmib, f"block in two MIBs: {m}"
            inmib.add(m)
            areas.add(blocks[m][1])
        assert len(areas) == 1, "MIB areas differ"
    inb = set()
    for name, edges in inp["bnds"]:
        assert blocks[name][0] in ("S", "F"), f"boundary on preplaced {name}"
        assert name not in inb
        inb.add(name)
        assert 1 <= len(edges) <= 2
        if len(edges) == 2:
            assert set(edges) not in ({"LEFT", "RIGHT"}, {"TOP", "BOTTOM"})


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    inp = parse_input(sys.argv[1])
    validate_input(inp)
    W, H, R = inp["W"], inp["H"], inp["R"]
    blocks = inp["blocks"]
    n = len(blocks)

    # 1. parse output
    tok = open(sys.argv[2]).read().split()
    try:
        if tok[0] != "WIRELENGTH":
            fail("parse: header")
        rep = float(tok[1])
        if tok[2] != "BLOCKS" or int(tok[3]) != n:
            fail("parse: block count")
        if len(tok) != 4 + 5 * n:
            fail("parse: token count")
        place = {}
        for i in range(n):
            r = tok[4 + 5 * i: 9 + 5 * i]
            name = r[0]
            # 2. names
            if name not in blocks:
                fail(f"unknown block {name}")
            if name in place:
                fail(f"duplicate block {name}")
            place[name] = tuple(int(v) for v in r[1:])
    except (ValueError, IndexError) as e:
        fail(f"parse: {e}")
    for name, (x, y, w, h) in place.items():
        if w <= 0 or h <= 0:
            fail(f"nonpositive size {name}")
    for name, (x, y, w, h) in place.items():
        b = blocks[name]
        if b[0] == "F" and (w, h) != (b[1], b[2]):
            fail(f"fixed shape {name}")
        if b[0] == "P" and (w, h, x, y) != b[1:]:
            fail(f"preplaced {name}")
    for name, (x, y, w, h) in place.items():
        b = blocks[name]
        if b[0] == "S" and 100 * abs(w * h - b[1]) > b[1]:
            fail(f"area {name}")
    for name, (x, y, w, h) in place.items():
        if blocks[name][0] == "S" and not (R * h >= w and R * w >= h):
            fail(f"aspect {name}")
    for name, (x, y, w, h) in place.items():
        if x < 0 or y < 0 or x + w > W or y + h > H:
            fail(f"outside {name}")
    items = list(place.items())
    for i in range(n):
        a, (x1, y1, w1, h1) = items[i]
        for j in range(i + 1, n):
            b, (x2, y2, w2, h2) = items[j]
            if max(x1, x2) < min(x1 + w1, x2 + w2) and \
                    max(y1, y2) < min(y1 + h1, y2 + h2):
                fail(f"Overlap: {a} {b}")

    def p2(name):
        if name in place:
            x, y, w, h = place[name]
            return 2 * x + w, 2 * y + h
        x, y = inp["terms"][name]
        return 2 * x, 2 * y

    wl2 = 0
    for a, b, c in inp["conns"]:
        (ax, ay), (bx, by) = p2(a), p2(b)
        wl2 += c * (abs(ax - bx) + abs(ay - by))
    wl = wl2 / 2
    T2 = inp["T2"]
    T = T2 / 2
    if abs(rep - wl) > 0.01:
        swl = 0.0
    elif wl2 <= T2:
        swl = 5.0
    elif 20 * wl2 <= 21 * T2:
        swl = 5.0 * (1.05 * T - wl) / (0.05 * T)
    else:
        swl = 0.0
    sat = 0
    for name, edges in inp["bnds"]:
        x, y, w, h = place[name]
        ok = {"LEFT": x == 0, "RIGHT": x + w == W, "BOTTOM": y == 0,
              "TOP": y + h == H}
        sat += all(ok[e] for e in edges)
    p = len(inp["bnds"])
    sb = 3.0 * sat / p if p else 3.0
    msat = sum(1 for mem in inp["mibs"]
               if len({place[m][2:] for m in mem}) == 1)
    q = len(inp["mibs"])
    sm = 2.0 * msat / q if q else 2.0
    tot = 10 + swl + sb + sm
    tot_area = sum(w * h for (_, _, w, h) in place.values())
    if "--quiet" not in sys.argv:
        print(f"Gate PASS; Reported WL {rep:.1f}; Computed WL {wl:.1f}; "
              f"Threshold T {T:.1f}; whitespace {1 - tot_area / (W * H):.4f}")
        print(f"WL {swl:.2f} Boundary {sb:.2f} ({sat}/{p}) MIB {sm:.2f} ({msat}/{q})")
    print(f"Total {tot:.2f} / 20")
    sys.exit(0)


if __name__ == "__main__":
    main()
