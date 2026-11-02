#!/usr/bin/env python3
"""Render a Lab3 floorplan as SVG (Python 3.10 stdlib only).

    python3 tools/plot.py <input> <output> <out.svg> [--nets]

Draws the outline, blocks coloured by kind (soft / fixed-shape / preplaced),
MIB groups as a shared dashed outline colour, boundary constraints as thick
bars on the required block edges (green = satisfied, red = not), terminals as
dots, and optionally connections as faint lines (--nets). Block names are
shown when n <= 60. The output file is read leniently so that infeasible
outputs can still be inspected; unparsable lines are skipped with a warning.
"""

import argparse
import sys
from xml.sax.saxutils import escape

KIND_FILL = {"SOFT": "#9ecae1", "FIXED": "#fdae6b", "PREPLACED": "#bdbdbd"}
KIND_LABEL = {"SOFT": "soft", "FIXED": "fixed-shape", "PREPLACED": "preplaced"}
MIB_COLORS = ["#7b3294", "#008837", "#e66101", "#0571b0", "#ca0020",
              "#5e3c99", "#1b7837", "#b35806", "#2166ac", "#d6604d",
              "#542788", "#01665e"]


def parse_input(path):
    with open(path) as f:
        tok = f.read().split()
    pos = 0

    def nxt():
        nonlocal pos
        if pos >= len(tok):
            raise ValueError("unexpected end of input file")
        pos += 1
        return tok[pos - 1]

    def kw(k):
        t = nxt()
        if t != k:
            raise ValueError(f"expected {k!r}, got {t!r}")

    inst = {"blocks": {}, "order": [], "terminals": {}, "conns": [],
            "mibs": [], "bounds": {}}
    kw("OUTLINE")
    inst["W"], inst["H"] = int(nxt()), int(nxt())
    kw("WLTHRESHOLD")
    inst["T"] = nxt()
    kw("ASPECTRATIO")
    inst["R"] = int(nxt())
    kw("BLOCKS")
    for _ in range(int(nxt())):
        kind, name = nxt(), nxt()
        b = {"kind": kind}
        if kind == "SOFT":
            b["area"] = int(nxt())
        elif kind in ("FIXED", "PREPLACED"):
            b["w"], b["h"] = int(nxt()), int(nxt())
            if kind == "PREPLACED":
                b["x"], b["y"] = int(nxt()), int(nxt())
        else:
            raise ValueError(f"unknown block kind {kind!r}")
        inst["blocks"][name] = b
        inst["order"].append(name)
    kw("TERMINALS")
    for _ in range(int(nxt())):
        kw("TERMINAL")
        name = nxt()
        inst["terminals"][name] = (int(nxt()), int(nxt()))
    kw("CONNECTIONS")
    for _ in range(int(nxt())):
        kw("CONN")
        inst["conns"].append((nxt(), nxt(), int(nxt())))
    kw("MIBS")
    for _ in range(int(nxt())):
        kw("MIB")
        k = int(nxt())
        inst["mibs"].append([nxt() for _ in range(k)])
    kw("BOUNDARIES")
    edges = {"LEFT", "RIGHT", "TOP", "BOTTOM"}
    for _ in range(int(nxt())):
        kw("BOUNDARY")
        name = nxt()
        es = [nxt()]
        if pos < len(tok) and tok[pos] in edges:
            es.append(nxt())
        inst["bounds"][name] = es
    return inst


def parse_output(path, inst):
    rects = {}
    wl = None
    with open(path) as f:
        lines = [ln.split() for ln in f]
    for i, t in enumerate(lines, 1):
        if not t:
            continue
        if t[0] == "WIRELENGTH" and len(t) == 2:
            wl = t[1]
            continue
        if t[0] == "BLOCKS" and len(t) == 2:
            continue
        if len(t) == 5 and t[0] in inst["blocks"]:
            try:
                rects[t[0]] = tuple(int(v) for v in t[1:])
                continue
            except ValueError:
                pass
        print(f"plot.py: warning: skipping output line {i}: {' '.join(t)}",
              file=sys.stderr)
    return rects, wl


def boundary_ok(edge, r, W, H):
    x, y, w, h = r
    return {"LEFT": x == 0, "RIGHT": x + w == W,
            "BOTTOM": y == 0, "TOP": y + h == H}[edge]


def render(inst, rects, wl, show_nets):
    W, H = inst["W"], inst["H"]
    n = len(inst["order"])
    draw_w = 900.0
    s = draw_w / W
    draw_h = H * s
    margin = 30.0
    legend_h = 70.0
    title_h = 28.0
    svg_w = draw_w + 2 * margin
    svg_h = draw_h + 2 * margin + legend_h + title_h
    oy = margin + title_h

    def X(x):
        return margin + x * s

    def Y(y):  # flip: input y grows upward
        return oy + (H - y) * s

    font = max(7.0, min(13.0, 0.35 * s * min(W, H) / max(4, n ** 0.5)))
    out = []
    out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w:.0f}" '
               f'height="{svg_h:.0f}" viewBox="0 0 {svg_w:.1f} {svg_h:.1f}" '
               f'font-family="Helvetica, Arial, sans-serif">')
    out.append('<defs><pattern id="hatch" width="6" height="6" '
               'patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
               '<line x1="0" y1="0" x2="0" y2="6" stroke="#555" '
               'stroke-width="1" stroke-opacity="0.35"/></pattern></defs>')
    out.append(f'<rect x="0" y="0" width="{svg_w:.1f}" height="{svg_h:.1f}" fill="#ffffff"/>')
    title = f"Lab3 floorplan  {W}x{H}  blocks {len(rects)}/{n}"
    if wl is not None:
        title += f"  reported WL {wl}"
    title += f"  T {inst['T']}"
    out.append(f'<text x="{margin}" y="{margin + 4}" font-size="14" '
               f'font-weight="bold">{escape(title)}</text>')

    # Outline.
    out.append(f'<rect x="{X(0):.2f}" y="{Y(H):.2f}" width="{W * s:.2f}" '
               f'height="{H * s:.2f}" fill="#fafafa" stroke="#000" stroke-width="2"/>')

    mib_of = {}
    for gi, grp in enumerate(inst["mibs"]):
        for name in grp:
            mib_of[name] = gi

    # Blocks.
    for name in inst["order"]:
        if name not in rects:
            continue
        b = inst["blocks"][name]
        x, y, w, h = rects[name]
        fill = KIND_FILL.get(b["kind"], "#ffffff")
        stroke, sw, dash = "#333", 1.0, ""
        if name in mib_of:
            stroke = MIB_COLORS[mib_of[name] % len(MIB_COLORS)]
            sw, dash = 2.5, ' stroke-dasharray="6,3"'
        out.append(f'<rect x="{X(x):.2f}" y="{Y(y + h):.2f}" width="{w * s:.2f}" '
                   f'height="{h * s:.2f}" fill="{fill}" fill-opacity="0.85" '
                   f'stroke="{stroke}" stroke-width="{sw}"{dash}>'
                   f'<title>{escape(name)} {b["kind"]} ({x},{y}) {w}x{h}</title></rect>')
        if b["kind"] == "PREPLACED":
            out.append(f'<rect x="{X(x):.2f}" y="{Y(y + h):.2f}" width="{w * s:.2f}" '
                       f'height="{h * s:.2f}" fill="url(#hatch)" stroke="none"/>')

    # Boundary constraints: bar along each required edge of the block.
    for name, es in inst["bounds"].items():
        if name not in rects:
            continue
        r = rects[name]
        x, y, w, h = r
        for e in es:
            col = "#1a9641" if boundary_ok(e, r, W, H) else "#d7191c"
            seg = {"LEFT": (x, y, x, y + h), "RIGHT": (x + w, y, x + w, y + h),
                   "BOTTOM": (x, y, x + w, y), "TOP": (x, y + h, x + w, y + h)}[e]
            out.append(f'<line x1="{X(seg[0]):.2f}" y1="{Y(seg[1]):.2f}" '
                       f'x2="{X(seg[2]):.2f}" y2="{Y(seg[3]):.2f}" stroke="{col}" '
                       f'stroke-width="5" stroke-linecap="square"/>')

    # Connections.
    if show_nets:
        def centre(nm):
            if nm in rects:
                x, y, w, h = rects[nm]
                return x + w / 2, y + h / 2
            if nm in inst["terminals"]:
                return inst["terminals"][nm]
            return None
        for a, b, wgt in inst["conns"]:
            pa, pb = centre(a), centre(b)
            if pa is None or pb is None:
                continue
            op = 0.12 + 0.05 * wgt
            out.append(f'<line x1="{X(pa[0]):.2f}" y1="{Y(pa[1]):.2f}" '
                       f'x2="{X(pb[0]):.2f}" y2="{Y(pb[1]):.2f}" stroke="#444" '
                       f'stroke-opacity="{op:.2f}" stroke-width="{0.5 + 0.15 * wgt:.2f}"/>')

    # Terminals.
    for name, (tx, ty) in inst["terminals"].items():
        out.append(f'<circle cx="{X(tx):.2f}" cy="{Y(ty):.2f}" r="4" fill="#000">'
                   f'<title>{escape(name)} ({tx},{ty})</title></circle>')
        if n <= 60:
            out.append(f'<text x="{X(tx) + 5:.2f}" y="{Y(ty) - 5:.2f}" '
                       f'font-size="{font * 0.85:.1f}" fill="#000">{escape(name)}</text>')

    # Labels.
    if n <= 60:
        for name in inst["order"]:
            if name not in rects:
                continue
            x, y, w, h = rects[name]
            out.append(f'<text x="{X(x + w / 2):.2f}" y="{Y(y + h / 2):.2f}" '
                       f'font-size="{font:.1f}" text-anchor="middle" '
                       f'dominant-baseline="central">{escape(name)}</text>')

    # Legend.
    ly = oy + draw_h + margin * 0.6
    lx = margin
    items = [("SOFT", None), ("FIXED", None), ("PREPLACED", None)]
    for kind, _ in items:
        out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="18" height="12" '
                   f'fill="{KIND_FILL[kind]}" stroke="#333"/>')
        if kind == "PREPLACED":
            out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="18" height="12" '
                       f'fill="url(#hatch)"/>')
        out.append(f'<text x="{lx + 24:.1f}" y="{ly + 10:.1f}" font-size="12">'
                   f'{KIND_LABEL[kind]}</text>')
        lx += 120
    out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="18" height="12" fill="none" '
               f'stroke="{MIB_COLORS[0]}" stroke-width="2.5" stroke-dasharray="6,3"/>')
    out.append(f'<text x="{lx + 24:.1f}" y="{ly + 10:.1f}" font-size="12">'
               f'MIB group (one colour per group)</text>')
    lx = margin
    ly += 22
    for col, txt in (("#1a9641", "boundary edge satisfied"),
                     ("#d7191c", "boundary edge violated")):
        out.append(f'<line x1="{lx:.1f}" y1="{ly + 6:.1f}" x2="{lx + 18:.1f}" '
                   f'y2="{ly + 6:.1f}" stroke="{col}" stroke-width="5"/>')
        out.append(f'<text x="{lx + 24:.1f}" y="{ly + 10:.1f}" font-size="12">{txt}</text>')
        lx += 200
    out.append(f'<circle cx="{lx + 9:.1f}" cy="{ly + 6:.1f}" r="4" fill="#000"/>')
    out.append(f'<text x="{lx + 24:.1f}" y="{ly + 10:.1f}" font-size="12">terminal</text>')
    lx += 100
    if show_nets:
        out.append(f'<line x1="{lx:.1f}" y1="{ly + 6:.1f}" x2="{lx + 18:.1f}" '
                   f'y2="{ly + 6:.1f}" stroke="#444" stroke-opacity="0.4"/>')
        out.append(f'<text x="{lx + 24:.1f}" y="{ly + 10:.1f}" font-size="12">'
                   f'connection (opacity ~ weight)</text>')

    out.append("</svg>")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description="Render a Lab3 floorplan as SVG.")
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("svg")
    ap.add_argument("--nets", action="store_true", help="draw connections")
    a = ap.parse_args()
    try:
        inst = parse_input(a.input)
    except (OSError, ValueError) as e:
        print(f"plot.py: input error: {e}", file=sys.stderr)
        return 2
    try:
        rects, wl = parse_output(a.output, inst)
    except OSError as e:
        print(f"plot.py: output error: {e}", file=sys.stderr)
        return 2
    with open(a.svg, "w") as f:
        f.write(render(inst, rects, wl, a.nets))
    return 0


if __name__ == "__main__":
    sys.exit(main())
