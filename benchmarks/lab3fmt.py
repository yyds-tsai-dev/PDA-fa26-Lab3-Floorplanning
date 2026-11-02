"""Shared helpers for the benchmark converters: Lab3 writer, terminal projection,
multi-pin net decomposition, connection dedupe, and input sanity checks.

Lab3 input format: docs/design.md section 3. Python 3.10, stdlib only.
"""

import math

ASPECT_R = 3
WEIGHT_MIN, WEIGHT_MAX = 1, 10


def project_to_border(x, y, W, H):
    """Clamp (x, y) into [0,W]x[0,H], then snap to the nearest outline side.

    Returns integer coordinates on the border."""
    x = min(max(x, 0.0), float(W))
    y = min(max(y, 0.0), float(H))
    d = [(x, 0), (W - x, 1), (y, 2), (H - y, 3)]  # left, right, bottom, top
    side = min(d)[1]
    if side == 0:
        return 0, min(max(int(round(y)), 0), H)
    if side == 1:
        return W, min(max(int(round(y)), 0), H)
    if side == 2:
        return min(max(int(round(x)), 0), W), 0
    return min(max(int(round(x)), 0), W), H


def net_to_pairs(pins, is_block):
    """Decompose one multi-pin net into weighted 2-pin pairs.

    Degree <= 5: clique, weight 1 per pair. Degree > 5: star from every pin to
    the net's first block pin, weight 1. Terminal-terminal pairs and self pairs
    are skipped. A net without any block pin yields nothing."""
    pins = list(dict.fromkeys(pins))  # drop repeated pins, keep order
    blocks = [p for p in pins if is_block(p)]
    if not blocks:
        return []
    out = []
    if len(pins) <= 5:
        for i in range(len(pins)):
            for j in range(i + 1, len(pins)):
                a, b = pins[i], pins[j]
                if is_block(a) or is_block(b):
                    out.append((a, b, 1))
    else:
        c = blocks[0]
        for p in pins:
            if p != c:
                out.append((c, p, 1))
    return out


def dedupe_connections(pairs, is_block):
    """Sum weights of identical unordered pairs, drop self and terminal-terminal
    pairs, round and clamp weights to [1, 10]. Returns a sorted list."""
    acc = {}
    for a, b, w in pairs:
        if a == b or not (is_block(a) or is_block(b)):
            continue
        key = (a, b) if a < b else (b, a)
        acc[key] = acc.get(key, 0) + w
    out = []
    for (a, b), w in sorted(acc.items()):
        wi = min(max(int(round(w)), WEIGHT_MIN), WEIGHT_MAX)
        # block first for readability
        if not is_block(a):
            a, b = b, a
        out.append((a, b, wi))
    return out


class Case:
    """In-memory Lab3 case. blocks: list of tuples
    ('SOFT', name, area) | ('FIXED', name, w, h) | ('PREPLACED', name, w, h, x, y)."""

    def __init__(self, W, H, T=0.0, R=ASPECT_R):
        self.W, self.H, self.T, self.R = int(W), int(H), float(T), int(R)
        self.blocks = []
        self.terminals = []      # (name, x, y)
        self.conns = []          # (a, b, w)
        self.mibs = []           # [names]
        self.boundaries = []     # (name, [edges])

    def block_area(self, b):
        return b[2] if b[0] == 'SOFT' else b[2] * b[3]

    def total_area(self):
        return sum(self.block_area(b) for b in self.blocks)

    def whitespace(self):
        return 1.0 - self.total_area() / (self.W * self.H)

    def write(self, path):
        L = [f'OUTLINE {self.W} {self.H}',
             f'WLTHRESHOLD {self.T:.1f}',
             f'ASPECTRATIO {self.R}',
             f'BLOCKS {len(self.blocks)}']
        for b in self.blocks:
            L.append(' '.join(str(v) for v in b))
        L.append(f'TERMINALS {len(self.terminals)}')
        for n, x, y in self.terminals:
            L.append(f'TERMINAL {n} {x} {y}')
        L.append(f'CONNECTIONS {len(self.conns)}')
        for a, b, w in self.conns:
            L.append(f'CONN {a} {b} {w}')
        L.append(f'MIBS {len(self.mibs)}')
        for g in self.mibs:
            L.append(f'MIB {len(g)} ' + ' '.join(g))
        L.append(f'BOUNDARIES {len(self.boundaries)}')
        for n, e in self.boundaries:
            L.append(f'BOUNDARY {n} ' + ' '.join(e))
        with open(path, 'w') as f:
            f.write('\n'.join(L) + '\n')


def read_case(path):
    """Parse a Lab3 .in file written by Case.write (one record per line) back
    into a Case; used by the sanity checker."""
    with open(path) as f:
        lines = [l.split() for l in f if l.strip()]
    pos = 0

    def take(key, nmin):
        nonlocal pos
        p = lines[pos]
        pos += 1
        if p[0] != key or len(p) < nmin:
            raise ValueError(f'line {pos}: expected {key}, got {p}')
        return p

    W, H = map(int, take('OUTLINE', 3)[1:3])
    T = float(take('WLTHRESHOLD', 2)[1])
    c = Case(W, H, T, int(take('ASPECTRATIO', 2)[1]))
    for _ in range(int(take('BLOCKS', 2)[1])):
        p = lines[pos]
        pos += 1
        want = {'SOFT': 3, 'FIXED': 4, 'PREPLACED': 6}.get(p[0])
        if want is None or len(p) != want:
            raise ValueError(f'line {pos}: bad block {p}')
        c.blocks.append((p[0], p[1]) + tuple(int(v) for v in p[2:]))
    for _ in range(int(take('TERMINALS', 2)[1])):
        p = take('TERMINAL', 4)
        c.terminals.append((p[1], int(p[2]), int(p[3])))
    for _ in range(int(take('CONNECTIONS', 2)[1])):
        p = take('CONN', 4)
        c.conns.append((p[1], p[2], int(p[3])))
    for _ in range(int(take('MIBS', 2)[1])):
        p = take('MIB', 2)
        if len(p) != 2 + int(p[1]):
            raise ValueError(f'line {pos}: MIB count mismatch')
        c.mibs.append(p[2:])
    for _ in range(int(take('BOUNDARIES', 2)[1])):
        p = take('BOUNDARY', 3)
        c.boundaries.append((p[1], p[2:]))
    if pos != len(lines):
        raise ValueError('trailing lines')
    return c


def sanity(c):
    """Return a list of violations of docs/design.md section 3 rules (empty = ok)."""
    err = []
    W, H = c.W, c.H
    if W <= 0 or H <= 0:
        err.append('bad outline')
    if c.T < 0 or (c.T * 2) != int(c.T * 2):
        err.append('bad T')
    names = [b[1] for b in c.blocks] + [t[0] for t in c.terminals]
    if len(set(names)) != len(names):
        err.append('duplicate names')
    kind = {b[1]: b for b in c.blocks}
    tset = {t[0] for t in c.terminals}
    for b in c.blocks:
        if b[0] == 'SOFT' and b[2] < 10000:
            err.append(f'soft area < 10000: {b[1]}')
        if b[0] in ('FIXED', 'PREPLACED') and (b[2] <= 0 or b[3] <= 0):
            err.append(f'bad dims: {b[1]}')
        if b[0] == 'FIXED' and (b[2] > W or b[3] > H):
            err.append(f'fixed block larger than outline: {b[1]}')
    pre = [b for b in c.blocks if b[0] == 'PREPLACED']
    for b in pre:
        _, n, w, h, x, y = b
        if x <= 0 or y <= 0 or x + w >= W or y + h >= H:
            err.append(f'preplaced not strictly inside outline: {n}')
    for i in range(len(pre)):
        for j in range(i + 1, len(pre)):
            _, n1, w1, h1, x1, y1 = pre[i]
            _, n2, w2, h2, x2, y2 = pre[j]
            if max(x1, x2) < min(x1 + w1, x2 + w2) and max(y1, y2) < min(y1 + h1, y2 + h2):
                err.append(f'preplaced overlap: {n1} {n2}')
    if c.total_area() >= W * H:
        err.append('total block area >= outline area')
    for n, x, y in c.terminals:
        on = ((x in (0, W) and 0 <= y <= H) or (y in (0, H) and 0 <= x <= W))
        if not on:
            err.append(f'terminal off border: {n}')
    seen = set()
    for a, b, w in c.conns:
        if a not in kind and b not in kind:
            err.append(f'conn without block: {a} {b}')
        for e in (a, b):
            if e not in kind and e not in tset:
                err.append(f'conn unknown name: {e}')
        if a == b:
            err.append(f'self conn: {a}')
        key = frozenset((a, b))
        if key in seen:
            err.append(f'duplicate conn: {a} {b}')
        seen.add(key)
        if not (WEIGHT_MIN <= w <= WEIGHT_MAX):
            err.append(f'bad weight: {a} {b} {w}')
    inmib = set()
    for g in c.mibs:
        if len(g) < 2:
            err.append('MIB with < 2 members')
        areas = set()
        for n in g:
            b = kind.get(n)
            if b is None or b[0] != 'SOFT':
                err.append(f'MIB member not soft: {n}')
                continue
            areas.add(b[2])
            if n in inmib:
                err.append(f'block in two MIBs: {n}')
            inmib.add(n)
        if len(areas) > 1:
            err.append(f'MIB areas differ: {g}')
    inb = set()
    for n, e in c.boundaries:
        b = kind.get(n)
        if b is None or b[0] == 'PREPLACED':
            err.append(f'boundary on non-soft/fixed: {n}')
        if n in inb:
            err.append(f'block in two BOUNDARY records: {n}')
        inb.add(n)
        if not (1 <= len(e) <= 2) or any(x not in ('LEFT', 'RIGHT', 'TOP', 'BOTTOM') for x in e):
            err.append(f'bad edges: {n} {e}')
        if set(e) in ({'LEFT', 'RIGHT'}, {'TOP', 'BOTTOM'}) or len(set(e)) != len(e):
            err.append(f'opposite/duplicate edges: {n} {e}')
    return err


def square_side(total_area, whitespace):
    return int(math.ceil(math.sqrt(total_area / (1.0 - whitespace))))
