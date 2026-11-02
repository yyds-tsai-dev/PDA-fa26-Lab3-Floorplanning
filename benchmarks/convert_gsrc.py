#!/usr/bin/env python3
"""Convert GSRC n100/n200/n300 (UCSC/UCLA bookshelf .blocks/.nets/.pl) to Lab3.

Usage: python3 -I benchmarks/convert_gsrc.py <gsrc_dir> <out_dir>
  <gsrc_dir> holds HARD/nXXX.{blocks,nets,pl} and SOFT/nXXX.{blocks,nets,pl}
  as downloaded from http://vlsicad.eecs.umich.edu/BK/GSRCbench/.

For each n and each of HARD/SOFT, writes nXXX_{hard,soft}_ws{15,10}.in.
See benchmarks/README.md for the rules.
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lab3fmt  # noqa: E402

CASES = ('n100', 'n200', 'n300')
WHITESPACE = ((15, 0.15), (10, 0.10))


def strip_comments(path):
    with open(path) as f:
        for line in f:
            line = line.split('#', 1)[0].strip()
            if line:
                yield line


def read_blocks(path):
    """Return (blocks, terminals): blocks = [(name, kind, w, h, area)]."""
    blocks, terms = [], []
    for line in strip_comments(path):
        p = line.replace(',', ' ').replace('(', ' ').replace(')', ' ').split()
        if p[0] in ('UCSC', 'UCLA') or ':' in line:
            continue
        if p[1] == 'terminal':
            terms.append(p[0])
        elif p[1] == 'hardrectilinear':
            nv = int(p[2])
            if nv != 4:
                raise ValueError(f'{path}: non-rectangular block {p[0]}')
            xs = [int(v) for v in p[3::2]]
            ys = [int(v) for v in p[4::2]]
            w, h = max(xs) - min(xs), max(ys) - min(ys)
            blocks.append((p[0], 'hard', w, h, w * h))
        elif p[1] == 'softrectangular':
            a = float(p[2])
            if a != int(a):
                raise ValueError(f'{path}: non-integer soft area {p[0]}')
            blocks.append((p[0], 'soft', None, None, int(a)))
        else:
            raise ValueError(f'{path}: unknown line {line}')
    return blocks, terms


def read_nets(path):
    nets, cur, nnets = [], None, None
    for line in strip_comments(path):
        p = line.split()
        if p[0] in ('UCSC', 'UCLA'):
            continue
        if p[0] == 'NumNets':
            nnets = int(p[-1])
        elif p[0] == 'NumPins':
            continue
        elif p[0] == 'NetDegree':
            cur = []
            nets.append(cur)
        else:
            cur.append(p[0])
    if nnets is not None and nnets != len(nets):
        raise ValueError(f'{path}: NumNets {nnets} != {len(nets)}')
    return nets


def read_pl(path):
    pos = {}
    for line in strip_comments(path):
        p = line.split()
        if p[0] in ('UCSC', 'UCLA'):
            continue
        pos[p[0]] = (float(p[1]), float(p[2]))
    return pos


def convert(blocks, terms, nets, pl, scale, ws, soft):
    s = scale
    total = sum(b[4] for b in blocks) * s * s
    W = H = lab3fmt.square_side(total, ws)
    c = lab3fmt.Case(W, H)
    for n, kind, w, h, a in blocks:
        if soft:
            c.blocks.append(('SOFT', n, a * s * s))
        else:
            c.blocks.append(('FIXED', n, w * s, h * s))
    # pad frame = bounding box of the pads in the .pl file
    px = [pl[t][0] for t in terms]
    py = [pl[t][1] for t in terms]
    fx0, fx1, fy0, fy1 = min(px), max(px), min(py), max(py)
    for t in terms:
        x, y = pl[t]
        sx = (x - fx0) / (fx1 - fx0) * W
        sy = (y - fy0) / (fy1 - fy0) * H
        c.terminals.append((t,) + lab3fmt.project_to_border(sx, sy, W, H))
    bnames = {b[0] for b in blocks}
    tnames = set(terms)
    is_block = bnames.__contains__
    pairs = []
    for net in nets:
        for p in net:
            if p not in bnames and p not in tnames:
                raise ValueError(f'unknown pin {p}')
        pairs += lab3fmt.net_to_pairs(net, is_block)
    c.conns = lab3fmt.dedupe_connections(pairs, is_block)
    return c, (fx0, fy0, fx1, fy1)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    for n in CASES:
        data = {}
        for d in ('HARD', 'SOFT'):
            blocks, terms = read_blocks(os.path.join(src, d, n + '.blocks'))
            nets = read_nets(os.path.join(src, d, n + '.nets'))
            pl = read_pl(os.path.join(src, d, n + '.pl'))
            data[d] = (blocks, terms, nets, pl)
        # one integer linear scale per n, shared by the hard and soft versions,
        # so that every soft area is >= 10000
        min_soft = min(b[4] for b in data['SOFT'][0])
        scale = int(math.ceil(math.sqrt(10000.0 / min_soft)))
        for d in ('HARD', 'SOFT'):
            blocks, terms, nets, pl = data[d]
            for tag, ws in WHITESPACE:
                c, frame = convert(blocks, terms, nets, pl, scale, ws, d == 'SOFT')
                fx0, fy0, fx1, fy1 = frame
                interior = sum(1 for t in terms
                               if fx0 < pl[t][0] < fx1 and fy0 < pl[t][1] < fy1)
                err = lab3fmt.sanity(c)
                name = f'{n}_{d.lower()}_ws{tag}'
                if err:
                    sys.exit(f'{name}: ' + '; '.join(err[:10]))
                path = os.path.join(out, name + '.in')
                c.write(path)
                print(f'{path}: n={len(c.blocks)} scale={scale} W=H={c.W} '
                      f'terms={len(c.terminals)} (interior pads {interior}) '
                      f'conns={len(c.conns)} nets={len(nets)} ws={c.whitespace():.4f}')


if __name__ == '__main__':
    main()
