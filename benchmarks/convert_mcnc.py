#!/usr/bin/env python3
"""Convert MCNC ami33 / ami49 (Lab2 .block/.nets format) to Lab3 input.

Usage: python3 -I benchmarks/convert_mcnc.py <src_dir_with_ami33_ami49> <out_dir>
  <src_dir> holds ami33/ami33.block, ami33/ami33.nets, ami49/ami49.block, ...
  (in this repo: reference/Lab2).

Outputs ami33_fixed.in, ami49_fixed.in (all FIXED, original outline; ami49_fixed
is the t49 yardstick of docs/design.md section 9) and ami49_soft.in (every block
SOFT with area w*h, same outline). See benchmarks/README.md for the rules.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lab3fmt  # noqa: E402


def read_block(path):
    with open(path, newline=None) as f:
        lines = [l.split() for l in f.read().replace('\r', '').split('\n')]
    lines = [l for l in lines if l]
    W = H = None
    blocks, terms = [], []
    for p in lines:
        if p[0] == 'Outline:':
            W, H = int(p[1]), int(p[2])
        elif p[0] in ('NumBlocks:', 'NumTerminals:'):
            continue
        elif len(p) == 4 and p[1] == 'terminal':
            terms.append((p[0], float(p[2]), float(p[3])))
        elif len(p) == 3:
            blocks.append((p[0], int(p[1]), int(p[2])))
        else:
            raise ValueError(f'{path}: unexpected line {p}')
    return W, H, blocks, terms


def read_nets(path):
    with open(path) as f:
        tok = f.read().replace('\r', '').split()
    nets = []
    i = 0
    nnets = None
    while i < len(tok):
        if tok[i] == 'NumNets:':
            nnets = int(tok[i + 1])
            i += 2
        elif tok[i] == 'NetDegree:':
            d = int(tok[i + 1])
            nets.append(tok[i + 2:i + 2 + d])
            i += 2 + d
        else:
            raise ValueError(f'{path}: unexpected token {tok[i]}')
    if nnets is not None and nnets != len(nets):
        raise ValueError(f'{path}: NumNets {nnets} != {len(nets)}')
    return nets


def build(W, H, blocks, terms, nets, soft):
    c = lab3fmt.Case(W, H)
    bnames = {b[0] for b in blocks}
    for n, w, h in blocks:
        c.blocks.append(('SOFT', n, w * h) if soft else ('FIXED', n, w, h))
    tnames = {t[0] for t in terms}
    for n, x, y in terms:
        px, py = lab3fmt.project_to_border(x, y, W, H)
        c.terminals.append((n, px, py))
    is_block = bnames.__contains__
    pairs = []
    dropped = 0
    for net in nets:
        for p in net:
            if p not in bnames and p not in tnames:
                raise ValueError(f'unknown pin {p}')
        pr = lab3fmt.net_to_pairs(net, is_block)
        if not pr:
            dropped += 1
        pairs += pr
    c.conns = lab3fmt.dedupe_connections(pairs, is_block)
    return c, dropped


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    for name, variants in (('ami33', ('fixed',)), ('ami49', ('fixed', 'soft'))):
        W, H, blocks, terms = read_block(os.path.join(src, name, name + '.block'))
        nets = read_nets(os.path.join(src, name, name + '.nets'))
        moved = sum(1 for _, x, y in terms
                    if lab3fmt.project_to_border(x, y, W, H) != (x, y))
        for v in variants:
            c, dropped = build(W, H, blocks, terms, nets, soft=(v == 'soft'))
            err = lab3fmt.sanity(c)
            if err:
                sys.exit(f'{name}_{v}: ' + '; '.join(err[:10]))
            path = os.path.join(out, f'{name}_{v}.in')
            c.write(path)
            print(f'{path}: n={len(c.blocks)} terms={len(c.terminals)} '
                  f'(projected {moved}) conns={len(c.conns)} nets={len(nets)} '
                  f'(no-block nets dropped {dropped}) ws={c.whitespace():.4f}')


if __name__ == '__main__':
    main()
