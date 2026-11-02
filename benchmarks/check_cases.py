#!/usr/bin/env python3
"""Re-parse converted Lab3 inputs, run the section-3 sanity checks, and print a
Markdown table (used for benchmarks/README.md).

Usage: python3 -I benchmarks/check_cases.py <file.in> ...
Exit 0 if every file passes, 1 otherwise.
"""

import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lab3fmt  # noqa: E402


def main():
    bad = 0
    print('| file | n (soft/fixed/pre) | terminals | conns | MIBs | boundaries (corners) | outline | whitespace |')
    print('|---|---|---|---|---|---|---|---|')
    for path in sys.argv[1:]:
        try:
            c = lab3fmt.read_case(path)
            err = lab3fmt.sanity(c)
        except Exception as e:  # noqa: BLE001
            err = [f'parse error: {e}']
            c = None
        if err:
            bad += 1
            print(f'FAIL {path}: ' + '; '.join(err[:10]), file=sys.stderr)
            continue
        k = collections.Counter(b[0] for b in c.blocks)
        corners = sum(1 for _, e in c.boundaries if len(e) == 2)
        print(f'| {os.path.relpath(path)} | {len(c.blocks)} ({k["SOFT"]}/{k["FIXED"]}/{k["PREPLACED"]}) '
              f'| {len(c.terminals)} | {len(c.conns)} | {len(c.mibs)} '
              f'| {len(c.boundaries)} ({corners}) | {c.W}x{c.H} | {100 * c.whitespace():.2f}% |')
    print(f'{len(sys.argv) - 1 - bad} ok, {bad} failed', file=sys.stderr)
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
