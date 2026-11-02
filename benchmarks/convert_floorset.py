#!/usr/bin/env python3
"""Convert FloorSet-Lite validation samples (Intel, CC BY 4.0) to Lab3 input.

Usage: python -I benchmarks/convert_floorset.py <LiteTensorDataTest_dir> <out_dir> [n ...]
  <LiteTensorDataTest_dir> holds config_<n>/litedata_1.pth and litelabel_1.pth
  (extracted from LiteTensorDataTest.tar.gz on huggingface.co/datasets/IntelLabs/FloorSet).
  Default sample sizes: 21, 28, ..., 119 (15 samples). Needs numpy, not torch.

The .pth files are read with a restricted unpickler that only rebuilds plain
tensors into numpy arrays (no torch import, no arbitrary globals).
Data layout (FloorSet lite_dataset_test.py / iccad2026_evaluate.py):
  litedata  = [[blocks(n,6), b2b(e,3), p2b(e,3), pins(r,2)]]
      blocks columns: [area, fixed, preplaced, mib_id, cluster_id, boundary_code]
      b2b rows: [block_i, block_j, weight]; p2b rows: [pin, block, weight]
      boundary_code bitmask: 1=left, 2=right, 4=top, 8=bottom
  litelabel = [[metrics(8), fp_sol(n,5,2)]]  fp_sol = closed rectangle polygons
See benchmarks/README.md for the conversion rules.
"""

import collections
import math
import os
import pickle
import sys
import zipfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lab3fmt  # noqa: E402

WHITESPACE = 0.10
DEFAULT_SIZES = list(range(21, 121, 7))

_DT = {'FloatStorage': np.float32, 'DoubleStorage': np.float64,
       'LongStorage': np.int64, 'IntStorage': np.int32, 'BoolStorage': np.bool_,
       'ShortStorage': np.int16, 'ByteStorage': np.uint8, 'CharStorage': np.int8,
       'HalfStorage': np.float16}


class _StorageType:
    def __init__(self, name):
        self.name = name


def _rebuild_tensor(storage, offset, size, stride, *unused):
    if len(size) == 0:
        return storage[offset].item()
    return np.lib.stride_tricks.as_strided(
        storage[offset:], shape=size,
        strides=[x * storage.itemsize for x in stride]).copy()


def load_pth(path):
    """Load a torch zip-format .pth containing only tensors/lists into numpy."""
    z = zipfile.ZipFile(path)
    pkl = [n for n in z.namelist() if n.endswith('/data.pkl')]
    if len(pkl) != 1:
        raise ValueError(f'{path}: not a torch zip archive')
    prefix = pkl[0][:-len('data.pkl')]
    if any(n.endswith('/byteorder') for n in z.namelist()):
        bo = z.read(prefix + 'byteorder').decode().strip()
        if bo != 'little':
            raise ValueError(f'{path}: byteorder {bo}')

    class Unpickler(pickle.Unpickler):
        def find_class(self, mod, name):
            if mod == 'torch._utils' and name == '_rebuild_tensor_v2':
                return _rebuild_tensor
            if mod == 'torch' and name in _DT:
                return _StorageType(name)
            if mod == 'collections' and name == 'OrderedDict':
                return collections.OrderedDict
            raise pickle.UnpicklingError(f'blocked global {mod}.{name}')

        def persistent_load(self, pid):
            tag, st, key, loc, numel = pid
            if tag != 'storage' or not isinstance(st, _StorageType):
                raise pickle.UnpicklingError(f'bad persistent id {pid!r}')
            arr = np.frombuffer(z.read(f'{prefix}data/{key}'), dtype=_DT[st.name])
            if arr.size != numel:
                raise ValueError(f'{path}: storage {key} size mismatch')
            return arr

    with z.open(pkl[0]) as f:
        return Unpickler(f).load()


EDGE_BITS = ((1, 'LEFT'), (2, 'RIGHT'), (4, 'TOP'), (8, 'BOTTOM'))


def convert(data_dir, n):
    d = load_pth(os.path.join(data_dir, f'config_{n}', 'litedata_1.pth'))
    lab = load_pth(os.path.join(data_dir, f'config_{n}', 'litelabel_1.pth'))
    blk, b2b, p2b, pins = d[0]
    sol = lab[0][1]
    nb = blk.shape[0]
    if nb != n or sol.shape[0] != nb:
        raise ValueError(f'config_{n}: block count {nb} / {sol.shape[0]}')
    area = blk[:, 0]
    fixed = blk[:, 1] != 0
    pre = blk[:, 2] != 0
    mib = blk[:, 3].astype(int)
    bcode = blk[:, 5].astype(int)
    x0, x1 = sol[:, :, 0].min(1), sol[:, :, 0].max(1)
    y0, y1 = sol[:, :, 1].min(1), sol[:, :, 1].max(1)
    for arr in (area, x0, x1, y0, y1):
        if not np.all(arr == np.round(arr)):
            raise ValueError(f'config_{n}: non-integer data')
    if x0.min() != 0 or y0.min() != 0:
        raise ValueError(f'config_{n}: fp_sol bbox min corner is not (0,0)')
    bw, bh = int(x1.max()), int(y1.max())
    soft = ~(fixed | pre)

    # one integer linear scale so every soft area target is >= 10000
    s = int(math.ceil(math.sqrt(10000.0 / area[soft].min())))
    names = [f'b{i}' for i in range(nb)]
    c_blocks = []
    for i in range(nb):
        w, h = int(x1[i] - x0[i]) * s, int(y1[i] - y0[i]) * s
        if pre[i]:
            c_blocks.append(('PREPLACED', names[i], w, h, int(x0[i]) * s, int(y0[i]) * s))
        elif fixed[i]:
            c_blocks.append(('FIXED', names[i], w, h))
        else:
            c_blocks.append(('SOFT', names[i], int(area[i]) * s * s))

    # outline: same-aspect enlargement of the fp_sol bbox to 10% whitespace,
    # with the bbox centered in it (integer offset); preplaced blocks must lie
    # strictly inside, so grow the outline by 2 until none touches an edge
    tmp = lab3fmt.Case(1, 1)
    tmp.blocks = c_blocks
    A = tmp.total_area()
    Ws, Hs = bw * s, bh * s
    f = math.sqrt(A / ((1.0 - WHITESPACE) * Ws * Hs))
    if f > 1.0:
        W, H = int(math.ceil(Ws * f)), int(math.ceil(Hs * f))
    else:
        W, H = Ws, Hs
    pres = [b for b in c_blocks if b[0] == 'PREPLACED']
    while True:
        dx, dy = (W - Ws) // 2, (H - Hs) // 2
        okx = all(b[4] + dx > 0 and b[4] + dx + b[2] < W for b in pres)
        oky = all(b[5] + dy > 0 and b[5] + dy + b[3] < H for b in pres)
        if okx and oky:
            break
        if not okx:
            W += 2
        if not oky:
            H += 2
    c_blocks = [b[:4] + (b[4] + dx, b[5] + dy) if b[0] == 'PREPLACED' else b
                for b in c_blocks]
    c = lab3fmt.Case(W, H)
    c.blocks = c_blocks

    # terminals: every pin, scaled and shifted by (dx, dy), then projected onto the outline border
    for k in range(pins.shape[0]):
        px, py = float(pins[k, 0]) * s + dx, float(pins[k, 1]) * s + dy
        c.terminals.append((f'p{k}',) + lab3fmt.project_to_border(px, py, W, H))

    # connections: weight = round(w / w_min) per case, summed over duplicate
    # pairs, clamped to [1, 10]
    wmin = float(min(b2b[:, 2].min(), p2b[:, 2].min()))
    pairs = []
    for i, j, w in b2b:
        pairs.append((names[int(i)], names[int(j)], round(float(w) / wmin)))
    for k, i, w in p2b:
        pairs.append((names[int(i)], f'p{int(k)}', round(float(w) / wmin)))
    bset = set(names)
    c.conns = lab3fmt.dedupe_connections(pairs, bset.__contains__)

    # MIB groups: soft members only, with the group's most common area target
    dropped_mib = 0
    for g in sorted(set(mib[mib > 0].tolist())):
        mem = [i for i in range(nb) if mib[i] == g]
        sm = [i for i in mem if soft[i]]
        if sm:
            common = collections.Counter(int(area[i]) for i in sm).most_common(1)[0][0]
            sm = [i for i in sm if int(area[i]) == common]
        dropped_mib += len(mem) - len(sm)
        if len(sm) >= 2:
            c.mibs.append([names[i] for i in sm])

    # boundaries: not on preplaced, no opposite-edge combinations
    dropped_bnd = 0
    for i in range(nb):
        code = int(bcode[i])
        if code == 0:
            continue
        edges = [e for bit, e in EDGE_BITS if code & bit]
        if pre[i] or (code & 3) == 3 or (code & 12) == 12 or len(edges) > 2:
            dropped_bnd += 1
            continue
        c.boundaries.append((names[i], edges))

    info = dict(scale=s, bbox=(Ws, Hs), dropped_mib=dropped_mib,
                dropped_bnd=dropped_bnd, wmax_ratio=float(
                    max(b2b[:, 2].max(), p2b[:, 2].max()) / wmin))
    return c, info


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    src, out = sys.argv[1], sys.argv[2]
    sizes = [int(v) for v in sys.argv[3:]] or DEFAULT_SIZES
    os.makedirs(out, exist_ok=True)
    for n in sizes:
        c, info = convert(src, n)
        err = lab3fmt.sanity(c)
        name = f'fs{n:03d}'
        if err:
            sys.exit(f'{name}: ' + '; '.join(err[:10]))
        path = os.path.join(out, name + '.in')
        c.write(path)
        kinds = collections.Counter(b[0] for b in c.blocks)
        print(f'{path}: n={len(c.blocks)} soft/fixed/pre={kinds["SOFT"]}/'
              f'{kinds["FIXED"]}/{kinds["PREPLACED"]} scale={info["scale"]} '
              f'outline={c.W}x{c.H} terms={len(c.terminals)} conns={len(c.conns)} '
              f'mibs={len(c.mibs)} (members dropped {info["dropped_mib"]}) '
              f'bnd={len(c.boundaries)} (dropped {info["dropped_bnd"]}) '
              f'wratio_max={info["wmax_ratio"]:.0f} ws={c.whitespace():.4f}')


if __name__ == '__main__':
    main()
