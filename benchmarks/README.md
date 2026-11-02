# Tuning benchmarks (TA-internal)

Public floorplanning benchmarks converted to the Lab3 input format (docs/design.md §3).
They are used **only** to tune the TA reference solver (docs/design.md §9 step 1, ADR-0002).
Never release them to students.

Every file has `WLTHRESHOLD 0.0` (placeholder) and `ASPECTRATIO 3`. None of these cases has a
known construction that passes the gate with full boundary/MIB points; they are tuning inputs,
not generator cases.

## Sources and licenses

| Suite | Source | License / terms |
|---|---|---|
| MCNC ami33, ami49 | `reference/Lab2/ami33`, `reference/Lab2/ami49` (Lab2 `.block`/`.nets` format, originally MCNC'93 building-block benchmarks) | MCNC benchmarks, freely distributed for research and teaching; no formal license file |
| GSRC n100/n200/n300 | http://vlsicad.eecs.umich.edu/BK/GSRCbench/ (`HARD/` and `SOFT/`, UCSC/UCLA bookshelf `.blocks`/`.nets`/`.pl`, created 2000 by H. Wu) | GSRC Bookshelf, freely distributed for research; no formal license file |
| FloorSet-Lite validation | https://github.com/IntelLabs/FloorSet, data `LiteTensorDataTest.tar.gz` from https://huggingface.co/datasets/IntelLabs/FloorSet (100 samples, configs 21–120) | CC BY 4.0 (Intel Labs). Attribution: FloorSet, Intel Labs |

Downloads used (each archive extracted into its own new empty directory):

```
# GSRC
mkdir -p gsrc/HARD gsrc/SOFT
for d in HARD SOFT; do for n in n100 n200 n300; do for e in blocks nets pl; do
  curl -sSL -o gsrc/$d/$n.$e https://vlsicad.eecs.umich.edu/BK/GSRCbench/$d/$n.$e
done; done; done
# FloorSet-Lite validation (69 MB)
curl -sSL -o LiteTensorDataTest.tar.gz \
  https://huggingface.co/datasets/IntelLabs/FloorSet/resolve/main/LiteTensorDataTest.tar.gz
mkdir fsdata && tar xzf LiteTensorDataTest.tar.gz -C fsdata
```

## Regenerating

```
python3 -I benchmarks/convert_mcnc.py reference/Lab2 benchmarks/mcnc
python3 -I benchmarks/convert_gsrc.py <gsrc_dir> benchmarks/gsrc
uv run --no-project --with numpy python -I benchmarks/convert_floorset.py <fsdata>/LiteTensorDataTest benchmarks/floorset
python3 -I benchmarks/check_cases.py benchmarks/*/*.in      # sanity checks + table below
```

`lab3fmt.py` holds the shared writer, parser, projection, net decomposition and sanity checks.
MCNC/GSRC converters are stdlib-only. The FloorSet converter needs numpy but not torch: it reads
the `.pth` zip archives with a restricted unpickler that only rebuilds plain tensors.

## Conversion rules

Common to all suites:

- **Terminals** are moved onto the outline border: clamp into `[0,W]×[0,H]`, then snap to the
  nearest side, rounding to integers.
- **Multi-pin nets → 2-pin connections** (MCNC and GSRC): a net of degree `d ≤ 5` becomes a clique
  with weight 1 per pair; `d > 5` becomes a star from every pin to the net's first *block* pin,
  weight 1. Terminal–terminal pairs are skipped; a net without a block pin is dropped (none occur).
- **Dedupe**: identical unordered pairs have their weights summed, then the weight is rounded and
  clamped to `[1, 10]`. Self-pairs and terminal–terminal pairs are removed.
- **Scaling**: one integer linear factor `s = ceil(sqrt(10000 / min soft area))` per case, applied
  to every dimension and coordinate (areas × `s²`). Integer inputs therefore stay exact; no
  rounding of areas, fixed dimensions or preplaced positions happens (FloorSet's centering shift is an integer too).

### MCNC (`mcnc/`)

- `ami33_fixed`, `ami49_fixed`: every block FIXED with its original `w h`, original outline
  (`ami49_fixed` is the t₄₉ difficulty yardstick of §9 step 2). `s = 1`.
- `ami49_soft`: every block SOFT with `area = w·h`, same outline (min area 63504, so `s = 1`).
- Original terminal coordinates lie on a pad frame larger than the Lab2 outline (24 of 40 ami33 and
  12 of 22 ami49 terminals are outside); projection puts them on the outline border, so several
  terminals end up on corners.
- No MIBs, no boundaries.

### GSRC (`gsrc/`)

- `HARD` blocks → FIXED (`w·s`, `h·s`); `SOFT` blocks → SOFT with `area·s²`. GSRC's own aspect
  bounds are ignored; Lab3's `R = 3` applies. Hard and soft versions of one `n` share the same `s`
  (from the soft minimum area): n100 `s = 6`, n200/n300 `s = 9`.
- Outline: square with side `ceil(sqrt(A / (1 − ws)))` for `ws` = 15% and 10% (actual whitespace in
  the table).
- Pads: the `.pl` pad frame is `[0,800]²` for all cases, and every pad lies on its border; pads are
  scaled to `W/800`, `H/800` and then projected (projection only fixes rounding).
- No MIBs, no boundaries. n300 exceeds the 200-block bound of the graded cases on purpose.

### FloorSet-Lite (`floorset/`)

Data layout, from FloorSet `lite_dataset_test.py` and `iccad2026contest/iccad2026_evaluate.py`:
`litedata_1.pth = [[blocks(n,6), b2b(e,3), p2b(e,3), pins(r,2)]]`, block columns
`[area, fixed, preplaced, mib_id, cluster_id, boundary]`, boundary bitmask `1=left, 2=right,
4=top, 8=bottom`; `litelabel_1.pth = [[metrics(8), fp_sol(n,5,2)]]` (closed rectangle polygons).

- Selected samples: configs 21, 28, …, 119 (every 7 blocks, 15 files `fsNNN.in`, NNN = blocks).
  All 100 validation samples convert and pass the sanity checks with the same script.
- Kinds: `preplaced ≠ 0` → PREPLACED with `fp_sol` rectangle; else `fixed ≠ 0` → FIXED with `fp_sol`
  dimensions; else SOFT with the area target. All raw data are integers; `s` per case 10–25.
- Outline: the `fp_sol` bounding box (its min corner is (0,0) in every sample; this is checked),
  enlarged with the same aspect ratio to 10% whitespace: `W = ceil(s·bw·f)`, `H = ceil(s·bh·f)`,
  `f = sqrt(A / (0.9·s²·bw·bh))`. The scaled bbox is **centered** in the outline with integer
  offsets `dx = (W − s·bw)//2`, `dy = (H − s·bh)//2`, and every preplaced position is shifted by
  `(dx, dy)`. Preplaced blocks must lie strictly inside the outline (docs/design.md §3); if one
  still touches an edge after centering, `W` and/or `H` is increased by 2 (re-centering) until none
  does. (For the current 15 samples centering alone suffices, so outlines are unchanged.)
  (`fp_sol` itself has 2–4% whitespace.)
- Terminals: every pin `p<k>`, scaled by `s` and shifted by `(dx, dy)`, then projected. FloorSet pins sit on a frame larger
  than the block bbox, so many are clamped to the right/top side of the new outline.
- Connections: raw weights are small floats that are integer multiples of a per-case minimum
  `w_min`; weight = `round(w / w_min)`, summed over duplicate pairs, clamped to `[1, 10]`. Ratios up
  to 59 occur, so heavy nets are clipped at 10. b2b → block–block CONN, p2b → block–terminal CONN.
- MIB: per group, keep only SOFT members with the group's most common area target; drop groups
  with fewer than 2 members. (Every sample has one MIB group; fixed/preplaced members are removed.)
- Boundary: bitmask → edges; dropped on PREPLACED blocks and for opposite-edge codes (none occur).
  FloorSet boundaries refer to the solution's bounding box, while Lab3's refer to the outline; after
  enlargement `fp_sol` no longer touches RIGHT/TOP, so full boundary points are not guaranteed to be
  achievable.
- Grouping (cluster_id) is ignored.

## Files

All 30 files pass `check_cases.py` and are accepted by the verifier's input parser
(`lab3::parse_input_file` from `verifier/lab3_io.h`, compiled into a throwaway harness; the
`verifier` binary itself was not built yet when the files were made).

| file | n (soft/fixed/pre) | terminals | conns | MIBs | boundaries (corners) | outline | whitespace |
|---|---|---|---|---|---|---|---|
| benchmarks/mcnc/ami33_fixed.in | 33 (0/33/0) | 40 | 124 | 0 | 0 (0) | 1205x1095 | 12.36% |
| benchmarks/mcnc/ami49_fixed.in | 49 (0/49/0) | 22 | 260 | 0 | 0 (0) | 5336x7673 | 13.43% |
| benchmarks/mcnc/ami49_soft.in | 49 (49/0/0) | 22 | 260 | 0 | 0 (0) | 5336x7673 | 13.43% |
| benchmarks/gsrc/n100_hard_ws10.in | 100 (0/100/0) | 334 | 949 | 0 | 0 (0) | 2680x2680 | 10.03% |
| benchmarks/gsrc/n100_hard_ws15.in | 100 (0/100/0) | 334 | 949 | 0 | 0 (0) | 2758x2758 | 15.05% |
| benchmarks/gsrc/n100_soft_ws10.in | 100 (100/0/0) | 334 | 949 | 0 | 0 (0) | 2680x2680 | 10.03% |
| benchmarks/gsrc/n100_soft_ws15.in | 100 (100/0/0) | 334 | 949 | 0 | 0 (0) | 2758x2758 | 15.05% |
| benchmarks/gsrc/n200_hard_ws10.in | 200 (0/200/0) | 564 | 2172 | 0 | 0 (0) | 3977x3977 | 10.02% |
| benchmarks/gsrc/n200_hard_ws15.in | 200 (0/200/0) | 564 | 2172 | 0 | 0 (0) | 4092x4092 | 15.01% |
| benchmarks/gsrc/n200_soft_ws10.in | 200 (200/0/0) | 564 | 2172 | 0 | 0 (0) | 3977x3977 | 10.02% |
| benchmarks/gsrc/n200_soft_ws15.in | 200 (200/0/0) | 564 | 2172 | 0 | 0 (0) | 4092x4092 | 15.01% |
| benchmarks/gsrc/n300_hard_ws10.in | 300 (0/300/0) | 569 | 2664 | 0 | 0 (0) | 4959x4959 | 10.02% |
| benchmarks/gsrc/n300_hard_ws15.in | 300 (0/300/0) | 569 | 2664 | 0 | 0 (0) | 5103x5103 | 15.03% |
| benchmarks/gsrc/n300_soft_ws10.in | 300 (300/0/0) | 569 | 2664 | 0 | 0 (0) | 4959x4959 | 10.02% |
| benchmarks/gsrc/n300_soft_ws15.in | 300 (300/0/0) | 569 | 2664 | 0 | 0 (0) | 5103x5103 | 15.03% |
| benchmarks/floorset/fs021.in | 21 (18/2/1) | 68 | 126 | 1 | 10 (4) | 1104x671 | 10.07% |
| benchmarks/floorset/fs028.in | 28 (25/2/1) | 104 | 254 | 1 | 17 (4) | 2061x1981 | 10.04% |
| benchmarks/floorset/fs035.in | 35 (31/3/1) | 135 | 293 | 1 | 18 (3) | 2448x1944 | 10.06% |
| benchmarks/floorset/fs042.in | 42 (36/3/3) | 119 | 486 | 1 | 17 (4) | 1472x1185 | 10.11% |
| benchmarks/floorset/fs049.in | 49 (41/7/1) | 143 | 638 | 1 | 18 (4) | 2651x2295 | 10.05% |
| benchmarks/floorset/fs056.in | 56 (46/7/3) | 170 | 372 | 1 | 19 (4) | 2390x3118 | 10.03% |
| benchmarks/floorset/fs063.in | 63 (56/6/1) | 183 | 457 | 1 | 22 (4) | 2846x3305 | 10.02% |
| benchmarks/floorset/fs070.in | 70 (59/7/4) | 197 | 783 | 1 | 23 (3) | 2297x2680 | 10.03% |
| benchmarks/floorset/fs077.in | 77 (69/6/2) | 167 | 511 | 1 | 26 (4) | 2794x3300 | 10.03% |
| benchmarks/floorset/fs084.in | 84 (68/12/4) | 172 | 969 | 1 | 21 (4) | 2245x3193 | 10.00% |
| benchmarks/floorset/fs091.in | 91 (71/13/7) | 185 | 1972 | 1 | 25 (3) | 3114x4486 | 10.02% |
| benchmarks/floorset/fs098.in | 98 (88/7/3) | 328 | 3658 | 1 | 28 (4) | 3637x5262 | 10.02% |
| benchmarks/floorset/fs105.in | 105 (90/14/1) | 295 | 1713 | 1 | 31 (4) | 3012x5784 | 10.02% |
| benchmarks/floorset/fs112.in | 112 (101/10/1) | 270 | 4832 | 1 | 32 (3) | 3169x5956 | 10.01% |
| benchmarks/floorset/fs119.in | 119 (104/8/7) | 424 | 2471 | 1 | 31 (4) | 2938x5255 | 10.03% |
