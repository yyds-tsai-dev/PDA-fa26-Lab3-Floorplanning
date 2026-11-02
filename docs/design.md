# Lab3 Design Spec (TA-internal)

Source of truth for every Lab3 deliverable. Vocabulary follows `CONTEXT.md`; decisions with history are in `docs/adr/`. The student-facing LaTeX spec is written from this file.

## 1. Problem

Fixed-outline floorplanning derived from ICCAD 2026 CAD Contest Problem C (Intel FloorSet). Place every block inside the outline without overlap. Soft blocks choose their own integer width/height. Fixed-shape blocks have given dimensions. Preplaced blocks have given dimensions and position. Rotation is not allowed. Optimize wirelength and satisfy MIB and boundary constraints for partial points.

## 2. Command line

```
make            # produces ./Lab3
./Lab3 <input> <output>
./verifier <input> <output>
```

Runtime limit 120 s wall-clock per case, measured by the TA on the course workstation. The program must be single-threaded.

## 3. Input format

Tokens are separated by whitespace. Sections appear in exactly this order, and each `<count>` line is followed by exactly that many records. There are no comments in real input files.

```
OUTLINE <W> <H>
WLTHRESHOLD <T>
ASPECTRATIO <R>
BLOCKS <n>
SOFT <name> <area>
FIXED <name> <w> <h>
PREPLACED <name> <w> <h> <x> <y>
...                                   (n block records, the three kinds in any order)
TERMINALS <r>
TERMINAL <name> <x> <y>
CONNECTIONS <m>
CONN <name1> <name2> <weight>
MIBS <q>
MIB <k> <name1> ... <namek>
BOUNDARIES <p>
BOUNDARY <name> <EDGE> [<EDGE>]
```

- All numbers are integers, except `T`, which is a non-negative multiple of 0.5 printed with one decimal (for example `123456.5`).
- `R` is an integer (3 in all released cases). A soft block must satisfy `R*h >= w` and `R*w >= h`.
- Names contain no whitespace, are unique across blocks and terminals, and are never `LEFT`, `RIGHT`, `TOP`, or `BOTTOM`.
- `area >= 10000` for every soft block.
- Bounds: every input coordinate and dimension is in `[0, 100000]` (32-bit safe); soft `area <= W·H <= 10^10`; `n <= 200`, terminals `r <= 100`, connections `m <= 1000`, MIB groups `q <= 20`, boundaries `p <= 100`. Products (`w·h`, `100·|w·h − area|`), total area and WL need 64-bit (WL can exceed 2^31). The verifier parses all integers as int64. Basis (`gen.py`): largest case (hidden4, `--whitespace 0.05`, seeds 1-3) has `W,H ≈ 6000–7500`, 446–459 connections, 72 terminals, 12 MIB groups, 40 boundaries; outline scales as `sqrt(n·480²)`, and the min-soft-area upscale factor is 1 in practice.
- `WLTHRESHOLD`: `T > 0` in all released cases.
- Token-based input: the verifier reads input as whitespace-separated tokens (line breaks not significant). Output is parsed line by line (blank lines skipped; one record per line).
- Whitespace ratio (hidden-case definition, on area targets): `1 − (Σ soft area + Σ fixed/preplaced w·h) / (W·H)`.
- Every preplaced block lies strictly inside the outline (it touches no outline edge), and no two preplaced blocks overlap.
- A terminal is a point on the outline border: `x ∈ {0, W}` with `0 <= y <= H`, or `y ∈ {0, H}` with `0 <= x <= W`.
- `CONN`: at least one endpoint is a block; never terminal–terminal; never a block to itself; each unordered pair appears at most once; `1 <= weight <= 10`.
- `MIB`: `k >= 2`; members are soft blocks with identical area targets; a block belongs to at most one MIB group.
- `BOUNDARY`: `EDGE ∈ {LEFT, RIGHT, TOP, BOTTOM}`. Two edges must be adjacent (a corner); `LEFT RIGHT` and `TOP BOTTOM` never occur. Only soft or fixed-shape blocks carry boundary constraints. A block appears in at most one `BOUNDARY` record.
- Generator guarantee: every case has a known solution that passes the gate and satisfies every MIB group and boundary constraint.

## 4. Output format

```
WIRELENGTH <wl>
BLOCKS <n>
<name> <x> <y> <w> <h>
...                                   (n lines, any order)
```

- `<wl>` is the decimal wirelength (recommended: one decimal place).
- `<x> <y>` is the lower-left corner; all four numbers are integers.
- Every input block appears exactly once. There are no unknown names and no extra tokens.

## 5. Wirelength

Block centroid: `cx = x + w/2`, `cy = y + h/2`. For each connection with weight `c` between endpoints `a` and `b` (a terminal's point is its own coordinates):

```
WL = Σ c · (|cx_a − cx_b| + |cy_a − cy_b|)
```

With integer coordinates, `WL` is a multiple of 0.5. The verifier computes `2·WL` exactly in 64-bit integers, using doubled coordinates.

## 6. Verifier

### Gate (10 points; any failure ⇒ case score 0)

Checked in this order. The verifier reports the first violation with specifics (for example `Overlap: b12 b37`).

1. Output parses: header tokens, block count equals `n`, each line has 5 tokens, integers parse, no trailing tokens.
2. Every input block name appears exactly once, and there are no unknown names.
3. `w > 0`, `h > 0`.
4. Fixed-shape: `(w,h)` equals the input. Preplaced: `(x,y,w,h)` equals the input.
5. Soft area tolerance: `100·|w·h − area| <= area` (64-bit).
6. Soft aspect ratio: `R·h >= w` and `R·w >= h`.
7. Inside outline: `x >= 0`, `y >= 0`, `x+w <= W`, `y+h <= H`.
8. Overlap-free: no pair has `max(x1,x2) < min(x1+w1,x2+w2)` and `max(y1,y2) < min(y1+h1,y2+h2)`. Touching is allowed.

The runtime limit and the single-thread rule are part of the gate but are enforced by the TA grading script, not the verifier.

### Scored constraints (only on a feasible output)

- **WL (5)**: if `|reported − WL| > 0.01`, the score is 0. Otherwise it is `5` if `WL <= T`, `5·(1.05T − WL)/(0.05T)` if `T < WL <= 1.05T`, and `0` if `WL > 1.05T`.
- **Boundary (3)**: `3 · satisfied / p`. An edge is touched when LEFT `x == 0`, RIGHT `x+w == W`, BOTTOM `y == 0`, TOP `y+h == H`. A block is satisfied only when all of its listed edges are touched.
- **MIB (2)**: `2 · satisfied / q`. A group is satisfied when all members have identical `(w,h)`.
- If `p == 0` or `q == 0`, that item gives full points. Released cases always have `p, q > 0`.

### Report

The verifier prints each gate check as PASS/FAIL, then `Reported WL`, `Computed WL`, `Threshold T`, the item scores with 2 decimals, and `Total x.xx / 20`. It exits 0 on a feasible output and 1 otherwise. The grading verifier and the released verifier are the same binary. It does not check runtime. The release build is stripped, with libstdc++/libgcc linked statically (see §10).

## 7. Test cases

One public case and four hidden cases. Every case contains every constraint type. Difficulty rises with block count, whitespace ratio, and constraint density.

| Case | blocks (soft / fixed / preplaced) | terminals | MIB groups | boundary (corners) |
|---|---|---|---|---|
| public | 50 (30/14/6) | 24 | 4 | 10 (2) |
| hidden1 | 80 (48/22/10) | 36 | 6 | 16 (3) |
| hidden2 | 120 (72/34/14) | 48 | 8 | 24 (4) |
| hidden3 | 160 (96/45/19) | 60 | 16 (MIB-heavy, k = 4–6) | 32 (4) |
| hidden4 | 200 (120/56/24) | 72 | 20 (MIB-heavy, k = 5–6) | 40 (4) |

- MIB-heavy (decided 2026-11-01, see docs/research/student-proxy.md "Gate hardening sweep"): hidden3 and hidden4 are generated with `--mib-groups 16 --mib-k-min 4 --mib-k-max 6` and `--mib-groups 20 --mib-k-min 5 --mib-k-max 6`; the other cases keep the preset defaults (group size 2–4). It is the one structure lever that is representation-fair (the sequence-pair reference and the coordinate-repair student proxy both pass the gate) and still costs an AI-written solver WL and boundary points. Realised sizes at gen seeds 9004/9006: hidden3 7×4, 1×5, 8×6 (81 of 96 soft blocks), hidden4 13×5, 7×6 (107 of 120); both constructions pass the verifier with full boundary and MIB. hidden4 was first generated with seed 9005; that instance was rejected under the §9.2 start-level rule (§9.7) and seed 9006 is used.
- Generator seeds of the graded cases: public 9001, hidden1 9002, hidden2 9003, hidden3 9004, hidden4 9006. Calibrated whitespace and T: §9.7 and cases/graded/README.md.
- At most 4 corner constraints exist per case (a corner block contains the outline corner point). Corners are counted inside the boundary total.
- Connections: about 3–4 per block on average, mostly between spatial neighbours in the generator's construction, plus about 10% long-range connections. Each terminal connects to 1–2 nearby blocks. Weights are 1–10.
- Preplaced blocks are spread through the interior, not only along the border.
- Whitespace is calibrated (§9): each hidden case at the tightest value the reference solver handles reliably. Public is the loosest case.
- The spec announces hidden ranges: at most 200 blocks, whitespace at least the calibrated minimum, and every constraint type present.

## 8. Generator

Python 3.10, stdlib plus numpy only. Deterministic for a given seed.

1. Choose the outline and recursively dissect it into `n` cells with zero whitespace. Keep cell aspect ratios within about 2.5, so that shrunk blocks satisfy R = 3. MIB groups come from splitting a region into k equal strips, which gives identical cells.
2. Shrink each cell to a block, so that total whitespace hits the target. Anchor each block inside its cell; boundary-constrained blocks are anchored to their required outline edges or corner. Identical cells in an MIB group are shrunk identically.
3. Assign kinds. Preplaced blocks are chosen only from interior cells, never touch an outline edge, and keep their constructed position. Fixed-shape blocks keep their constructed `(w,h)`. Soft blocks get `area = w·h` of the constructed block. Scale coordinates so that every soft area is at least 10000.
4. Place terminals on the border near boundary-adjacent blocks. Generate connections as in §7.
5. Emit `<case>.in` and the construction as `<case>.construction.out`. The construction must pass the verifier with full boundary and MIB points. `WLTHRESHOLD` is a placeholder until calibration rewrites it.

## 9. Calibration

1. **Reference solver frozen first.** It is tuned only on public benchmarks (MCNC ami33/ami49, GSRC, FloorSet-Lite converted) and on a generated dev set with different seeds. It is never tuned on the five graded cases (ADR-0002). **Freeze record (2026-11-02 03:39 Asia/Taipei):** `solver/fp.cpp` md5 `2950f613c0b80d54fd9b5680b2d8537f` (Phase G = Phase F behaviour; Phase H measured compiler flags and second-based knobs on the workstation, all a wash at 10 seeds, docs/research/adhd-solver-round3.md). The calibrated binary is the workstation's default build, `make` in `solver/` on vda04 with g++ 11.2.0 and `-O2 -Wall -Wextra -std=c++17` (md5 `23f5eaf0d47702ed7674515f502ee0e4`, kept in the repo as `solver/fp.vda04`); the NVL4 build of the same source is `solver/fp` (md5 `b3121dfba8a8fdf939bcd23739f97023`). The default build was chosen over the optional `make fp-ref` (`-O2 -march=native -ffp-contract=off`, +4% moves/s, WL inside seed noise) so that anyone running `make` on the workstation reproduces the calibrated binary byte for byte.
2. **Hidden whitespace (tight)**: hidden cases are meant to be hard to make feasible, so that a quickly AI-written solver fails the gate. For each hidden case, lower whitespace in 0.5% steps and keep the lowest value at which the frozen strong reference solver, under the student runtime budget, reaches gate + full boundary + full MIB on all 10 seeds with time-to-feasible ≤ 60 s (half the limit, as margin for host load and seeds outside the sample; measured on the workstation itself, item 7). Each hidden case gets its own tightest value; larger cases may end up slightly looser, since block count already adds difficulty. The search starts at 15% whitespace and steps down to 3% over a pre-generated grid. An instance that passes only at the 15% start level has no calibrated whitespace: it is rejected, the next generator seed is used, and the rejected instance and its numbers are recorded (first applied to hidden4 seed 9005, §9.7). Public is exempt and looser: its whitespace is the largest hidden value plus 3 percentage points, rounded up to the 0.5% grid and capped at 12%, evaluated with `--fixed-ws` (10 seeds, same T rule).
3. (Dropped) The ami49 t₄₉ yardstick: ami49 without rotation is feasible in about 0.08 s even in the weak configuration (E2), so it gives no useful lower bound.
4. **T**: the median WL over 10 seeds under the student runtime limit (120 s wall, single-threaded; the solver's internal budget leaves margin for I/O so every run itself passes the limit), rounded up to a multiple of 0.5. Then `WLTHRESHOLD := T`.
5. **Sanity**: a feasibility-only run (WL weight 0) must score WL > 1.05T, i.e. 0 WL points.
6. **Student proxy**: Codex writes a solver from the student spec alone. Its scores show the expected student distribution and expose ambiguities in the spec. If it beats the reference WL, the reference solver is not strong enough: improve it, refreeze, and recalibrate.
7. **Calibration host = the course workstation** (140.113.202.33, decided 2026-11-01 once the TA account existed). Measured 2026-11-01: Intel Xeon E5-2630 v4 @ 2.20 GHz (fixed clock), Rocky Linux 8.10, g++ 11.2.0, 20 physical cores with hyper-threading (cpu N and N+20 are siblings) shared by about 19 users. An uncontended physical core runs the reference at 340–360k final-anneal moves/s on hidden2_ws0.0500 (NVL4, same case, same binary flags, 2 seeds: 486k–518k), so one workstation core is about 0.7× an NVL4 core; a job whose hyper-thread sibling is busy drops further to about 60% of that (195–216k), and the first probe (157k moves/s, ttf 16 s, WL 1.117× the NVL4 median) was taken with a busy sibling. Both effects make a T measured on NVL4 unfair; the second is contention, which the TA controls at grading time (§11) and students experience when the host is busy. Calibration runs one job per physical core with the sibling idle (`--cores 0-19`, `--jobs ≤ 20`, or `10-19` with 10 jobs when another batch runs), at low host load, and records loadavg per level. The workstation has no numpy on any python (3.6.8 / 3.9 / 3.11), so cases are pre-generated locally for the whole whitespace grid and shipped; `tools/calibrate.py --pregen DIR` reads them instead of calling the generator, and the tool is run with `python3.11` (python3 is 3.6.8). Earlier NVL4 numbers (G0 preview: hidden1 3.5%, hidden2 5.0%, hidden3 4.5%, hidden4 6.5%) are superseded by the workstation calibration below.

   **Calibration result (2026-11-02, workstation):**

   | case | gen seed | whitespace | T | ref ttf median / max (s) | ref runs over 1.05T | T / construction WL | loadavg at level |
   |---|---|---|---|---|---|---|---|
   | public | 9001 | 11.0% | 286355.5 | 0.39 / 9.2 | 2 of 10 | 0.956 | 6.65 |
   | hidden1 | 9002 | 5.0% | 645484.0 | 3.5 / 35.1 | 1 of 10 | 1.075 | 20.1 |
   | hidden2 | 9003 | 5.0% | 879026.5 | 3.8 / 49.8 | 4 of 10 | 1.122 | 20.1 |
   | hidden3 | 9004 | 8.0% | 1520007.0 | 8.3 / 42.4 | 2 of 10 | 1.141 | 10.0 |
   | hidden4 | 9006 | 6.5% | 2022729.0 | 11.8 / 44.9 | 2 of 10 | 1.210 | 14.4 |

   The calibration ran on 2026-11-02 from 03:40 to 05:42 Asia/Taipei on the workstation (vda04) with the frozen default-build binary (md5 `23f5eaf0d47702ed7674515f502ee0e4`), using `python3.11 tools/calibrate.py --pregen` with `--time 110 --ttf-max 60`, 10 seeds per level, and the search starting at 15% whitespace and stepping down. Each job ran on its own physical core with the hyper-thread sibling idle (`--cores 0-9` and `10-19`). The hidden1 and hidden2 levels were measured while two 10-job batches ran, i.e. 20 jobs on 20 cores, which explains their loadavg of 20.1; during the hidden4 run other users were present beside our 10 jobs. hidden3 was generated with `--mib-groups 16 --mib-k-min 4 --mib-k-max 6` and hidden4 with `--mib-groups 20 --mib-k-min 5 --mib-k-max 6`. At the next hidden4 level, 6.0%, the reference passed on 9 of 10 seeds: all 10 were feasible, but one took 83 s to first feasibility. The sanity run (`--set s_wl0=1`) exceeded 1.05T on all five cases. Public follows the item 2 rule: the largest hidden whitespace (8.0%) plus 3 percentage points gives 11.0%, evaluated with `--fixed-ws`. The graded files and their md5 sums are in cases/graded/ (see its README.md).

   hidden4 was first generated with seed 9005, and that instance was rejected under the item 2 start-level rule. On the workstation it passed at 15.0% on 10 of 10 seeds twice (seeds 1-10, ttf max 32.8 s; seeds 11-20, ttf max 15.7 s), but only 6 of 10 and 5 of 10 at 14.5%; at 14.0% and 13.0% all 10 seeds were feasible but the ttf maximum was 74.6 s and 80.0 s; at 12% it passed on 7 of 10 and at 10% on 5 of 10. On NVL4 (about 1.4x faster per core) at 10%, seed 9005 passed on 0 of 10 while seeds 9006 to 9010 each passed on 10 of 10 (ttf max 29-50 s); at 12%, seed 9005 passed on 4 of 10 and the others on 10 of 10. The round-6 student proxy (coordinate repair) passes the gate on the 9005 instance at 10% on 5 of 5 seeds (full boundary on 3 of 5, WL 2.2-3.0x the construction). The 9005 instance therefore exposes a weakness of the sequence-pair reference, not a generally infeasible instance. Its cause was not identified: the block, MIB and preplaced statistics look like those of the other seeds. The rejected files are kept in scratch only.

   Reading: on hidden1 to hidden3, T is 7.5-14% above the generator's planted construction, so at 5-8% whitespace on the workstation the reference does not recover the planted packing. On hidden2, four of the ten reference seeds themselves score 0 WL points (max/median 1.107). T is therefore beatable, but the reference's own seed spread is wide. The user kept T as the median, with 0 points above 1.05T.
8. **Load discipline.** The solver is deadline-scheduled, so host load changes its search path, not only its wall time: in the Phase F A/Bs (40 jobs on a host at load 85–130) bench wall reached 125.9 s with the solver's own clock at 107.8 s, from launch latency before `main()`; a clean re-run gave 108.8 s. Calibration and the T measurement therefore run with `--jobs ≤ 20` on dedicated cores (`tools/calibrate.py --cores`), and the time-to-feasible test uses the solver's own FIRSTFEAS time. On the workstation the cores are not dedicated (other users are not pinned), so the tool records `/proc/loadavg` before every whitespace level and calibration runs off-peak.

## 10. Release package

- `Lab3.pdf`: built from LaTeX (xelatex), in English, compact (at most 6 pages, no bibliography). It includes a one-page toy example (input, output, figure, verifier report) and names ICCAD 2026 Problem C / FloorSet as the origin.
- `verifier`: stripped; libstdc++/libgcc static, glibc dynamic (glibc ≥ 2.25). Fully static if glibc-static becomes available.
- `public.in`.
- No reference solution of any case.

## 11. Administration

- Submission: a `Student_ID/` folder with source (`.cpp`/`.h`), `Makefile`, and `readme.txt` (how to build and run), packed with `tar cvf Student_ID.tar Student_ID` and uploaded to e3. Any other format scores 0. Failure to compile or run on the workstation scores 0.
- Release 2026-11-10 00:00. Deadline 2026-11-23 23:59 (Asia/Taipei). Late submissions lose 20% per day (as in reference/p2_floorplanning.pdf); any part of a day counts as a full day (after 11/23 23:59 up to 11/24 23:59 is −20%, up to 11/25 23:59 is −40%, and so on). The package goes to the head TA and the professor for review on 2026-11-04.
- Plagiarism: code similarity is checked, and both parties score 0 on confirmation.
- Grading runs: on the course workstation, one submission per physical core with its hyper-thread sibling idle (`taskset -c <k>` with k in 0–19 and cpus 20–39 unused), at most 20 concurrent, at low host load, so every run gets the same compute the calibration had (§9.7). A run that exceeds 120 s is re-run once alone before it is failed.
- No report. No stated C++ standard or library restriction. AI tools are allowed.

## 12. Repository layout

```
spec/        LaTeX source of Lab3.pdf, figures
verifier/    verifier.cpp, Makefile (static build)
generator/   gen.py
solver/      reference solver (C++, TA-internal)
cases/       public.in, hidden*.in, dev/ (tuning set), *.construction.out
benchmarks/  converters + converted MCNC/GSRC/FloorSet cases
tools/       grade.sh (timeout + single-thread check), plot.py (SVG), calibrate.py
docs/        design.md, adr/, research/, papers
```
