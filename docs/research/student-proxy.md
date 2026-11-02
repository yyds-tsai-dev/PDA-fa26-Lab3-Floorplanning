# Student proxy results (design.md §9.6)

Blind solvers written from spec/Lab3.pdf + the public case + the verifier binary only, run on our cases at 110 s, one seed, verifier-checked. WL is relative to the generator's construction.

| proxy | time spent | dev_tight (16) gate | WL/constr median (dev_tight) | dev hidden3/4 gate | notes |
|---|---|---|---|---|---|
| Codex v1 (2026-10-29 19:10) | 5 min | 10/16 | – (1.4–2.3 on feasible) | hidden4 0/4 | coordinate SA, no legalizer |
| Codex v2, half-finished (20:21–21:00, quota hit) | ~35 min | 15/16 | 2.4 | hidden3 2/2, hidden4 0/2 | SP SA + overlap-repair with GLS weights |
| Opus round 1 (22:13–23:27) | 75 min | **16/16**, full boundary/MIB | **1.51** (1.17–1.85) | hidden2/3/4 4/4 | repair starts from the analytical placement with a decaying WL weight; stall restarts; soft-boundary fallback after 80% of budget |
| Opus round 2 (23:30–01:10, "round 4" in its REPORT.md) | +100 min (≈3 h total) | **16/16**, full boundary/MIB | **1.36** (1.19–1.93) | hidden2/3/4 4/4 | multi-start in tight mode: fresh global placement + legalisation per attempt, attempt ends 12 s after its last improvement, best legal kept (3–7 attempts per run); its own A/B: −2.9% tuning, −2.2% hold-out |
| Opus round 3 (01:15–02:20, "round 5" in its REPORT.md) | +65 min (≈4 h total) | **16/16**, full boundary/MIB | **1.40** (1.18–1.89) | hidden2/3/4 4/4 | early abort of attempts >0.5% behind the best; its own A/B −1.5% tuning, hold-out neutral. Our single-seed median moved 1.36 → 1.40, inside its ±5% seed noise: no measurable gain this round |
| Opus round 4 (02:20–03:20, "round 6" in its REPORT.md) | +60 min (≈5 h total) | not re-evaluated on our cases | – | – | iterated local search: every second attempt permutes the 8% of blocks around a random block in the best legal placement and re-legalises; its own A/B −1.8% tuning (one lucky case; 0.997 without it), hold-out 0.996. Plateau: two rounds neutral |
| Opus round 5 (03:15–04:00, "round 7" in its REPORT.md) | +45 min (≈6 h total) | unchanged binary (round 6) | – | – | density-based global placement and "safe first" legalisation tried and rejected (1.013 / 1.011 vs round 6); three rounds without gain |
| Opus round 6 (04:00–04:45, "round 8" in its REPORT.md) | +45 min (≈7 h total) | unchanged binary (round 6) | – | – | order search after legality with repair-based re-legalisation: ~1 s per move at 5% whitespace, 1.006 / 1.014 vs round 6, rejected; legal-only swaps fit 0.65% of the time. Four rounds without gain |
| Opus round 7 (04:41–05:15, "round 9" in its REPORT.md) | +35 min (≈7.5 h total) | unchanged binary (round 6) | – | – | windowed re-pack legaliser (3k–47k moves/s) neutral: 1.006 / 0.995 vs round 6. The proxy declares a plateau: local reordering gains ≤0.5% while seed spread is 2–6%; a higher ceiling needs an analytical placer + constraint-graph/LP legalisation ("several days"). Rounds paused here pending the user's decision |

Reference solver at the same date: Phase E dev_tight × 10 seeds 144/160 feasible, median 1.276; Phase F step 1 (GLS-first) × 5 seeds 80/80, median 1.187; Phase F step 2 standalone homotopy 1.10–1.13 but 74–79/80 feasible (not adopted); **Phase F final (hybrid GLS-first → homotopy, 2026-10-30 02:40) × 10 seeds 159/160 feasible, median 1.061**, ttf median 4.1 s. Student round 3 is 1.32× Phase F.

Reading: after ~2 h of AI-assisted work the student proxy already passes the gate on every tight case (including hidden4_ws0.06_s1101, which the reference still fails on some seeds) with WL about 1.27× the reference's. After ~3 h it is at 1.14× the reference's step-1 median (1.36 vs 1.187; single seed per case, its own seed spread is ±3–5%). With T = reference median and 0 points at 1.05T it would still score 15/20 on hidden cases, but the margin is shrinking by roughly 10% per round; the reference must keep improving (hybrid homotopy, Phase G) to keep the WL points discriminating. Worst student case is hidden4_ws0.06_s1101 (1.93). The 24 h run continues in 2 h rounds; later rounds are appended here.

## Student round 3 binary at the G0 calibrated whitespace (2026-10-30 03:20)

Cases: the G0 lowest-passing instances (hidden1 ws 3.5%, hidden2 5.0%, hidden3 4.5%, hidden4 6.5%; gen seed 1001), 110 s, seeds default/3/7, T = Phase F median of 10 seeds.

| case | ws | gate (3 seeds) | boundary / MIB | student WL median / T | WL points |
|---|---|---|---|---|---|
| hidden1 | 3.5% | 3/3 | full | 1.23 | 0 |
| hidden2 | 5.0% | 3/3 | full | 1.31 | 0 |
| hidden3 | 4.5% | 3/3 | full | 1.41 | 0 |
| hidden4 | 6.5% | 3/3 | full | 1.35 | 0 |

Reading: the tight whitespace does not stop this proxy (≈4 h of AI-assisted work): it passes the gate on every seed down to 3.5% whitespace and scores 15/20 per hidden case. The gate is discriminating only against weaker attempts (Codex v1/v2). The WL margin (1.23–1.41 × T) is what separates it from the reference; keeping that margin is Phase G's job.

## Gate hardening sweep (2026-11-01)

Question: which test-case *structures* (not instances) make the gate harder for the round-6 student binary while Phase F still passes. `generator/gen.py` gained optional, backward-compatible knobs (`--size-tail`, `--big-cells/--big-frac/--big-kind`, `--pre-arrange line/cross/ring/largest`, `--pre-fill`, `--outline-aspect`, `--mib-k-min/max`); defaults are byte-identical to the previous output. Screening: 28 structures at hidden4 counts x ws {0.05, 0.07} x 3 generator seeds x {Phase F, Phase E, student}, 110 s, verifier-checked (504 runs; scratch `gate/sweep.csv`, `gate/summary.md`).

- Full passes per ws (of 84): Phase F 79 at 0.07 / 59 at 0.05; Phase E 54 / 21; student 76 / 69. At 5% whitespace the student passes more often than the reference.
- Walls, rings, crosses and large or many preplaced obstacles (`--pre-arrange`, `--pre-fill`, `--preplaced 60`) do not hurt the student (3/3 everywhere) and hurt the sequence-pair reference (S12 ring + 60 preplaced at 0.05: F 1/3, student 3/3; S20 spread + 60 preplaced at 0.05: F 0/3, student 2/3). Coordinate repair handles obstacles natively; the SP search rejects preplaced-order moves.
- Size heterogeneity is the only lever that beats the student where Phase F (seed 1) passes: `--size-tail 1.5` at 0.07 gives F 3/3, student 1/3 full (its repair stalls at 34-48 s with 17-19% overlap left), Phase E 0/3; `--big-cells 8 --big-frac 0.3 --big-kind fixed` at 0.07 gives F 2/3, student 1/3, E 1/3. Fairness flag: Phase E fails every heavy-tail structure, so the lever mainly measures legaliser strength. On 10 seeds Phase F itself is only 7/10 at ws 0.09 on `--size-tail 1.5` (gate/cal.csv), so the strict section 9 rule gives it no calibrated whitespace.
- MIB-heavy structures (16 groups of 4-6 or 20 groups of 5-6) are representation-fair (Phase E 3/3 at 0.07) and cost the student 25-77 s to first legal (F 4.5-13 s) with WL 1.76x the construction (F 1.06-1.11x); gate still passes at 0.07.
- Section 9 emulation (gen seed 1001, 10 reference seeds, ttf <= 60 s, ws 0.09 down to 0.045 in 0.005 steps; student 5 seeds per level; scratch `gate/cal.csv`, `gate/cal_table.md`): `--size-tail 1.5` and `--big-cells 8 --big-frac 0.3 --big-kind fixed` have no calibrated ws (Phase F 7/10 and 9/10 already at 0.09); `--size-tail 1.0` calibrates only at the 0.09 start level (student 5/5). `--mib-groups 20 --mib-k-min 5 --mib-k-max 6` calibrates at ws 0.075 (F ttf 17.8 s median / 52.5 s max, WL 1.27x construction); the student there is 3/5 full (5/5 gate; two runs fell back to soft boundary with 8/40 and 25/40), ttf 46.4/50.7 s on the full runs, WL 1.71x. `--mib-groups 16 --mib-k-min 4 --mib-k-max 6` calibrates at 0.07; student 5/5 full, ttf 25.0/26.1 s, WL 1.81x (F 8.3/47.5 s, 1.19x). `--big-cells 12 --big-frac 0.5 --big-kind pre` calibrates at 0.08; student 5/5 full, ttf 18.8/71.1 s, WL 1.47x (F 7.0/20.3 s, 1.00x).
- Decision (2026-11-01, confirmed by the user the same evening): no second search. The gate stays as calibrated by section 9 and WL points carry the discrimination; MIB-heavy is used on two hidden cases, hidden3 (`--mib-groups 16 --mib-k-min 4 --mib-k-max 6`) and hidden4 (`--mib-groups 20 --mib-k-min 5 --mib-k-max 6`, gen seed 9006 after seed 9005 was rejected under the section 9.2 start-level rule), see design.md section 7. MIB-heavy (16-20 groups of 4-6) is the one structure lever that may be used (representation-fair, costs the student boundary points on 2/5 seeds and 1.7-1.8x WL at the calibrated ws; it changes the section 7 MIB column). Size tails, big fixed blocks, walls, rings and preplaced-heavy mixes are not used. The generator knobs stay (defaults byte-identical).

## Round-6 binary on the workstation-calibrated graded cases (2026-11-02)

- Setup: the round-6 student proxy was built with g++ 11.2 on the course workstation (md5 `2f1d106fa19c4ef248e935716d8ffd6e`) and run on the five graded cases (cases/graded/) with a 110 s budget on 5 seeds: the default plus `LAB3_SEED=3,7,11,13`. The proxy's third command-line argument is the time budget, not a seed, so seeds are set only through `LAB3_SEED`. Each job ran on its own physical core; the runs started after calibration, at host load 13.
- Totals: gate 25/25; full boundary 24/25 (the only miss is hidden3 seed 7 with 23/32); full MIB 25/25.

| case | ws | gate | full boundary | WL/T min / median / max | score |
|---|---|---|---|---|---|
| public | 11.0% | 5/5 | 5/5 | 1.01 / 1.16 / 1.23 | 15/20 on four seeds; 18.76 on seed 7 (WL 1.012T) |
| hidden1 | 5.0% | 5/5 | 5/5 | 1.09 / 1.22 / 1.30 | 15/20 |
| hidden2 | 5.0% | 5/5 | 5/5 | 1.21 / 1.21 / 1.45 | 15/20 |
| hidden3 | 8.0% | 5/5 | 4/5 | 1.28 / 1.32 / 1.40 | 15/20 on four seeds; 14.16 on seed 7 |
| hidden4 | 6.5% | 5/5 | 5/5 | 1.44 / 1.52 / 2.29 (2.29 is the default seed) | 15/20 |

- Reading: the proxy passes every gate down to 5% whitespace, and at 6.5% on the MIB-heavy hidden4, so the gate separates only weak attempts. Its WL margin is 1.2-1.5x T on the hidden cases (one public seed reached 1.01T), so the proxy scores about 15/20 per case (mean about 15.1/20 over the five cases) against the reference's 18.2-18.8 mean (about 18.5).
- Rejected hidden4 instance (gen seed 9005): the proxy passes the gate there at 10% whitespace on 5 of 5 seeds (full boundary on 3 of 5, WL 2.2-3.0x the construction), while the reference passes on 0 of 10 seeds at 10% on NVL4 and on 5 of 10 on the workstation. See design.md section 9 item 7.
