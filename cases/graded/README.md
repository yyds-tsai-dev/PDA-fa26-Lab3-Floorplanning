# Graded cases (calibrated 2026-11-02 on the course workstation)

Do not release the hidden `.in` files or any `.construction.out`. Students get only
`public_ws0.1100.in`, shipped under the name `public.in`, together with `spec/Lab3.pdf` and the
workstation-built verifier.

Calibration: design.md section 9 (search from 15% whitespace down in 0.5% steps; keep the
lowest level where the frozen reference, `solver/fp.vda04` md5 `23f5eaf0d47702ed7674515f502ee0e4`,
is feasible with full boundary and MIB on 10/10 seeds with time-to-feasible <= 60 s at `--time 110`;
T = median wirelength of those 10 runs, rounded up to 0.5; the `WLTHRESHOLD` line holds T).
Runs were pinned one per physical core (`taskset -c 0-9` / `10-19`, hyper-thread siblings idle),
`python3.11 tools/calibrate.py --pregen`. Reports, logs and the student-proxy results are in
`calibration/` here (summarised in docs/design.md section 9.7).

| file | generator seed | whitespace | T (WLTHRESHOLD) | reference ttf median / max (s) | reference runs over 1.05T | T / construction WL | md5 |
|---|---|---|---|---|---|---|---|
| public_ws0.1100.in | 9001 | 11.0% | 286355.5 | 0.39 / 9.2 | 2 / 10 | 0.956 | 240932e9ed9b482dbede59de51831b90 |
| hidden1_ws0.0500.in | 9002 | 5.0% | 645484.0 | 3.5 / 35.1 | 1 / 10 | 1.075 | 1599486faf494f8d9c1ea650b707e9ce |
| hidden2_ws0.0500.in | 9003 | 5.0% | 879026.5 | 3.8 / 49.8 | 4 / 10 | 1.122 | 6ee8dd19d008714d26537db1132582b5 |
| hidden3_ws0.0800.in | 9004 (`--mib-groups 16 --mib-k-min 4 --mib-k-max 6`) | 8.0% | 1520007.0 | 8.3 / 42.4 | 2 / 10 | 1.141 | f265110f3c440867a19c95430b0998da |
| hidden4_ws0.0650.in | 9006 (`--mib-groups 20 --mib-k-min 5 --mib-k-max 6`) | 6.5% | 2022729.0 | 11.8 / 44.9 | 2 / 10 | 1.210 | 2f8f54dde092a454137e1a3c7d54ee41 |

Public whitespace = largest hidden value (8.0%) + 3 points = 11.0% (`--fixed-ws`, same T rule).

hidden4 was first generated with seed 9005. That instance passed only at the 15% start level
(10/10 at 15.0%, 5-6/10 at 14.5%, 10/10 at 14.0% and 13.0% but ttf up to 80 s, 7/10 at 12%,
5/10 at 10%; on NVL4 0/10 at 10% while seeds 9006-9010 are 10/10 there), so it was rejected under
the section 9.2 rule and the next seed, 9006, was used. The rejected files are kept in `rejected/`
(together with the first public calibration at 12%, `public_ws0.1200`, which the 15% hidden4 value implied)
and their reports in `calibration/hidden4_g9005/` and `calibration/public_r1/`.

Every `.construction.out` passes the verifier with full boundary and MIB on its `.in`.
