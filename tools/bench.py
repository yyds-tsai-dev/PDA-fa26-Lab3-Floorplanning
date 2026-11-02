#!/usr/bin/env python3
"""Benchmark the TA reference solver (Python 3.10, stdlib only).

    python3 tools/bench.py --solver solver/fp --cases cases/dev/*.in \\
        --seeds 1-10 --time 24 --mode strong --jobs 20 --csv out.csv

Runs every (case, seed) job in parallel, each pinned to its own core with
`taskset -c`, then runs the verifier on each output and records:
case, seed, mode, solver feasible/ttf/wl/boundary/mib/moves (from the solver's
RESULT line), verifier feasible/WL/boundary/MIB/total, WL ratio against
<case>.construction.out when present, wall time, and a disagreement flag
(solver claims feasible but the verifier does not, or the reverse).
Prints a per-case summary at the end.

Refuses any case under cases/graded/ (ADR-0002: the reference solver is never
tuned on the graded cases).
"""
import argparse
import concurrent.futures as cf
import csv
import hashlib
import os
import queue
import re
import statistics
import subprocess
import sys
import tempfile
import time

MAX_JOBS = 40


def parse_seeds(spec):
    seeds = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            seeds.extend(range(int(a), int(b) + 1))
        else:
            seeds.append(int(part))
    return seeds


def is_graded(path):
    parts = os.path.realpath(path).split(os.sep)
    for i in range(len(parts) - 1):
        if parts[i] == "cases" and parts[i + 1] == "graded":
            return True
    return False


def default_cores(jobs):
    # Upper range of the host's CPUs that this process may use (shared host).
    avail = sorted(os.sched_getaffinity(0))
    return avail[-jobs:][::-1]


RESULT_RE = re.compile(r"^RESULT (.*)$", re.M)


def parse_result(stderr):
    m = None
    for m in RESULT_RE.finditer(stderr):
        pass
    if not m:
        return {}
    out = {}
    for kv in m.group(1).split():
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def run_verifier(verifier, inp, outp):
    try:
        p = subprocess.run([verifier, inp, outp], capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return {"v_feasible": 0, "v_note": "verifier timeout"}
    txt = p.stdout
    r = {"v_feasible": 1 if "Result: FEASIBLE" in txt else 0}
    m = re.search(r"Computed WL\s+(\S+)", txt)
    r["v_wl"] = float(m.group(1)) if m else ""
    m = re.search(r"Boundary\s+\S+ / 3\s+\((\d+)/(\d+) satisfied\)", txt)
    r["v_boundary"] = f"{m.group(1)}/{m.group(2)}" if m else ""
    m = re.search(r"MIB\s+\S+ / 2\s+\((\d+)/(\d+) satisfied\)", txt)
    r["v_mib"] = f"{m.group(1)}/{m.group(2)}" if m else ""
    m = re.search(r"Total\s+(\S+) / 20", txt)
    r["v_total"] = float(m.group(1)) if m else 0.0
    if "WL mismatch" in txt:
        r["v_note"] = "WL mismatch"
    if not r["v_feasible"]:
        m = re.search(r"Result: INFEASIBLE \((.*)\)", txt)
        r["v_note"] = m.group(1) if m else "infeasible"
    return r


def construction_wl(verifier, inp):
    base = inp[:-3] if inp.endswith(".in") else inp
    cons = base + ".construction.out"
    if not os.path.exists(cons):
        return None
    r = run_verifier(verifier, inp, cons)
    return r.get("v_wl") if r.get("v_feasible") else None


def file_info(path):
    with open(path, "rb") as f:
        sha = hashlib.sha256(f.read()).hexdigest()[:12]
    mt = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(os.path.getmtime(path)))
    return sha, mt


def pct(vals, q):
    if not vals:
        return None
    v = sorted(vals)
    k = (len(v) - 1) * q
    lo = int(k)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (k - lo)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solver", required=True)
    ap.add_argument("--cases", nargs="+", required=True)
    ap.add_argument("--seeds", default="1-10")
    ap.add_argument("--time", type=float, default=24.0)
    ap.add_argument("--mode", choices=["weak", "strong"], default="strong")
    ap.add_argument("--jobs", type=int, default=20)
    ap.add_argument("--csv", default=None)
    ap.add_argument("--verifier", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier", "verifier"))
    ap.add_argument("--outdir", default=None, help="where solver outputs/logs go (default: a new temp dir)")
    ap.add_argument("--cores", default=None, help="comma list of cores to pin to (default: highest --jobs cores)")
    ap.add_argument("--extra", default="", help="extra solver arguments, one string")
    a = ap.parse_args()

    bad = [c for c in a.cases if is_graded(c)]
    if bad:
        print("REFUSED: graded cases must never be used for solver tuning (ADR-0002): " + " ".join(bad), file=sys.stderr)
        return 2
    if a.jobs < 1 or a.jobs > MAX_JOBS:
        print(f"--jobs must be in [1, {MAX_JOBS}] (shared host)", file=sys.stderr)
        return 2
    for c in a.cases:
        if not os.path.isfile(c):
            print(f"no such case: {c}", file=sys.stderr)
            return 2
    solver = os.path.abspath(a.solver)
    verifier = os.path.abspath(a.verifier)
    cores = [int(x) for x in a.cores.split(",")] if a.cores else default_cores(a.jobs)
    if len(cores) < a.jobs:
        print("not enough cores for --jobs", file=sys.stderr)
        return 2
    outdir = a.outdir or tempfile.mkdtemp(prefix="lab3_bench.")
    os.makedirs(outdir, exist_ok=True)
    seeds = parse_seeds(a.seeds)
    print(f"bench: {len(a.cases)} cases x {len(seeds)} seeds, time {a.time}s, mode {a.mode}, jobs {a.jobs}, "
          f"cores {cores[:a.jobs]}, outdir {outdir}", flush=True)

    info = {c: file_info(c) for c in a.cases}
    cons = {c: construction_wl(verifier, c) for c in a.cases}

    pool_cores = queue.Queue()
    for c in cores[:a.jobs]:
        pool_cores.put(c)

    def job(case, seed):
        core = pool_cores.get()
        try:
            name = os.path.basename(case)[:-3] if case.endswith(".in") else os.path.basename(case)
            outp = os.path.join(outdir, f"{name}.{a.mode}.s{seed}.out")
            logp = outp[:-4] + ".log"
            cmd = ["taskset", "-c", str(core), solver, case, outp, "--time", str(a.time), "--seed", str(seed),
                   "--mode", a.mode] + (a.extra.split() if a.extra else [])
            t0 = time.monotonic()
            try:
                p = subprocess.run(cmd, capture_output=True, text=True, timeout=a.time + 60)
                err, rc = p.stderr, p.returncode
            except subprocess.TimeoutExpired as e:
                err, rc = (e.stderr or b"").decode() if isinstance(e.stderr, bytes) else (e.stderr or ""), -9
            wall = time.monotonic() - t0
        finally:
            pool_cores.put(core)
        with open(logp, "w") as f:
            f.write(err)
        res = parse_result(err)
        row = {"case": name, "seed": seed, "mode": a.mode, "rc": rc, "wall": round(wall, 2),
               "feasible": int(res.get("feasible", 0)), "ttf": float(res.get("ttf", -1)),
               "wl": res.get("wl", ""), "boundary": res.get("boundary", ""), "mib": res.get("mib", ""),
               "moves": res.get("moves", ""), "sha": info[case][0], "mtime": info[case][1], "out": outp}
        if os.path.exists(outp):
            row.update(run_verifier(verifier, case, outp))
        else:
            row.update({"v_feasible": 0, "v_note": "no output"})
        cw = cons[case]
        row["construction_wl"] = cw if cw is not None else ""
        row["wl_ratio"] = round(row["v_wl"] / cw, 4) if (cw and row.get("v_feasible") and row.get("v_wl")) else ""
        row["disagree"] = int(row["feasible"] != row["v_feasible"])
        return case, row

    rows = []
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(job, c, s) for c in a.cases for s in seeds]
        for f in cf.as_completed(futs):
            case, row = f.result()
            rows.append(row)
            print(f"  {row['case']} s{row['seed']}: feasible={row['feasible']}/{row['v_feasible']} "
                  f"ttf={row['ttf']:.2f} wl={row.get('v_wl', '')} ratio={row['wl_ratio']} "
                  f"total={row.get('v_total', '')} {row.get('v_note', '')}", flush=True)

    rows.sort(key=lambda r: (r["case"], r["seed"]))
    fields = ["case", "seed", "mode", "feasible", "v_feasible", "disagree", "ttf", "wl", "v_wl", "construction_wl",
              "wl_ratio", "boundary", "v_boundary", "mib", "v_mib", "v_total", "moves", "wall", "rc", "v_note",
              "sha", "mtime", "out"]
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in fields})

    print()
    print(f"{'case':<22}{'runs':>5}{'feas':>7}{'ttf_med':>9}{'ttf_p90':>9}{'ratio_med':>10}{'total_mean':>11}"
          f"{'Mmoves_med':>11}{'disagree':>9}")
    by = {}
    for r in rows:
        by.setdefault(r["case"], []).append(r)
    tot_dis = 0
    for case in sorted(by):
        rs = by[case]
        feas = [r for r in rs if r["v_feasible"]]
        ttfs = [r["ttf"] for r in rs if r["ttf"] >= 0]
        ratios = [r["wl_ratio"] for r in feas if r["wl_ratio"] != ""]
        tots = [r.get("v_total", 0.0) or 0.0 for r in rs]
        mv = [int(r["moves"]) for r in rs if r["moves"]]
        dis = sum(r["disagree"] for r in rs)
        tot_dis += dis
        f = lambda v, fmt: (fmt % v) if v is not None else "-"
        print(f"{case:<22}{len(rs):>5}{len(feas) / len(rs):>7.2f}{f(pct(ttfs, 0.5), '%.2f'):>9}"
              f"{f(pct(ttfs, 0.9), '%.2f'):>9}{f(statistics.median(ratios) if ratios else None, '%.4f'):>10}"
              f"{statistics.mean(tots):>11.2f}{f(statistics.median(mv) / 1e6 if mv else None, '%.2f'):>11}{dis:>9}")
    print(f"total disagreements (solver vs verifier feasibility): {tot_dis}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
