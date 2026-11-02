#!/usr/bin/env python3
"""Calibrate whitespace and WLTHRESHOLD for one case (docs/design.md section 9; Python 3.7+, stdlib only).

    python3 tools/calibrate.py --solver solver/fp --preset hidden1 --seed 5001 --ws-start 0.12 \\
        --out-dir cases/graded --cores 100-139 --jobs 10

    # course workstation (no numpy, so no generator): use cases pre-generated elsewhere with gen.py
    python3.11 tools/calibrate.py --solver solver/fp --verifier verifier/verifier --preset hidden1 --seed 9002 \\
        --ws-start 0.12 --ws-min 0.03 --pregen pregen --out-dir cal/hidden1 --cores 20-29 --jobs 10

--pregen DIR replaces the generator call by a lookup of DIR/<preset>_ws<WS>.in and
DIR/<preset>_ws<WS>.construction.out (WS with 4 decimals, e.g. hidden2_ws0.0500); both are copied into
--out-dir and only the copies are modified.  --seed is then recorded only (it must match the pregen seed).

1. Whitespace search: generate the case at --ws-start, run --runs seeds (1..runs) pinned one per core
   (`--time T --mode strong`), verify every output.  PASS iff all runs are verifier-feasible with full boundary
   and full MIB and solver ttf <= --ttf-max.  On PASS lower ws by --ws-step and repeat; stop at the first FAIL or
   at --ws-min and keep the lowest passing ws.  --fixed-ws evaluates only --ws-start (public case).
2. T = median verifier WL of the runs at the chosen ws, rounded UP to a multiple of 0.5 (exact integer
   arithmetic on 2*WL).  The WLTHRESHOLD line of the chosen .in is rewritten and the verifier re-run.
3. Sanity (section 9.5): one run with the WL weight removed must give WL > 1.05T.  fp has no WL-weight knob,
   so the check is skipped unless --sanity-args supplies the extra solver arguments for such a run.
4. Writes DIR/calibrate_log.json (every ws tried) and DIR/calibrate_report.json, prints the report.

Unlike bench.py this tool may write under cases/graded/ (it creates the graded cases).  Exit 0 on success,
1 if --ws-start itself fails, 2 on usage errors.
"""
import argparse
import concurrent.futures as cf
import hashlib
import json
import math
import os
import queue
import re
import shutil
import statistics
import subprocess
import sys
import time

MAX_JOBS = 40
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
RESULT_RE = re.compile(r"^RESULT (.*)$", re.M)


def parse_cores(spec):
    out = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def parse_result(stderr):
    m = None
    for m in RESULT_RE.finditer(stderr):
        pass
    if not m:
        return {}
    return dict(kv.split("=", 1) for kv in m.group(1).split() if "=" in kv)


def run_verifier(verifier, inp, outp):
    try:
        p = subprocess.run([verifier, inp, outp], capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return {"v_feasible": 0, "v_note": "verifier timeout"}
    txt = p.stdout
    r = {"v_feasible": 1 if "Result: FEASIBLE" in txt else 0}
    m = re.search(r"Computed WL\s+(\S+)", txt)
    r["v_wl"] = float(m.group(1)) if m else None
    m = re.search(r"Boundary\s+\S+ / 3\s+\((\d+)/(\d+) satisfied\)", txt)
    r["bnd"] = [int(m.group(1)), int(m.group(2))] if m else None
    m = re.search(r"MIB\s+\S+ / 2\s+\((\d+)/(\d+) satisfied\)", txt)
    r["mib"] = [int(m.group(1)), int(m.group(2))] if m else None
    m = re.search(r"WL score\s+(\S+) / 5", txt)
    r["v_wl_score"] = float(m.group(1)) if m else None
    m = re.search(r"Total\s+(\S+) / 20", txt)
    r["v_total"] = float(m.group(1)) if m else 0.0
    if not r["v_feasible"]:
        m = re.search(r"Result: INFEASIBLE \((.*)\)", txt)
        r["v_note"] = m.group(1) if m else "infeasible"
    return r


def pct(vals, q):
    v = sorted(vals)
    k = (len(v) - 1) * q
    lo = int(k)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (k - lo)


class Runner:
    def __init__(self, a, cores):
        self.a = a
        self.solver = os.path.abspath(a.solver)
        self.verifier = os.path.abspath(a.verifier)
        self.q = queue.Queue()
        for c in cores[:a.jobs]:
            self.q.put(c)

    def one(self, inp, outp, seed, extra=()):
        core = self.q.get()
        try:
            cmd = ["taskset", "-c", str(core), self.solver, inp, outp, "--time", str(self.a.time),
                   "--seed", str(seed), "--mode", "strong"] + list(extra)
            t0 = time.monotonic()
            try:
                p = subprocess.run(cmd, capture_output=True, text=True, timeout=self.a.time + 60)
                err, rc = p.stderr, p.returncode
            except subprocess.TimeoutExpired as e:
                err = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or "")
                rc = -9
            wall = time.monotonic() - t0
        finally:
            self.q.put(core)
        with open(outp[:-4] + ".log", "w") as f:
            f.write(err)
        res = parse_result(err)
        row = {"seed": seed, "rc": rc, "wall": round(wall, 2), "out": outp,
               "feasible": int(res.get("feasible", 0)), "ttf": float(res.get("ttf", -1))}
        if os.path.exists(outp):
            row.update(run_verifier(self.verifier, inp, outp))
        else:
            row.update({"v_feasible": 0, "v_note": "no output"})
        return row

    def batch(self, inp, outdir, seeds):
        os.makedirs(outdir, exist_ok=True)
        with cf.ThreadPoolExecutor(max_workers=self.a.jobs) as ex:
            futs = [ex.submit(self.one, inp, os.path.join(outdir, f"s{s}.out"), s) for s in seeds]
            rows = [f.result() for f in futs]
        rows.sort(key=lambda r: r["seed"])
        return rows


def full(pair):
    return bool(pair) and pair[0] == pair[1]


def run_ok(r, ttf_max):
    return (r.get("v_feasible") == 1 and r.get("feasible") == 1 and full(r.get("bnd")) and full(r.get("mib"))
            and 0 <= r["ttf"] <= ttf_max)


def generate(a, ws, base):
    cmd = [sys.executable, os.path.join(ROOT, "generator", "gen.py"), "--preset", a.preset,
           "--whitespace", f"{ws:g}", "--seed", str(a.seed), "--out", base]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0 or not os.path.exists(base + ".in"):
        raise RuntimeError("generator failed: " + p.stderr + p.stdout)


def use_pregen(a, ws, base):
    src = os.path.join(a.pregen, f"{a.preset}_ws{ws:.4f}")
    for ext in (".in", ".construction.out"):
        if not os.path.exists(src + ext):
            raise RuntimeError("pregen file missing: " + src + ext)
        shutil.copyfile(src + ext, base + ext)


def loadavg():
    try:
        with open("/proc/loadavg") as f:
            return [float(x) for x in f.read().split()[:3]]
    except (OSError, ValueError):
        return None


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solver", required=True)
    ap.add_argument("--preset", required=True, choices=["public", "hidden1", "hidden2", "hidden3", "hidden4"])
    ap.add_argument("--seed", type=int, required=True, help="generator seed")
    ap.add_argument("--ws-start", type=float, required=True)
    ap.add_argument("--ws-step", type=float, default=0.005)
    ap.add_argument("--ws-min", type=float, default=0.03)
    ap.add_argument("--runs", type=int, default=10)
    ap.add_argument("--time", type=float, default=110.0)
    ap.add_argument("--ttf-max", type=float, default=60.0)
    ap.add_argument("--cores", default=None, help="core list, e.g. 100-109 or 3,5,7 (default: highest --jobs cores)")
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--fixed-ws", action="store_true")
    ap.add_argument("--sanity-args", default="", help="extra solver args for the WL-weight-0 run (none exist in fp)")
    ap.add_argument("--verifier", default=os.path.join(ROOT, "verifier", "verifier"))
    ap.add_argument("--pregen", default=None, help="dir of pre-generated <preset>_ws<WS>.in/.construction.out")
    a = ap.parse_args()

    if a.jobs < 1 or a.jobs > MAX_JOBS:
        print(f"--jobs must be in [1, {MAX_JOBS}] (shared host)", file=sys.stderr)
        return 2
    if a.runs < 1 or a.ws_step <= 0:
        print("bad --runs / --ws-step", file=sys.stderr)
        return 2
    cores = parse_cores(a.cores) if a.cores else sorted(os.sched_getaffinity(0))[-a.jobs:][::-1]
    if len(cores) < a.jobs:
        print("not enough cores for --jobs", file=sys.stderr)
        return 2
    if a.pregen:
        a.pregen = os.path.abspath(a.pregen)
        if not os.path.isdir(a.pregen):
            print("--pregen dir not found: " + a.pregen, file=sys.stderr)
            return 2
    solver_md5 = md5_file(a.solver)
    out_dir = os.path.abspath(a.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    runner = Runner(a, cores)
    seeds = list(range(1, a.runs + 1))
    print(f"calibrate: {a.preset} gseed {a.seed} ws {a.ws_start} step {a.ws_step} min {a.ws_min} runs {a.runs} "
          f"time {a.time}s ttf-max {a.ttf_max}s jobs {a.jobs} cores {cores[:a.jobs]} out {out_dir}", flush=True)

    log, best = [], None
    ws_i = round(a.ws_start * 1e6)          # integer micro-units: no float drift across steps
    step_i = round(a.ws_step * 1e6)
    min_i = round(a.ws_min * 1e6)
    while True:
        ws = ws_i / 1e6
        tag = f"{a.preset}_ws{ws:.4f}"
        base = os.path.join(out_dir, tag)
        if a.pregen:
            use_pregen(a, ws, base)
        else:
            generate(a, ws, base)
        load = loadavg()
        rows = runner.batch(base + ".in", os.path.join(out_dir, tag + "_runs"), seeds)
        ok = [run_ok(r, a.ttf_max) for r in rows]
        ttfs = [r["ttf"] for r in rows if r["ttf"] >= 0]
        ent = {"ws": ws, "pass": all(ok), "runs": len(rows), "ok_runs": sum(ok),
               "feasible": sum(r.get("v_feasible", 0) for r in rows),
               "full_boundary": sum(full(r.get("bnd")) for r in rows),
               "full_mib": sum(full(r.get("mib")) for r in rows),
               "ttf_ok": sum(1 for r in rows if 0 <= r["ttf"] <= a.ttf_max),
               "ttf_max": max(ttfs) if ttfs else None, "base": base,
               "loadavg": load, "cpu_count": os.cpu_count()}
        log.append(ent)
        print(f"  ws {ws:.4f}: {'PASS' if ent['pass'] else 'FAIL'} ok {ent['ok_runs']}/{ent['runs']} "
              f"feas {ent['feasible']} bnd {ent['full_boundary']} mib {ent['full_mib']} "
              f"ttf_max {ent['ttf_max']} load {load}", flush=True)
        with open(os.path.join(out_dir, "calibrate_log.json"), "w") as f:
            json.dump(log, f, indent=1)
        if not ent["pass"]:
            break
        best = (ws, base, rows)
        if a.fixed_ws or ws_i - step_i < min_i:
            break
        ws_i -= step_i

    if best is None:
        print(f"FAIL: start ws {a.ws_start} does not pass", file=sys.stderr)
        rep = {"preset": a.preset, "gen_seed": a.seed, "ok": False, "reason": "ws-start failed",
               "solver": os.path.abspath(a.solver), "solver_md5": solver_md5, "pregen": a.pregen, "log": log}
        with open(os.path.join(out_dir, "calibrate_report.json"), "w") as f:
            json.dump(rep, f, indent=1)
        return 1

    ws, base, rows = best
    inp = base + ".in"
    # T from the verifier's WL, in exact integers of 2*WL.
    d = sorted(round(r["v_wl"] * 2) for r in rows)
    n = len(d)
    m4 = d[n // 2] * 2 if n % 2 else d[n // 2 - 1] + d[n // 2]     # 4 * median(WL)
    t2 = -(-m4 // 2)                                                # ceil(2 * median) = T / 0.5
    T = t2 / 2
    with open(inp) as f:
        text = f.read()
    text, cnt = re.subn(r"^WLTHRESHOLD.*$", f"WLTHRESHOLD {T:.1f}", text, count=1, flags=re.M)
    if cnt != 1:
        print("no WLTHRESHOLD line in " + inp, file=sys.stderr)
        return 2
    with open(inp, "w") as f:
        f.write(text)
    for r in rows:
        r.update({k: v for k, v in run_verifier(runner.verifier, inp, r["out"]).items() if k != "v_wl"})
    wls = [r["v_wl"] for r in rows]
    med = statistics.median(wls)

    # Sanity: WL weight 0.
    sanity = {"status": "skipped", "note": "fp has no WL-weight Config key (--set); pass --sanity-args to enable"}
    if a.sanity_args.strip():
        sr = runner.one(inp, os.path.join(out_dir, f"{a.preset}_ws{ws:.4f}_sanity.out"), 1, a.sanity_args.split())
        passed = bool(sr.get("v_feasible")) and sr.get("v_wl") is not None and sr["v_wl"] > 1.05 * T
        sanity = {"status": "pass" if passed else "fail", "wl": sr.get("v_wl"), "threshold_1.05T": 1.05 * T,
                  "args": a.sanity_args}

    rep = {"preset": a.preset, "gen_seed": a.seed, "ok": True, "ws": ws, "T": T, "input": inp,
           "wl": {"min": min(wls), "median": med, "p90": pct(wls, 0.9), "max": max(wls),
                  "max_over_median": max(wls) / med if med else None},
           "runs": len(rows), "full_wl_runs": sum(1 for w in wls if round(w * 2) <= t2),
           "over_1.05T_runs": sum(1 for w in wls if w > 1.05 * T),
           "wl_score_mean": statistics.mean(r["v_wl_score"] for r in rows if r.get("v_wl_score") is not None),
           "total_mean": statistics.mean(r["v_total"] for r in rows), "total_min": min(r["v_total"] for r in rows),
           "ttf": {"median": statistics.median(r["ttf"] for r in rows), "max": max(r["ttf"] for r in rows)},
           "sanity": sanity, "solver": os.path.abspath(a.solver), "solver_md5": solver_md5, "pregen": a.pregen,
           "per_run": [{k: r.get(k) for k in ("seed", "ttf", "v_wl", "v_wl_score", "v_total", "bnd", "mib")}
                       for r in rows],
           "log": log}
    with open(os.path.join(out_dir, "calibrate_report.json"), "w") as f:
        json.dump(rep, f, indent=1)
    print()
    print(f"case {a.preset} (gen seed {a.seed}): ws {ws:.4f}  T {T:.1f}  input {inp}")
    print(f"  WL min {rep['wl']['min']:.1f}  median {med:.1f}  p90 {rep['wl']['p90']:.1f}  max {rep['wl']['max']:.1f}"
          f"  max/median {rep['wl']['max_over_median']:.4f}")
    print(f"  full WL (WL<=T): {rep['full_wl_runs']}/{len(rows)}   >1.05T: {rep['over_1.05T_runs']}/{len(rows)}"
          f"   total mean {rep['total_mean']:.2f} min {rep['total_min']:.2f}")
    print(f"  ttf median {rep['ttf']['median']:.2f}  max {rep['ttf']['max']:.2f}")
    print(f"  sanity: {sanity['status']}" + (f" (WL {sanity['wl']}, 1.05T {1.05 * T:.1f})" if 'wl' in sanity else
                                             f" - {sanity['note']}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
