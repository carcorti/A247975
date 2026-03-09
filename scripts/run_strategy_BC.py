#!/usr/bin/env python3
"""
run_strategy_BC.py
==================
Orchestrator for the extended minimality-verification campaign
of Sun's Conjecture 4.1(i) / OEIS A247975.

Strategy B+C (recommended by the revision plan):
  B1 – 50 smallest exceptional cases  (a(m) just above 10^8)
       Fills the key gap: no minimality check yet for ANY exceptional case.
  B2 – 100 random non-exceptional cases  (a(m) < 10^8)
       Extends the random sample from 30 → 130 cases.
  C  – 10–15 cases with a(m) > 10^9, incl. m=14740 and m=16167
       Demonstrates the verifier works for the most extreme cases.

Requires:
  verify_min    (compiled from verify_min.c; must be in the same directory
                 or on PATH — see Makefile / compile instructions)
  fase6_m100000_b2e9.csv   (Phase-6 output; path configurable below)

Usage:
  python3 run_strategy_BC.py [--csv PATH] [--bin PATH] [--seed N]
                             [--timeout-b1 S] [--timeout-b2 S] [--timeout-c S]
                             [--out-dir PATH]

Output:
  <out_dir>/verification_results.csv   — one row per case, machine-readable
  <out_dir>/verification_log.txt       — human-readable log with timings
  <out_dir>/verification_summary.txt   — final tally for the paper
"""

import argparse
import csv
import os
import random
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# ── Default configuration ────────────────────────────────────────────────────
DEFAULT_CSV     = "fase6_m100000_b2e9.csv"
DEFAULT_BIN     = "./verify_min"
DEFAULT_SEED    = 42
DEFAULT_TIMEOUT_B1 = 120   # seconds; B1 cases have a_m ~ 10^8–1.2×10^8 → ~10s each
DEFAULT_TIMEOUT_B2 =  60   # seconds; B2 cases have a_m < 10^8 → <10s each
DEFAULT_TIMEOUT_C  = 1800  # seconds; C cases have a_m > 10^9 → up to 30 min each
DEFAULT_OUT_DIR = "."

# ── How many cases per group ─────────────────────────────────────────────────
N_B1 = 100    # smallest exceptional
N_B2 = 150   # random non-exceptional
N_C  = 40    # large cases (a_m > 10^9); mandatory: m=14740 and m=16167

# ── Thresholds ────────────────────────────────────────────────────────────────
EXCEPTIONAL_THRESHOLD = 100_000_000   # a(m) > 10^8
B2_MAX_AM             =  99_999_999   # a(m) < 10^8 for group B2
C_MIN_AM              = 1_000_000_000  # a(m) > 10^9 for group C
# Mandatory cases for group C (from Table 1 of the paper)
C_MANDATORY = {14740, 16167}

# ─────────────────────────────────────────────────────────────────────────────

def parse_csv(path: str) -> list[dict]:
    """Parse the Phase-6 CSV.  Returns list of dicts for resolved cases only."""
    records = []
    with open(path, newline="") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("m,"):
                continue
            parts = line.split(",")
            if len(parts) < 4:
                continue
            if parts[2] == "NOT_FOUND":
                continue   # skip unresolved cases
            try:
                records.append({
                    "m":    int(parts[0]),
                    "pm":   int(parts[1]),
                    "am":   int(parts[2]),
                    "pam":  int(parts[3]),
                })
            except ValueError:
                continue
    return records


def select_cases(records: list[dict], seed: int) -> dict[str, list[dict]]:
    """Select cases according to Strategy B+C."""
    rng = random.Random(seed)

    # Group B1: smallest N_B1 exceptional cases
    exceptional = [r for r in records if r["am"] > EXCEPTIONAL_THRESHOLD]
    exceptional.sort(key=lambda r: r["am"])
    b1 = exceptional[:N_B1]

    # Group B2: random non-exceptional cases (a_m < 10^8)
    b1_ms = {r["m"] for r in b1}
    non_exc = [r for r in records
               if r["am"] <= B2_MAX_AM and r["m"] not in b1_ms]
    b2 = rng.sample(non_exc, min(N_B2, len(non_exc)))

    # Group C: large cases (a_m > 10^9), mandatory first, then top by a_m
    large = [r for r in records if r["am"] >= C_MIN_AM]
    large.sort(key=lambda r: r["am"], reverse=True)

    mandatory = [r for r in large if r["m"] in C_MANDATORY]
    optional  = [r for r in large if r["m"] not in C_MANDATORY]
    # Fill up to N_C, mandatory cases always included
    c_list = mandatory + optional
    c_list = c_list[:N_C]   # take at most N_C

    return {"B1": b1, "B2": b2, "C": c_list}


def estimate_seconds(am: int) -> float:
    """Rough wall-clock estimate (seconds) on Ryzen 9 7940HS.
    Linear scaling: empirical ~7s per 10^8 primes (from calibration runs)."""
    return max(1.0, 7.0 * am / 1e8)


def run_verifier(binary: str, record: dict, timeout: int) -> dict:
    """Run verify_min for one record; return a result dict."""
    m, pm, am, pam = record["m"], record["pm"], record["am"], record["pam"]
    cmd = [binary, str(m), str(pm), str(am), str(pam)]

    t_start = time.monotonic()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        elapsed = time.monotonic() - t_start
        stdout = proc.stdout.strip()
        stderr = proc.stderr.strip()   # progress lines from verifier

        # Parse the one-line result from stdout.
        # verify_min output format: STATUS m first_n first_p elapsed
        #   parts[0] = STATUS
        #   parts[1] = m       (redundant, already known)
        #   parts[2] = first_n (the solution index found)
        #   parts[3] = first_p (the prime p_{first_n})
        #   parts[4] = elapsed (reported by the C binary)
        parts = stdout.split()
        status = parts[0] if parts else "ERROR"
        if status == "VERIFIED":
            return {"status": "VERIFIED", "m": m, "am": am,
                    "first_n": int(parts[2]) if len(parts) > 2 else -1,
                    "first_p": int(parts[3]) if len(parts) > 3 else -1,
                    "elapsed": elapsed, "stderr": stderr}
        elif status == "MISMATCH":
            return {"status": "MISMATCH", "m": m, "am": am,
                    "first_n": int(parts[2]) if len(parts) > 2 else -1,
                    "first_p": int(parts[3]) if len(parts) > 3 else -1,
                    "elapsed": elapsed, "stderr": stderr,
                    "note": f"EARLIER solution found at n={parts[2] if len(parts) > 2 else '?'}"}
        elif status == "NOT_FOUND":
            return {"status": "NOT_FOUND", "m": m, "am": am,
                    "first_n": -1, "first_p": -1,
                    "elapsed": elapsed, "stderr": stderr}
        else:
            return {"status": "ERROR", "m": m, "am": am,
                    "first_n": -1, "first_p": -1,
                    "elapsed": elapsed,
                    "note": f"Unexpected output: {stdout!r}",
                    "stderr": stderr}

    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t_start
        return {"status": "TIMEOUT", "m": m, "am": am,
                "first_n": -1, "first_p": -1,
                "elapsed": elapsed, "note": f"Timed out after {timeout}s"}
    except Exception as e:
        elapsed = time.monotonic() - t_start
        return {"status": "ERROR", "m": m, "am": am,
                "first_n": -1, "first_p": -1,
                "elapsed": elapsed, "note": str(e)}


def fmt_time(seconds: float) -> str:
    if seconds < 120:
        return f"{seconds:.1f}s"
    return f"{seconds/60:.1f}min"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--csv",        default=DEFAULT_CSV,        help="Path to Phase-6 CSV")
    ap.add_argument("--bin",        default=DEFAULT_BIN,        help="Path to verify_min binary")
    ap.add_argument("--seed",       type=int, default=DEFAULT_SEED, help="Random seed")
    ap.add_argument("--timeout-b1", type=int, default=DEFAULT_TIMEOUT_B1,
                    help="Per-case timeout (s) for group B1")
    ap.add_argument("--timeout-b2", type=int, default=DEFAULT_TIMEOUT_B2,
                    help="Per-case timeout (s) for group B2")
    ap.add_argument("--timeout-c",  type=int, default=DEFAULT_TIMEOUT_C,
                    help="Per-case timeout (s) for group C")
    ap.add_argument("--out-dir",    default=DEFAULT_OUT_DIR,    help="Output directory")
    ap.add_argument("--dry-run",    action="store_true",
                    help="Print selection only, do not run verifier")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_csv  = out_dir / f"verification_results_{ts}.csv"
    log_file     = out_dir / f"verification_log_{ts}.txt"
    summary_file = out_dir / f"verification_summary_{ts}.txt"

    # ── Check binary ─────────────────────────────────────────────────────────
    if not args.dry_run:
        if not Path(args.bin).exists():
            print(f"ERROR: verifier binary not found: {args.bin}")
            print("Compile with:  gcc -O2 -o verify_min verify_min.c -lm")
            sys.exit(1)

    # ── Parse CSV ────────────────────────────────────────────────────────────
    print(f"Reading {args.csv} …")
    records = parse_csv(args.csv)
    print(f"  {len(records):,} resolved cases loaded.")

    # ── Select cases ─────────────────────────────────────────────────────────
    groups = select_cases(records, args.seed)
    timeouts = {"B1": args.timeout_b1, "B2": args.timeout_b2, "C": args.timeout_c}

    print()
    total_cases = sum(len(v) for v in groups.values())
    total_est   = 0.0
    for grp, cases in groups.items():
        est = sum(estimate_seconds(r["am"]) for r in cases)
        total_est += est
        am_min = min(r["am"] for r in cases)
        am_max = max(r["am"] for r in cases)
        print(f"  Group {grp}: {len(cases)} cases  "
              f"a(m) in [{am_min:,}, {am_max:,}]  "
              f"est. {fmt_time(est)}")
    print(f"\n  Total: {total_cases} cases  est. {fmt_time(total_est)}\n")

    if args.dry_run:
        print("=== DRY RUN: case list only ===")
        for grp, cases in groups.items():
            print(f"\n-- Group {grp} --")
            for r in cases:
                print(f"  m={r['m']:7d}  a_m={r['am']:>14,}  p_m={r['pm']:>10,}  "
                      f"est {fmt_time(estimate_seconds(r['am']))}")
        return

    # ── Run verifier ─────────────────────────────────────────────────────────
    all_results = []
    campaign_start = time.monotonic()

    with (open(results_csv, "w", newline="") as rf,
          open(log_file,    "w")             as lf):

        writer = csv.DictWriter(rf,
            fieldnames=["group","m","pm","am","pam","status",
                        "first_n","first_p","elapsed_s","note"])
        writer.writeheader()

        def log(msg):
            print(msg)
            lf.write(msg + "\n")
            lf.flush()

        log(f"=== A247975 minimality verification  {datetime.now().isoformat()} ===")
        log(f"Binary:  {args.bin}")
        log(f"CSV:     {args.csv}")
        log(f"Seed:    {args.seed}")
        log(f"Groups:  B1={len(groups['B1'])}  B2={len(groups['B2'])}  C={len(groups['C'])}")
        log("")

        case_idx = 0
        for grp, cases in groups.items():
            timeout = timeouts[grp]
            log(f"{'='*60}")
            log(f"Group {grp}  ({len(cases)} cases, timeout {timeout}s each)")
            log(f"{'='*60}")

            grp_ok = grp_mismatch = grp_timeout = grp_error = 0

            for r in cases:
                case_idx += 1
                elapsed_campaign = time.monotonic() - campaign_start
                est = estimate_seconds(r["am"])
                log(f"\n[{case_idx}/{total_cases}] m={r['m']}  "
                    f"a_m={r['am']:,}  p_m={r['pm']}  "
                    f"(est {fmt_time(est)}, campaign {fmt_time(elapsed_campaign)})")

                result = run_verifier(args.bin, r, timeout)
                result["group"] = grp
                result["pm"]    = r["pm"]
                result["pam"]   = r["pam"]
                result.setdefault("note", "")
                all_results.append(result)

                writer.writerow({
                    "group":     grp,
                    "m":         r["m"],
                    "pm":        r["pm"],
                    "am":        r["am"],
                    "pam":       r["pam"],
                    "status":    result["status"],
                    "first_n":   result["first_n"],
                    "first_p":   result["first_p"],
                    "elapsed_s": f"{result['elapsed']:.2f}",
                    "note":      result.get("note",""),
                })
                rf.flush()

                status = result["status"]
                if status == "VERIFIED":
                    grp_ok += 1
                    log(f"  ✓ VERIFIED  n={result['first_n']:,}  "
                        f"p_n={result['first_p']:,}  "
                        f"({fmt_time(result['elapsed'])})")
                elif status == "MISMATCH":
                    grp_mismatch += 1
                    log(f"  ✗ MISMATCH  earlier solution at n={result['first_n']:,}  "
                        f"({fmt_time(result['elapsed'])})")
                    log(f"    *** ACTION REQUIRED: check m={r['m']} in original dataset ***")
                elif status == "TIMEOUT":
                    grp_timeout += 1
                    log(f"  ⏱ TIMEOUT  (>{timeout}s)  "
                        f"consider increasing --timeout-{grp.lower()}")
                else:
                    grp_error += 1
                    log(f"  ? {status}  {result.get('note','')}")

            log(f"\nGroup {grp} done: "
                f"VERIFIED={grp_ok}  MISMATCH={grp_mismatch}  "
                f"TIMEOUT={grp_timeout}  ERROR/NOT_FOUND={grp_error}")

        # ── Summary ──────────────────────────────────────────────────────────
        total_elapsed = time.monotonic() - campaign_start
        n_verified  = sum(1 for r in all_results if r["status"] == "VERIFIED")
        n_mismatch  = sum(1 for r in all_results if r["status"] == "MISMATCH")
        n_timeout   = sum(1 for r in all_results if r["status"] == "TIMEOUT")
        n_other     = len(all_results) - n_verified - n_mismatch - n_timeout

        log("")
        log("=" * 60)
        log("FINAL SUMMARY")
        log("=" * 60)
        log(f"Total cases attempted : {len(all_results)}")
        log(f"  VERIFIED            : {n_verified}")
        log(f"  MISMATCH            : {n_mismatch}   (*** check immediately ***)")
        log(f"  TIMEOUT             : {n_timeout}")
        log(f"  OTHER               : {n_other}")
        log(f"Total wall time       : {fmt_time(total_elapsed)}")
        log("")

        # Group breakdown
        for grp in ("B1", "B2", "C"):
            grp_res = [r for r in all_results if r["group"] == grp]
            ok  = sum(1 for r in grp_res if r["status"] == "VERIFIED")
            tot = len(grp_res)
            log(f"  Group {grp}: {ok}/{tot} verified")

    # ── Write summary for paper ───────────────────────────────────────────────
    with open(summary_file, "w") as sf:
        sf.write("Independent minimality verification — extended campaign\n")
        sf.write(f"Date: {datetime.now().strftime('%Y-%m-%d')}\n")
        sf.write(f"Binary: {args.bin}  (verify_min.c, independent sieve)\n\n")
        sf.write(f"Groups selected (seed={args.seed}):\n")
        for grp, cases in groups.items():
            am_vals = [r["am"] for r in cases]
            sf.write(f"  {grp}: {len(cases)} cases, "
                     f"a(m) in [{min(am_vals):,}, {max(am_vals):,}]\n")
        sf.write(f"\nResults:\n")
        sf.write(f"  VERIFIED   : {n_verified} / {len(all_results)}\n")
        sf.write(f"  MISMATCH   : {n_mismatch}\n")
        sf.write(f"  TIMEOUT    : {n_timeout}\n")
        sf.write(f"  OTHER      : {n_other}\n")
        sf.write(f"  Wall time  : {fmt_time(total_elapsed)}\n\n")
        sf.write("Verification details (VERIFIED cases only, first 20):\n")
        sf.write(f"{'m':>8}  {'a(m)':>14}  {'p_m':>10}  {'group':>5}  {'seconds':>8}\n")
        sf.write("-" * 58 + "\n")
        verified_sorted = sorted(
            [r for r in all_results if r["status"] == "VERIFIED"],
            key=lambda r: r["am"], reverse=True)
        for r in verified_sorted[:20]:
            sf.write(f"{r['m']:>8}  {r['am']:>14,}  {r['pm']:>10,}  "
                     f"{r['group']:>5}  {r['elapsed']:>8.1f}\n")
        if len(verified_sorted) > 20:
            sf.write(f"  ... and {len(verified_sorted)-20} more (see verification_results_*.csv)\n")

    print()
    print(f"Results CSV  : {results_csv}")
    print(f"Log file     : {log_file}")
    print(f"Summary      : {summary_file}")
    print(f"Wall time    : {fmt_time(total_elapsed)}")
    print()
    if n_mismatch > 0:
        print("*** WARNING: MISMATCHES FOUND — review log immediately ***")
    elif n_verified == len(all_results):
        print("All cases VERIFIED. ✓")
    else:
        print(f"{n_verified}/{len(all_results)} verified. "
              f"Check log for timeouts/errors.")


if __name__ == "__main__":
    main()
