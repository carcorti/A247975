#!/usr/bin/env python3
# =============================================================================
# run_verify_ext_batch.py
# Batch orchestrator for verify_min_ext — Groups E2 and E3
#
# Runs verify_min_ext on a stratified sample of 28 Extended cases:
#   E2: 15 cases with a_m in [2e9, 1e10)  — covers the low Extended range
#   E3:  9 cases with a_m in [1e10, 8.2e10] — covers the high Extended range
#
# (Group E1, the 4 top Cramér-ratio cases, is run separately.)
#
# Usage:
#   python3 run_verify_ext_batch.py [--jobs N] [--verifier PATH]
#
#   --jobs N       Max parallel processes (default: 4)
#   --verifier P   Path to verify_min_ext binary (default: ./verify_min_ext)
#
# Output:
#   verification_min_ext_results.csv   — one row per case (appended in real time)
#   run_verify_ext_batch.log           — stderr from each child process
#
# Result codes written to CSV:
#   VERIFIED   m first_n p_first elapsed_s
#   MISMATCH   m first_n p_first elapsed_s
#   NOT_FOUND  m elapsed_s
#   ERROR      m <reason>
#
# =============================================================================

import subprocess
import sys
import os
import time
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

# ---------------------------------------------------------------------------
# Case definitions — generated from ext1_results.csv + ext2_results.csv
# Fields: (m, p_m, a_m, p_a_m, cramér_ratio, group)
# ---------------------------------------------------------------------------

CASES = [
    # --- Group E2: a_m in [2e9, 1e10) — 15 uniformly spaced cases ---
    (102283, 1331951,  2029876287,   47790310313,  1720.4, "E2"),
    ( 92335, 1191563,  2223036115,   52549379809,  2105.8, "E2"),
    ( 36726,  437083,  2412333964,   57230326859,  6249.0, "E2"),
    (114161, 1500061,  2744916180,   65491493099,  2064.7, "E2"),
    ( 20686,  233173,  2858183263,   68314916771, 13904.3, "E2"),
    ( 61007,  760367,  3353229683,   80707333381,  4988.3, "E2"),
    (100871, 1311857,  3893067251,   94307911091,  3349.8, "E2"),
    ( 36303,  431833,  4006056395,   97164894163, 10509.9, "E2"),
    ( 19234,  215317,  5206723528,  127712739697, 27442.4, "E2"),
    ( 97351, 1262453,  5389653243,  132394050307,  4820.0, "E2"),
    (109430, 1433101,  5807064035,  143100102157,  4573.5, "E2"),
    (108013, 1412779,  6369975757,  157587051887,  5088.4, "E2"),
    ( 76761,  975497,  6814926932,  169075281041,  7892.7, "E2"),
    (100528, 1307083,  7675012978,  191366007089,  6628.4, "E2"),
    ( 97251, 1260901,  7850344255,  195922614071,  7028.5, "E2"),

    # --- Group E3: a_m in [1e10, 8.2e10] — all 9 available cases ---
    (112360, 1474519, 11071138165,  280276166167,  8472.7, "E3"),
    ( 51816,  636277, 11824762849,  300167153191, 21022.3, "E3"),
    ( 55468,  685453, 12790824118,  325737648979, 21110.2, "E3"),
    (109857, 1438973, 14929726016,  382615207459, 11708.6, "E3"),
    ( 69546,  876191, 17485936480,  451005863261, 22550.3, "E3"),
    ( 43297,  523007, 20932270152,  543819574261, 45285.2, "E3"),
    ( 37949,  453181, 23612836020,  616425195857, 59012.3, "E3"),
    ( 48834,  596419, 23630637176,  616908416777, 44821.1, "E3"),
    ( 37253,  443873, 29834694588,  786117161153, 76088.3, "E3"),
]

# ---------------------------------------------------------------------------
# Worker: run one verify_min_ext invocation and return result dict
# ---------------------------------------------------------------------------

def run_one(case, verifier, logfile):
    m, pm, am, pam, cr, group = case
    cmd = [verifier, str(m), str(pm), str(am), str(pam)]
    t0 = time.monotonic()

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800    # 30 minutes hard limit per case
        )
        stdout = proc.stdout.strip()
        stderr = proc.stderr.strip()
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t0
        with open(logfile, "a") as f:
            f.write(f"\n--- TIMEOUT m={m} after {elapsed:.0f}s ---\n")
        return {"group": group, "m": m, "am": am, "cr": cr,
                "status": "ERROR", "detail": "TIMEOUT", "elapsed": elapsed}
    except Exception as e:
        elapsed = time.monotonic() - t0
        return {"group": group, "m": m, "am": am, "cr": cr,
                "status": "ERROR", "detail": str(e), "elapsed": elapsed}

    # Log stderr
    with open(logfile, "a") as f:
        f.write(f"\n--- m={m} group={group} ---\n{stderr}\n")

    elapsed = time.monotonic() - t0

    # Parse stdout: "STATUS m [n] [p] elapsed"
    tokens = stdout.split()
    if not tokens:
        return {"group": group, "m": m, "am": am, "cr": cr,
                "status": "ERROR", "detail": "empty stdout", "elapsed": elapsed}

    status = tokens[0]
    return {"group": group, "m": m, "am": am, "cr": cr,
            "status": status, "stdout": stdout, "elapsed": elapsed}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Batch E2/E3 verifier for A247975")
    parser.add_argument("--jobs",     type=int, default=4,
                        help="Max parallel processes (default: 4)")
    parser.add_argument("--verifier", default="./verify_min_ext",
                        help="Path to verify_min_ext binary")
    args = parser.parse_args()

    verifier = args.verifier
    if not os.path.isfile(verifier):
        print(f"ERROR: verifier not found at '{verifier}'", file=sys.stderr)
        sys.exit(1)

    outcsv = "verification_min_ext_results.csv"
    logfile = "run_verify_ext_batch.log"

    # Write CSV header (overwrite if exists)
    with open(outcsv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "m", "a_m", "cramér_ratio", "status",
                    "result_line", "elapsed_s"])

    with open(logfile, "w") as f:
        f.write(f"run_verify_ext_batch — started {datetime.now()}\n")
        f.write(f"verifier: {verifier}  jobs: {args.jobs}\n")

    print(f"[batch] Starting {len(CASES)} cases  "
          f"(jobs={args.jobs}, verifier={verifier})")
    print(f"[batch] Output CSV : {outcsv}")
    print(f"[batch] Log file   : {logfile}")
    print()

    t_start = time.monotonic()
    done = 0
    verified = 0
    errors = 0

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_one, c, verifier, logfile): c for c in CASES}

        for fut in as_completed(futures):
            res = fut.result()
            done += 1
            ok = (res["status"] == "VERIFIED")
            if ok:
                verified += 1
            if res["status"] in ("ERROR", "MISMATCH", "NOT_FOUND"):
                errors += 1

            # Append to CSV immediately (thread-safe via GIL for small writes)
            with open(outcsv, "a", newline="") as f:
                w = csv.writer(f)
                w.writerow([
                    res["group"],
                    res["m"],
                    res["am"],
                    f"{res['cr']:.1f}",
                    res["status"],
                    res.get("stdout", res.get("detail", "")),
                    f"{res['elapsed']:.2f}",
                ])

            elapsed_total = time.monotonic() - t_start
            flag = "✓" if ok else "✗ " + res["status"]
            print(f"[{done:2d}/{len(CASES)}]  m={res['m']:6d}  "
                  f"group={res['group']}  {flag}  "
                  f"{res['elapsed']:.1f}s  (total {elapsed_total:.0f}s)")
            sys.stdout.flush()

    elapsed_total = time.monotonic() - t_start
    print()
    print(f"[batch] Done: {verified}/{len(CASES)} VERIFIED, "
          f"{errors} errors/mismatches  —  total {elapsed_total:.0f}s")
    print(f"[batch] Results in: {outcsv}")

    sys.exit(0 if errors == 0 else 1)

if __name__ == "__main__":
    main()
