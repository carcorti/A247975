# A247975 — Computational Study of Sun's Conjecture 4.1(i)

**Repository for the paper:**  
*Exceptional Cases and Statistical Structure of Sequence A247975:  
a Computational Study of Sun's Conjecture 4.1(i)*  
Carlo Corti — in preparation for submission to *Mathematics of Computation* (AMS), 2026.

> **arXiv:** [to be added upon deposit]  
> **Zenodo DOI:** https://doi.org/10.5281/zenodo.18920372  
> **OEIS sequence:** [A247975](https://oeis.org/A247975)

---

## The Conjecture

Sun's Conjecture 4.1(i) (Zhi-Wei Sun, 2013/2019) states:

> For every positive integer *m*, there exists a positive integer *n* such that  
> (*m* + *n*) divides *p*_m² + *p*_n²,  
> where *p*_k denotes the *k*-th prime.

The sequence A247975 is defined as *a*(*m*) = the least such *n*.

---

## What This Repository Contains

### Source code

| File | Description |
|------|-------------|
| `src/sun41_v6b.c` | Main search program (Phase 6). Computes *a*(*m*) for *m* ≤ *M* with search bound *n* ≤ *N*. Uses a segmented sieve of Eratosthenes with a modular obstruction filter. **This is the program that produced the initial dataset (Phase 6, *n* ≤ 2 × 10⁹).** |
| `src/sun41_v6.c` | Earlier version (Phase 5 and below). Differs from v6b only in the type of the prime counter (`int` vs `int64_t`). Included for completeness and reproducibility of earlier phases. |
| `src/sun41_v12.c` – `src/sun41_v12h.c` | Extended-search programs (Ext-1 through Ext-8). Family of targeted segmented-sieve programs that pushed the search bound from *n* ≤ 2 × 10⁹ to *n* ≤ 1.6 × 10¹⁰, resolving 36 of the 49 initially unresolved cases. Peak RAM < 100 MB per run. |
| `src/verify_min.c` | Independent minimality verifier. Uses a completely separate sieve implementation with no shared code with `sun41_v6b.c`. Verified 290 cases (Groups B1, B2, C) by exhaustive search from *n* = 1. Subject to future extension for the 36 extended-search cases. |

### Python/SymPy verification

| File | Description |
|------|-------------|
| `scripts/run_strategy_BC.py` | Orchestrator for the independent minimality verification campaign (Groups B1, B2, C). Drives `verify_min` and collects results into CSV and log. |

### Data

| File | Description |
|------|-------------|
| `data/fase6_m100000_b2e9.csv` | Full output of the Phase 6 search: *m* = 1…100 000, bound *n* ≤ 2 × 10⁹. Columns: *m*, *p*_m, *a*(*m*), *p*_{*a*(*m*)}, *s* = *m* + *a*(*m*), obstruction rate. Contains 99 951 resolved values and 49 unresolved cases. |
| `data/Ext-1.csv` – `data/Ext-7.csv` | Output of the extended-search runs Ext-1 through Ext-7 (targeted search on the 49 unresolved cases). Each file contains the status (RESOLVED / NOT_FOUND) and, where resolved, the value *a*(*m*) found. |
| `data/b247975.txt` | B-file for OEIS A247975: consecutive values *a*(*n*) for *n* = 1…22 801 (format: `n a(n)`). Updated to reflect the 36 new values from the extended search; the consecutive run stops at *n* = 22 801 (the first remaining unresolved index). |
| `data/a247975.txt` | Extended dataset: all 99 987 resolved values in OEIS a-file format, covering *m* = 1…100 000 with 13 gaps. Updated to reflect the extended search. |
| `data/verification_results_20260305_075834.csv` | Results of the 290 independent minimality verifications (Groups B1, B2, C). Columns: group, *m*, *p*_m, *a*(*m*), *p*_{*a*(*m*)}, status, elapsed time. |
| `data/verification_log_20260305_075834.txt` | Detailed log of the verification run. |

### Logs

| File | Description |
|------|-------------|
| `logs/Ext-1.txt` – `logs/Ext-7.txt` | Execution logs for extended-search runs Ext-1 through Ext-7, including progress, timing, and per-case FOUND/OPEN status. |

### Paper

| File | Description |
|------|-------------|
| `paper/A247975_rev20.tex` | LaTeX source of the submitted manuscript (Rev 20, in preparation for submission to *Mathematics of Computation*). |

---

## Key Results

- **99 987** values of *a*(*m*) determined for *m* = 1…100 000 (final search bound *n* ≤ 1.6 × 10¹⁰).
- **13** values remain unresolved within the final bound; all have estimated Cramér ratios > 14 000.
  The 13 unresolved indices are: 22802, 23995, 37249, 37253, 37949, 43297, 48834, 53963, 66257, 69546, 76868, 82300, 98379.
  An arithmetic cluster {37 249, 37 253, 37 949} and two extreme cases (*m* = 22 802, *m* = 23 995, ratios > 65 000) are of particular interest.
- **434** exceptional resolved cases (*a*(*m*) > 10⁸).
- **Divisibility** independently verified for all 434 exceptional cases via SymPy (398 from the initial computation + 36 from the extended search).
- **Minimality** independently verified for 290 cases via a separate C sieve (`verify_min.c`): 100 exceptional cases, 150 randomly selected non-exceptional cases, and all 40 cases with *a*(*m*) > 10⁹ in the initial dataset. Independent minimality verification for the 36 extended-search cases is planned as future work.
- Global maximum Cramér-type ratio: *m* = 11 924, ratio ≈ 54 679 (new record; previous record *m* = 4703, ratio ≈ 19 111).
- Largest absolute value: *a*(55 468) = 12 790 824 118.
- Hill tail-index estimate: α̂ ≈ 0.713 (heavy-tailed distribution).

---

## Building and Running

### Requirements

- C compiler with C99 support (tested: GCC 11+, Clang 14+)
- OpenMP support (for the extended-search programs `sun41_v12*.c`)
- ~21 GB RAM for the full Phase 6 run (*M* = 100 000, *N* = 2 × 10⁹)
- < 100 MB RAM for each extended-search run (`sun41_v12*.c`)
- Python 3.8+ and SymPy ≥ 1.12 for the verification scripts

### Compile the main search program

```bash
gcc -O2 -std=c99 -o sun41_v6b src/sun41_v6b.c -lm
```

### Run (full Phase 6 parameters)

```bash
./sun41_v6b 100000 2000000000 51409198278
# Arguments: M_MAX  SEARCH_BOUND  SIEVE_LIMIT
# SIEVE_LIMIT should be ~9% above p_{SEARCH_BOUND} to avoid edge effects.
# Warning: requires ~21 GB RAM and ~10 minutes on a modern workstation.
```

### Run (small test, M=1000, bound=10^7)

```bash
./sun41_v6b 1000 10000000 230000000
# Completes in seconds; requires ~1.5 GB RAM.
```

### Compile and run the extended-search programs

```bash
gcc -O2 -std=c99 -fopenmp -o sun41_v12h src/sun41_v12h.c -lm
# Each sun41_v12*.c targets a specific subset of unresolved cases;
# see comments in each source file for the intended bound and case list.
```

### Compile and run the independent verifier

```bash
gcc -O2 -std=c99 -o verify_min src/verify_min.c -lm
# Then run via the Python orchestrator:
python3 scripts/run_strategy_BC.py
```

---

## Reproducing the Results

### Initial dataset (Phase 6)

The file `data/fase6_m100000_b2e9.csv` was produced by:

```bash
./sun41_v6b 100000 2000000000 51409198278
```

on an AMD Ryzen 9 7940HS workstation with 32 GB RAM.  
Elapsed time: 597 seconds (≈ 10 minutes). Peak RAM: ≈ 20.8 GB.  
The sieve uses deterministic primality (no probabilistic tests); output is fully reproducible.

Cross-consistency check: all 10 000 values for *m* ≤ 10 000 agree exactly between Phase 5 and Phase 6 outputs.

### Extended search (Ext-1 through Ext-8)

The extended-search runs used `sun41_v12.c`–`sun41_v12h.c` on an AMD Ryzen 9 7940HS (UM790 Pro, 16 OpenMP threads).  
Total elapsed time: approximately 129 minutes. Peak RAM per run: < 100 MB.  
Runs Ext-1 through Ext-7 each resolved at least one previously open case; their outputs (CSV and execution log) are archived in `data/Ext-1.csv`–`data/Ext-7.csv` and `logs/Ext-1.txt`–`logs/Ext-7.txt` respectively.
Run Ext-8 (`sun41_v12h.c, bound n ≤ 1.6 × 10¹⁰`) resolved no additional cases: all 13 remaining candidates exceeded the bound without a solution. As it produced no new values of a(m), its output is not archived separately; the source code `src/sun41_v12h.c` is included for completeness and to allow independent reproduction of this negative result.

---

## Citation

If you use this code or data, please cite:

```
Carlo Corti,
"Exceptional Cases and Statistical Structure of Sequence A247975:
 a Computational Study of Sun's Conjecture 4.1(i)",
in preparation for submission to Mathematics of Computation (AMS), 2026.
arXiv: [to be added upon deposit]
Zenodo DOI: https://doi.org/10.5281/zenodo.18920372
```

---

## License

The source code in `src/` and `scripts/` is released under the **MIT License** (see `LICENSE`).  
The data files in `data/` are released under **CC0 1.0 Universal** (public domain dedication).

---

## Acknowledgments

The author thanks Zhi-Wei Sun for the beautiful conjecture that motivated this work,  
and Chai Wah Wu for the original OEIS b-file (m ≤ 10 000) used as a cross-check.  
The computational implementation, including the C search program and verification sieve,  
was developed with extensive assistance from Anthropic's Claude language models  
(`claude-sonnet-4-6` and `claude-opus-4-6`); Claude also assisted substantially in the  
statistical analysis and in the preparation of the manuscript.  
All computations were executed and independently verified by the author,  
who bears full scientific responsibility for the results.
