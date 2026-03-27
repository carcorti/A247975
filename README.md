# A247975 — Computational Study of Sun's Conjecture 4.1(i)

**Repository for the paper:**  
*Exceptional Cases and Statistical Structure of Sequence A247975:  
a Computational Study of Sun's Conjecture 4.1(i)*  
Carlo Corti — manuscript in preparation, 2026.

> **arXiv:** [to be added upon deposit]  
> **Zenodo DOI:** [10.5281/zenodo.18920371](https://doi.org/10.5281/zenodo.18920371)  
> **OEIS sequence:** [A247975](https://oeis.org/A247975)

---

## The Conjecture

Sun's Conjecture 4.1(i) (Zhi-Wei Sun, 2013/2019) states:

> For every positive integer *m*, there exists a positive integer *n* such that  
> (*m* + *n*) divides *p*\_m² + *p*\_n²,  
> where *p*\_k denotes the *k*-th prime.

The sequence A247975 is defined as *a*(*m*) = the least such *n*.

---

## What This Repository Contains

### Source code (`src/`)

| File | Description |
|---|---|
| `sun_v6.c` | **Phase 1 search program.** Computes *a*(*m*) for *m* ≤ 120 000, search bound *n* ≤ 2 × 10⁹. Uses a global sieve of Eratosthenes with a modular filter. **This is the primary program that produced the Phase 1 dataset.** Requires ~21 GB RAM. |
| `sun_ext_v11.c` | **Extended search program.** Segmented sieve with no persistent prime table; peak RAM < 100 MB regardless of bound. Reads a list of unresolved *m*-values from file and searches up to a specified *n*-bound. Used in four sessions (S1–S4) to extend the bound from 2 × 10⁹ to 2 × 10¹¹, resolving 64 additional cases. |
| `verify_min_ext.c` | **Independent minimality verifier (Extended).** Segmented sieve with no shared code with the search programs. Uses `__uint128_t` arithmetic to handle *p*\_n² overflow at large *n*. Verified 29 Extended-search cases by exhaustive search from *n* = 1. |
| `sun41_v6b.c` | Earlier Phase 6 program (v1, *m* ≤ 100 000). Retained for reproducibility of v1 results. |
| `sun41_v6.c` | Earlier Phase 5 program (v1). Retained for completeness. |
| `verify_min.c` | Independent minimality verifier for Phase 1 cases (v1). Verified 290 cases (Groups B1, B2, C). |

### Python scripts (`scripts/`)

| File | Description |
|---|---|
| `verify_ext_sympy.py` | Verifies divisibility of all 64 Extended-search cases via `sympy.prime()`, independently of the C sieve. All 64 verified without error. |
| `run_verify_ext_batch.py` | Orchestrator for the Extended minimality verification campaign (Groups E2, E3). Drives `verify_min_ext` and collects results into `verification_min_ext_results.csv`. |
| `plot_a247975.py` | Generates the logarithmic scatterplot of A247975(n) for *n* = 1..120 000 (`a247975_2.png`), suitable for OEIS upload. Requires matplotlib and numpy. |
| `run_strategy_BC.py` | Orchestrator for the Phase 1 minimality verification campaign (Groups B1, B2, C) (v1). |

### Data (`data/`)

| File | Description |
|---|---|
| `phase1_results.csv` | **Phase 1 output:** 119 931 resolved values and 69 unresolved cases for *m* = 1…120 000. Columns: *m*, *p*\_m, *a*(*m*), *p*\_{*a*(*m*)}, *s* = *m* + *a*(*m*), Cramér ratio. Unresolved rows have `UNRESOLVED` in the *a*(*m*) column. |
| `ext1_results.csv` | Extended search S1 output (bound 2 × 10⁹ → 5 × 10¹⁰): 62 resolved, 7 still open. |
| `ext2_results.csv` | Extended search S2 output (bound 5 × 10¹⁰ → 10¹¹): 2 resolved, 5 still open. |
| `ext3_results.csv` | Extended search S3 output (bound 10¹¹ → 1.5 × 10¹¹): 0 resolved, 5 still open. |
| `ext4_results.csv` | Extended search S4 output (bound 1.5 × 10¹¹ → 2 × 10¹¹): 0 resolved. **These 5 are the definitive open cases.** |
| `a247975_v2.txt` | **OEIS a-file (v2):** 119 995 resolved values and 5 unresolved (marked `?`) for *n* = 1…120 000, in standard OEIS format (`n a(n)`). |
| `verification_ext_results.csv` | SymPy verification output for all 64 Extended-search cases. All 64 verified (status: `OK`). |
| `verification_min_ext_results.csv` | Independent sieve verification output for Extended campaigns E2 (15 cases) and E3 (9 cases). All 24 verified (status: `VERIFIED`). |
| `fase6_m100000_b2e9.csv` | v1 Phase 6 output (*m* ≤ 100 000). Retained for reproducibility. |
| `b247975.txt` | OEIS b-file (Chai Wah Wu, v1): consecutive *a*(*n*) for *n* = 1…11 923. |
| `a247975.txt` | OEIS a-file (v1): 99 987 resolved values for *n* = 1…100 000. |
| `verification_results_*.csv` | v1 Phase 1 minimality verification results (Groups B1, B2, C). |
| `verification_log_*.txt` | v1 verification run log. |

### Paper (`paper/`)

| File | Description |
|---|---|
| `A247975_rev25.tex` | **Current manuscript** (Rev 25, March 2026). LaTeX source, `amsart` class. |
| `A247975_rev15.tex` | v1 manuscript (Rev 15). Retained for reference. |

---

## Key Results (v2)

- **119 995** values of *a*(*m*) determined for *m* = 1…120 000.
- **Phase 1** (bound *n* ≤ 2 × 10⁹): 119 931 resolved, 69 unresolved.
- **Extended search** (bound *n* ≤ 2 × 10¹¹, four sessions): 64 of the 69 resolved.
- **5 definitive open cases** (*m* ∈ {37 249, 66 257, 76 868, 98 379, 117 862}), all with estimated Cramér ratio lower bound > 145 000; proposed as explicit open problems.
- **606** exceptional resolved cases (*a*(*m*) > 10⁸); all independently verified for divisibility via SymPy.
- **319** cases independently verified for minimality (290 Phase 1 + 29 Extended) via `verify_min.c` / `verify_min_ext.c`.
- **10 925** individual verifications total, zero discrepancies.

### Top 5 by Cramér-type ratio *C*(*m*) = *a*(*m*) / (*m* log *m*)

| Rank | *m* | *a*(*m*) | *p*\_*m* | *C*(*m*) |
|---:|---:|---:|---:|---:|
| 1 | 23 995 | 47 572 503 750 | 274 453 | 196 577 |
| 2 | 22 802 | 44 000 698 448 | 259 577 | 192 303 |
| 3 | 82 300 | 81 961 045 941 | 1 052 459 | 87 990 |
| 4 | 53 963 | 50 362 759 647 | 665 179 | 85 653 |
| 5 | 37 253 | 29 834 694 588 | 443 873 | 76 088 |

- Largest absolute value: *a*(82 300) = 81 961 045 941.
- Hill tail-index estimate: α̂ ≈ 0.713 (heavy-tailed distribution).
- Median Cramér ratio: ≈ 0.154.

---

## Building and Running

### Requirements

- C compiler with C11 support and OpenMP (tested: GCC 13, `-fopenmp`)
- ≥ 22 GB RAM for the full Phase 1 run (*m* ≤ 120 000, *n* ≤ 2 × 10⁹)
- < 100 MB RAM for the Extended search (`sun_ext_v11`)
- Python 3.8+, SymPy ≥ 1.12, matplotlib, numpy (for scripts)

### Phase 1 — `sun_v6.c`

```bash
# Compile (with PGO for best performance)
gcc -O2 -fopenmp -march=native -o sun_v6 src/sun_v6.c -lm

# Run (full parameters: m≤120000, n≤2e9)
./sun_v6 120000 2000000000 64101957761
# SIEVE_LIMIT = 64101957761 (~9% above p_{2e9} = 47055833459)
# Warning: requires ~22 GB RAM and ~9 minutes on 16-thread Ryzen 9 7940HS.

# Small test (m≤1000, n≤1e7)
./sun_v6 1000 10000000 230000000
```

### Extended search — `sun_ext_v11.c`

```bash
# Compile with PGO (recommended)
gcc -O3 -fopenmp -march=native -fprofile-generate \
    -o sun_ext_v11_pgo src/sun_ext_v11.c -lm
# [run training workload, then:]
gcc -O3 -fopenmp -march=native -fprofile-use \
    -o sun_ext_v11_pgo src/sun_ext_v11.c -lm

# Run Session S1 (from Phase 1 bound to 5e10)
./sun_ext_v11_pgo data/phase1_unresolved.txt 2000000001 50000000000 \
    > data/ext1_results.csv 2> ext1.log

# Regression test (mandatory after any recompilation)
echo "11924" > /tmp/check.txt
./sun_ext_v11_pgo /tmp/check.txt 2000000001 6200000000 \
    > /tmp/reg.csv 2>/tmp/reg.log
grep "FOUND" /tmp/reg.log
# Expected: FOUND m=11924 n=6119832581
```

### Independent verifier — `verify_min_ext.c`

```bash
# Compile
gcc -O2 -o verify_min_ext src/verify_min_ext.c -lm

# Smoke test (mandatory after recompilation)
./verify_min_ext 11924 127289 6119832581 151142662267
# Expected: VERIFIED 11924 6119832581 151142662267  (~149 s)
```

### SymPy verification

```bash
python3 scripts/verify_ext_sympy.py
# Reads ext1_results.csv + ext2_results.csv, verifies all 64 cases.
# Output: data/verification_ext_results.csv
```

### Generate OEIS plot

```bash
python3 scripts/plot_a247975.py
# Output: a247975_2.png  (upload to OEIS as supporting file)
```

---

## Reproducing the Results

### Phase 1

```bash
./sun_v6 120000 2000000000 64101957761
```

Produces `phase1_results.csv`. Platform: AMD Ryzen 9 7940HS (16 threads), ~9 min, ~22 GB RAM.  
Output is fully deterministic (no probabilistic primality tests).

### Extended search (4 sessions)

```bash
# S1
./sun_ext_v11_pgo data/phase1_unresolved.txt  2000000001  50000000000 > ext1_results.csv 2>ext1.log
# S2
./sun_ext_v11_pgo data/ext1_unresolved.txt   50000000001 100000000000 > ext2_results.csv 2>ext2.log
# S3
./sun_ext_v11_pgo data/ext2_unresolved.txt  100000000001 150000000000 > ext3_results.csv 2>ext3.log
# S4
./sun_ext_v11_pgo data/ext3_unresolved.txt  150000000001 200000000000 > ext4_results.csv 2>ext4.log
```

S1: ~49 min, 62 resolved. S2: ~69 min, 2 resolved. S3: ~103 min, 0. S4: ~135 min, 0.

### Critical lesson (bug OBS3-P1, March 2026)

> Any optimization that restructures an inner loop on mathematical computation  
> **must** be verified against known reference values before production use.  
> Static review — even multi-LLM — is not sufficient.  
> **Mandatory regression test after every recompilation of `sun_ext_v11.c`:**  
> expected `FOUND m=11924 n=6119832581`.

---

## Open Problems

The following 5 values of *m* ≤ 120 000 have no solution within the search bound *n* ≤ 2 × 10¹¹:

| *m* | *p*\_*m* | Obstruction rate | *C*(*m*) > |
|---:|---:|---:|---:|
| 37 249 | 443 837 | 76.2% | 510 126 |
| 66 257 | 831 301 | 76.9% | 271 910 |
| 76 868 | 977 087 | 77.0% | 231 280 |
| 98 379 | 1 276 747 | 77.2% | 176 831 |
| 117 862 | 1 552 879 | 77.4% | 145 316 |

For each of these, find the least positive integer *n* such that (*m* + *n*) | *p*\_*m*² + *p*\_*n*², or prove that none exists.

---

## Citation

If you use this code or data, please cite:

```
Carlo Corti,
"Exceptional Cases and Statistical Structure of Sequence A247975:
 a Computational Study of Sun's Conjecture 4.1(i)",
manuscript in preparation, 2026.
Zenodo DOI: 10.5281/zenodo.18920371
```

---

## License

The source code in `src/` and `scripts/` is released under the **MIT License** (see `LICENSE`).  
The data files in `data/` are released under **CC0 1.0 Universal** (public domain dedication).

---

## Acknowledgements

The author thanks Zhi-Wei Sun for the beautiful conjecture that motivated this work,
and Chai Wah Wu, whose OEIS b-file for A247975 (*m* ≤ 10 000) provided a valuable
cross-check for the Phase 1 results.  
Computational assistance: Anthropic's Claude language models (`claude-sonnet-4-6`, `claude-opus-4-6`).
