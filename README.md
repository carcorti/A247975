# A247975 — Computational Study of Sun's Conjecture 4.1(i)

**Repository for the paper:**  
*Exceptional Cases and Statistical Structure of Sequence A247975:  
a Computational Study of Sun's Conjecture 4.1(i)*  
Carlo Corti — submitted to the *Journal of Integer Sequences*, 2026.

> **arXiv:** [to be added upon deposit]  
> **Zenodo DOI:** [to be added upon release]  
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
| `src/sun41_v6b.c` | Main search program (Phase 6). Computes *a*(*m*) for *m* ≤ *M* with search bound *n* ≤ *N*. Uses a segmented sieve of Eratosthenes with a modular filter. **This is the program that produced the dataset.** |
| `src/sun41_v6.c` | Earlier version (Phase 5 and below). Differs from v6b only in the type of the prime counter (`int` vs `int64_t`). Included for completeness and reproducibility of earlier phases. |
| `src/verify_min.c` | Independent minimality verifier. Uses a completely separate sieve implementation with no shared code with `sun41_v6b.c`. Verified 290 cases (Groups B1, B2, C) by exhaustive search from *n* = 1. |

### Python/SymPy verification

| File | Description |
|------|-------------|
| `scripts/run_strategy_BC.py` | Orchestrator for the independent minimality verification campaign (Groups B1, B2, C). Drives `verify_min` and collects results. |

### Data

| File | Description |
|------|-------------|
| `data/fase6_m100000_b2e9.csv` | Full output of the Phase 6 search: 99 951 resolved values and 49 unresolved cases for *m* = 1…100 000, with *p*_m, *a*(*m*), *p*_{*a*(*m*)}, *s* = *m* + *a*(*m*), and obstruction rate. |
| `data/b247975.txt` | B-file for OEIS A247975: consecutive values *a*(*n*) for *n* = 1…11 923 (format: `n a(n)`). |
| `data/a247975.txt` | Extended dataset: all 99 951 resolved values in OEIS a-file format. |
| `data/verification_results_20260305_075834.csv` | Results of the 290 independent minimality verifications (Groups B1, B2, C). Columns: group, *m*, *p*_m, *a*(*m*), *p*_{*a*(*m*)}, status, elapsed time. |
| `data/verification_log_20260305_075834.txt` | Detailed log of the verification run. |

### Paper

| File | Description |
|------|-------------|
| `paper/A247975_rev15.tex` | LaTeX source of the submitted manuscript (Rev 15, final). |

---

## Key Results

- **99 951** values of *a*(*m*) determined for *m* = 1…100 000 (search bound *n* ≤ 2 × 10⁹).
- **49** values remain unresolved within this bound; all exhibit filter-detected obstruction rates of 75.1%–77.3%.
- **398** exceptional resolved cases (*a*(*m*) > 10⁸); **290** independently verified for minimality.
- Global maximum Cramér-type ratio: *m* = 4703, *a*(4703) = 760 027 770, ratio ≈ 19 111.
- Largest resolved value: *a*(14 740) = 1 994 463 433.
- Hill tail-index estimate: α̂ ≈ 0.713 (heavy-tailed distribution).

---

## Building and Running

### Requirements

- C compiler with C99 support (tested: GCC 11+, Clang 14+)
- ~21 GB RAM for the full Phase 6 run (*M* = 100 000, *N* = 2 × 10⁹)
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

### Compile and run the independent verifier

```bash
gcc -O2 -std=c99 -o verify_min src/verify_min.c -lm
# Then run via the Python orchestrator:
python3 scripts/run_strategy_BC.py
```

---

## Reproducing the Results

The dataset `data/fase6_m100000_b2e9.csv` was produced by:

```bash
./sun41_v6b 100000 2000000000 51409198278
```

on an AMD Ryzen 9 7940HS workstation with 32 GB RAM.  
Elapsed time: 597 seconds (≈ 10 minutes).  
The sieve uses deterministic primality (no probabilistic tests); output is fully reproducible.

Cross-consistency check: all 10 000 values for *m* ≤ 10 000 agree exactly between Phase 5 and Phase 6 outputs.

---

## Citation

If you use this code or data, please cite:

```
Carlo Corti,
"Exceptional Cases and Statistical Structure of Sequence A247975:
 a Computational Study of Sun's Conjecture 4.1(i)",
Journal of Integer Sequences, 2026.
arXiv: [to be added]
Zenodo DOI: [to be added]
```

---

## License

The source code in `src/` and `scripts/` is released under the **MIT License** (see `LICENSE`).  
The data files in `data/` are released under **CC0 1.0 Universal** (public domain dedication).

---

## Acknowledgments

The author thanks Zhi-Wei Sun for the conjecture that motivated this work,  
and Chai Wah Wu for the original OEIS b-file (m ≤ 10 000) used as a cross-check.  
Computational assistance by Anthropic's Claude language models (claude-sonnet-4-6, claude-opus-4-6).
