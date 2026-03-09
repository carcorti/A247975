/*
 * sun41_v6.c  --  Computational verification of Sun's Conjecture 4.1(i)
 *
 * Conjecture (Z.-W. Sun, 2014-09-29, arXiv:1309.1679):
 *   For every m in Z+, there exists n in Z+ such that (m+n) | (p_m^2 + p_n^2),
 *   where p_k denotes the k-th prime.
 *
 * Also verifies the companion part (ii):
 *   (m+n) | (p_{m^2} + p_{n^2}) for some n in Z+.
 *
 * OEIS references: A247975 (part i), A248354 (part ii).
 * Known hard case: m=4703, min_n = 760027770 (Sun's original note).
 *
 * Compile:
 *   gcc -O3 -march=znver4 -mtune=znver4 -funroll-loops -ffast-math \
 *       -fomit-frame-pointer -flto -fopenmp \
 *       src/sun41_v6.c -o sun41_v6 -lm
 *
 * Usage:
 *   ./sun41_v4 [M_MAX] [SEARCH_BOUND] [SIEVE_LIMIT]
 *   Defaults: M_MAX=10000, SEARCH_BOUND=10^8, SIEVE_LIMIT=2*10^8
 *
 * Output: CSV on stdout, progress/summary on stderr
 *
 * CHANGES v6 (dynamic prime table — enables sieve up to ~2×10⁹):
 *   v5 had MAX_PRIME_SLOTS = 12,000,000 (hardcoded), sufficient only for
 *   sieve ≤ 2×10⁸ (π(2×10⁸) ≈ 11.08M primes).  For Fase 2 of the
 *   experimental design (bound n ≤ 10⁸, requiring sieve ≤ 2.1×10⁹),
 *   the prime table needs ~98M slots = 784 MB.  A static array of that
 *   size would always occupy RAM even for small runs.
 *   Fix: prime table is now allocated dynamically after the sieve is built.
 *   The exact count is known at that point (g_prime_count after first pass),
 *   so we allocate precisely what is needed with a small safety margin.
 *   Two-pass sieve strategy:
 *     Pass 1: count primes only (no prime table written).
 *     Pass 2: fill prime table into freshly allocated array.
 *   Added: consistency check — if BOUND > π(SIEVE_LIMIT), BOUND is clamped
 *   and a warning is printed (instead of silent wrong results).
 *   RAM usage:
 *     sieve=2×10⁸ (Fase 1): bitset 12.5 MB + table  84 MB =  ~97 MB
 *     sieve=2.1×10⁹ (Fase 2): bitset 262 MB + table 784 MB = ~1.05 GB
 *
 * CHANGES v5 (parity mask for pm):
 *   In modular_obstruction_filter(), the check (pm % q == 0) was executed
 *   for every candidate s.  Since p_m ∈ {3,7,11,19} only for m ∈ {2,4,5,8},
 *   this check is almost always a no-op — but still costs 4 divisions per s.
 *   Fix: precompute a bitmask 'pm_qmask' (4 bits) once per m, where bit i is
 *   set iff p_m == FILTER_Q[i].  The inner loop skips the division entirely
 *   when bit i is clear (the common case for m > 8).
 *   This saves 4 integer divisions per s for ~99.96% of all m values.
 *   The wheel residue table (Q = 4389) suggested by the reviewer was evaluated
 *   but not implemented: it yields 42.9% static blocking vs. the current ~49%
 *   dynamic blocking with early exit, and requires a per-m table rebuild.
 *   Net gain would be marginal; parity mask gives better return per line of code.
 *
 * CHANGES v4 (modular wheel filter):
 *   Added a fast obstruction pre-filter based on small primes q ≡ 3 (mod 4).
 *   For each candidate s = m+n, before accessing the prime table, we test
 *   whether s has a prime factor q ∈ {3, 7, 11, 19} with odd valuation
 *   and q ∤ p_m.  Such a q makes p_n^2 ≡ -p_m^2 (mod q) unsolvable
 *   (since -1 is not a QR mod q).  This filter is mathematically exact
 *   (no false positives) and eliminates ~49% of candidates in practice,
 *   giving an empirical ~1.9x speedup before any prime-table access.
 *   See modular_obstruction_filter() for full documentation.
 *
 * CHANGES v3:
 *   1. Bug fix: filter used DEFAULT_SIEVE_LIMIT instead of runtime g_sieve_limit.
 *   2. has_partial_obstruction: residual cofactor checked for primality via bsieve.
 *   3. max_n_ii precomputed in main(), passed as parameter to find_min_n_4ii().
 *   4. Micro-opt: removed redundant inner modulo in divisibility check.
 *
 * CHANGES v2:
 *   1. Bitset sieve (~25 MB vs 200 MB).
 *   2. Loop restructured on s = m+n.
 *   3. Filter: prime s ≡ 3 (mod 4) → skip.
 *   4. Modular reduction before multiplication.
 *   5. Off-by-one fix in part (ii).
 *   6. has_obstruction() → has_partial_obstruction(), documented as heuristic.
 *   7. Progress monitor outside parallel region.
 *
 * Memory: ~25 MB bitset + ~90 MB prime table ≈ 115 MB total.
 *
 * Author:  [your name]
 * Date:    2026
 * License: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ------------------------------------------------------------------ */
/*  Configuration defaults                                             */
/* ------------------------------------------------------------------ */
#define DEFAULT_M_MAX        10000
#define DEFAULT_BOUND        100000000L
#define DEFAULT_SIEVE_LIMIT  200000000L
/* MAX_PRIME_SLOTS removed: prime table is now allocated dynamically.  */
/* Safety margin applied at allocation time: actual_count * 1.01 + 64 */

/* ------------------------------------------------------------------ */
/*  Global state                                                        */
/* ------------------------------------------------------------------ */
static uint8_t  *g_bsieve      = NULL;
static int64_t  *g_prime       = NULL;
static int       g_prime_count = 0;
static long      g_sieve_limit = DEFAULT_SIEVE_LIMIT;

/* ------------------------------------------------------------------ */
/*  Bitset sieve (one bit per odd integer >= 3)                        */
/*  Mapping: odd k = 2*i+3  <-->  bit i  (i = (k-3)/2)               */
/*  Bit value: 0 = prime, 1 = composite.                               */
/* ------------------------------------------------------------------ */
#define BSIEVE_IDX(k)    (((k) - 3) >> 1)
#define BSIEVE_BYTES(l)  ((((l) - 3) >> 1) / 8 + 2)

static inline int bsieve_is_prime(long k)
{
    if (k < 2)       return 0;
    if (k == 2)      return 1;
    if (k % 2 == 0)  return 0;
    long i = BSIEVE_IDX(k);
    return !((g_bsieve[i >> 3] >> (i & 7)) & 1);
}

static inline void bsieve_set_composite(long k)
{
    long i = BSIEVE_IDX(k);
    g_bsieve[i >> 3] |= (uint8_t)(1 << (i & 7));
}

static void build_sieve(long limit)
{
    size_t bsz = BSIEVE_BYTES(limit);
    g_bsieve = (uint8_t *)calloc(bsz, 1);
    if (!g_bsieve) {
        fprintf(stderr, "Memory allocation failed (bitset: %.1f MB)\n",
                (double)bsz / 1e6);
        exit(1);
    }
    fprintf(stderr, "  Bitset sieve: %.1f MB\n", (double)bsz / 1e6);

    /* Sieve odd composites */
    for (long i = 3; i * i <= limit; i += 2)
        if (bsieve_is_prime(i))
            for (long j = i * i; j <= limit; j += 2 * i)
                bsieve_set_composite(j);

    /* Pass 1: count primes to determine exact allocation size */
    int count = 1;  /* prime 2 */
    for (long k = 3; k <= limit; k += 2)
        if (bsieve_is_prime(k)) count++;

    /* Allocate prime table with 1% safety margin */
    size_t slots = (size_t)(count * 1.01) + 64;
    g_prime = (int64_t *)malloc(sizeof(int64_t) * slots);
    if (!g_prime) {
        fprintf(stderr, "Memory allocation failed (prime table: %.1f MB)\n",
                (double)(slots * sizeof(int64_t)) / 1e6);
        exit(1);
    }
    fprintf(stderr, "  Prime table: %d primes, %.1f MB allocated\n",
            count, (double)(slots * sizeof(int64_t)) / 1e6);

    /* Pass 2: fill prime table */
    g_prime[g_prime_count++] = 2;
    for (long k = 3; k <= limit; k += 2)
        if (bsieve_is_prime(k))
            g_prime[g_prime_count++] = k;
}

/* ------------------------------------------------------------------ */
/*  Modular obstruction filter (NEW in v4)                             */
/*                                                                      */
/*  Returns 1 if s is provably NOT a solution for ANY n with m+n=s,   */
/*  i.e. p_n^2 ≡ -p_m^2 (mod s) has no solution.                     */
/*  Returns 0 if s is NOT ruled out (candidate must be tested).        */
/*                                                                      */
/*  Mathematical basis:                                                 */
/*  For a prime power q^e || s (meaning q^e | s, q^{e+1} ∤ s):        */
/*    If q ≡ 3 (mod 4), e is odd, and q ∤ p_m:                        */
/*      then -p_m^2 is not a QR mod q, so p_n^2 ≡ -p_m^2 (mod q)     */
/*      has no solution, and hence no solution exists mod s.           */
/*                                                                      */
/*  Small primes tested: q ∈ {3, 7, 11, 19} (all ≡ 3 mod 4).         */
/*  These were chosen for maximum filtering with minimal overhead.     */
/*  Adding q=23 gives <1% extra gain; not worth the extra division.    */
/*                                                                      */
/*  This filter is also applied for prime s ≡ 3 (mod 4) (subsumes     */
/*  the v2 prime-s filter when s ∈ {3,7,11,19} or s is detected       */
/*  as prime by bsieve).  The bsieve prime check is retained for       */
/*  prime s outside the small set.                                     */
/*                                                                      */
/*  Empirical performance (measured over 10^5 consecutive s):          */
/*    ~49% of candidates eliminated before prime-table access.         */
/*    ~1.9x speedup on the inner loop.                                 */
/*    Zero false positives (mathematically exact).                     */
/*                                                                      */
/*  The v_q computation uses fast modular checks:                      */
/*    q  | s but q^2 ∤ s  →  v_q(s) = 1 (odd)  →  O(2 divisions)    */
/*    q^2 | s              →  full valuation loop  (rare)              */
/* ------------------------------------------------------------------ */

/* Small primes q ≡ 3 (mod 4) used for the filter */
#define N_FILTER_PRIMES 4
static const int64_t FILTER_Q[N_FILTER_PRIMES] = {3, 7, 11, 19};
/* Precomputed q^2 for fast first-level check */
static const int64_t FILTER_Q2[N_FILTER_PRIMES] = {9, 49, 121, 361};

/*
 * Returns the 2-adic valuation of s with respect to q, i.e. v_q(s) mod 2.
 * Returns 1 if v_q(s) is odd, 0 if even (including 0).
 * Fast path: check q | s and q^2 ∤ s first (handles ~83% of q|s cases).
 */
static inline int vq_is_odd(int64_t s, int qi)
{
    int64_t q  = FILTER_Q[qi];
    int64_t q2 = FILTER_Q2[qi];
    if (s % q  != 0) return 0;   /* q ∤ s → v_q = 0, even */
    if (s % q2 != 0) return 1;   /* v_q = 1, odd (most common case) */
    /* v_q >= 2: count fully */
    int e = 0;
    while (s % q == 0) { s /= q; e++; }
    return e & 1;
}

/*
 * Precompute pm_qmask for a given p_m.
 * Bit i is set iff p_m == FILTER_Q[i].
 * For m > 8 (i.e. p_m > 19), the result is always 0.
 * Call once per m; pass to modular_obstruction_filter().
 */
static inline uint8_t compute_pm_qmask(int64_t pm)
{
    uint8_t mask = 0;
    for (int i = 0; i < N_FILTER_PRIMES; i++)
        if (pm == FILTER_Q[i]) mask |= (uint8_t)(1 << i);
    return mask;
}

/*
 * pm_qmask: precomputed bitmask (bit i set iff p_m == FILTER_Q[i]).
 * Avoids 4 divisions per s in the common case where p_m ∉ {3,7,11,19}.
 */
static inline int modular_obstruction_filter(int64_t s, int64_t pm,
                                              uint8_t pm_qmask)
{
    /* 1. Fast check: small primes q ≡ 3 (mod 4) with odd valuation */
    for (int i = 0; i < N_FILTER_PRIMES; i++) {
        if (pm_qmask & (1 << i)) continue;  /* p_m == q: skip (p_m^2 ≡ 0 mod q) */
        if (vq_is_odd(s, i))
            return 1;                        /* obstruction confirmed */
    }
    /* 2. Retained from v2/v3: prime s ≡ 3 (mod 4) with s ≠ p_m */
    if (s <= g_sieve_limit && bsieve_is_prime((long)s)
            && (s % 4 == 3) && (s != pm))
        return 1;
    return 0;
}

/* ------------------------------------------------------------------ */
/*  Core search part (i)                                                */
/*                                                                      */
/*  Find smallest n in [1, bound] s.t. (m+n) | (p_m^2 + p_n^2).      */
/*  Returns -1 if not found within bound.                              */
/* ------------------------------------------------------------------ */
static inline int64_t find_min_n_4i(int m, int64_t bound)
{
    int64_t pm  = g_prime[m - 1];
    int64_t pm2 = pm * pm;
    int64_t lim = (bound < (int64_t)g_prime_count) ? bound
                                                    : (int64_t)g_prime_count - 1;
    /* Precompute pm_qmask once for this m (avoids pm%q in inner loop) */
    uint8_t pm_qmask = compute_pm_qmask(pm);

    for (int64_t s = (int64_t)m + 1; s <= (int64_t)m + lim; s++) {
        /* Modular obstruction filter: ~49% of s rejected here,
           before any prime-table access. */
        if (modular_obstruction_filter(s, pm, pm_qmask))
            continue;

        int64_t n  = s - m;
        int64_t pn = g_prime[n - 1];

        /* Reduce before multiplying: pnr < s, pnr^2 < s^2 < 4e16 < 2^63 */
        int64_t r   = pm2 % s;
        int64_t pnr = pn  % s;
        if ((r + pnr * pnr) % s == 0)
            return n;
    }
    return -1;
}

/* ------------------------------------------------------------------ */
/*  Core search part (ii)                                               */
/*  Find smallest n s.t. (m+n) | (p_{m^2} + p_{n^2}).                */
/*  max_n_ii: precomputed, largest n s.t. n^2 < g_prime_count.        */
/* ------------------------------------------------------------------ */
static inline int64_t find_min_n_4ii(int m, int64_t bound, int64_t max_n_ii)
{
    int64_t m2 = (int64_t)m * m;
    if (m2 >= (int64_t)g_prime_count) return -2;
    int64_t pm2_val = g_prime[m2 - 1];

    int64_t lim = (bound < max_n_ii) ? bound : max_n_ii;

    for (int64_t n = 1; n <= lim; n++) {
        int64_t n2 = n * n;
        if (n2 >= (int64_t)g_prime_count) break;
        int64_t pn2_val = g_prime[n2 - 1];
        int64_t s = (int64_t)m + n;
        if ((pm2_val % s + pn2_val % s) % s == 0)
            return n;
    }
    return -1;
}

/* ------------------------------------------------------------------ */
/*  Heuristic partial obstruction test (for NOT-FOUND analysis only)   */
/*                                                                      */
/*  Tests prime factors up to p_{999} = 7919.                         */
/*  Residual cofactor: checked for primality via bsieve when possible. */
/*  LIMITATIONS: heuristic, not exhaustive. See full comment in v3.    */
/* ------------------------------------------------------------------ */
static int has_partial_obstruction(int64_t s, int64_t pm)
{
    if (s <= 0) return 0;
    int64_t tmp = s;
    for (int i = 0; i < 1000 && g_prime[i] * g_prime[i] <= tmp; i++) {
        int64_t p = g_prime[i];
        if (tmp % p == 0) {
            int e = 0;
            while (tmp % p == 0) { tmp /= p; e++; }
            if (p != pm && p % 4 == 3 && e % 2 == 1)
                return 1;
        }
    }
    if (tmp > 1 && tmp != pm && tmp % 4 == 3) {
        if (tmp <= g_sieve_limit && bsieve_is_prime((long)tmp))
            return 1;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/*  Main                                                                */
/* ------------------------------------------------------------------ */
int main(int argc, char *argv[])
{
    int     M_MAX = DEFAULT_M_MAX;
    int64_t BOUND = DEFAULT_BOUND;

    if (argc > 1) M_MAX         = atoi(argv[1]);
    if (argc > 2) BOUND         = atoll(argv[2]);
    if (argc > 3) g_sieve_limit = atol(argv[3]);

    /* ---- build sieve and prime table ---- */
    time_t t_start = time(NULL);
    fprintf(stderr, "Building bitset sieve up to %ld ...\n", g_sieve_limit);
    build_sieve(g_sieve_limit);
    fprintf(stderr, "Done: %d primes (largest: %ld), %.1f s\n",
            g_prime_count, g_prime[g_prime_count - 1],
            difftime(time(NULL), t_start));
    if (g_prime_count >= 100000)
        fprintf(stderr, "Spot check: p(1000)=%ld  p(10000)=%ld  p(100000)=%ld\n",
                g_prime[999], g_prime[9999], g_prime[99999]);
    fprintf(stderr, "\n");

    if (BOUND > (int64_t)g_prime_count) {
        fprintf(stderr,
            "WARNING: BOUND=%ld > pi(SIEVE_LIMIT)=%d. "
            "Clamping BOUND to %d.\n"
            "  To search n up to %ld, rerun with SIEVE_LIMIT >= p_{%ld} ~ %.3e\n",
            (long)BOUND, g_prime_count, g_prime_count,
            (long)BOUND, (long)BOUND,
            (double)BOUND * log((double)BOUND));
        BOUND = (int64_t)g_prime_count;
    }

    /* Precompute max_n for part (ii) */
    int64_t max_n_ii = 1;
    while ((max_n_ii + 1) * (max_n_ii + 1) < (int64_t)g_prime_count)
        max_n_ii++;

    /* ---- allocate result arrays ---- */
    int64_t *res_i  = (int64_t *)calloc(M_MAX + 1, sizeof(int64_t));
    int64_t *res_ii = (int64_t *)calloc(M_MAX + 1, sizeof(int64_t));
    if (!res_i || !res_ii) { fprintf(stderr, "Out of memory\n"); exit(1); }

    /* ---- Part (i): parallel search ---- */
    fprintf(stderr, "=== Part (i): m=1..%d, n-bound=%ld ===\n",
            M_MAX, (long)BOUND);

    int found_i = 0, not_found_i = 0;
    int64_t max_n_i = 0, max_m_i = 0;

    #pragma omp parallel for schedule(dynamic, 10) \
            reduction(+:found_i, not_found_i)
    for (int m = 1; m <= M_MAX; m++) {
        if (m - 1 >= g_prime_count) continue;
        int64_t n = find_min_n_4i(m, BOUND);
        res_i[m] = n;
        if (n >= 0) {
            found_i++;
            #pragma omp critical
            { if (n > max_n_i) { max_n_i = n; max_m_i = m; } }
        } else {
            not_found_i++;
        }
    }

    /* ---- Part (ii): parallel search ---- */
    fprintf(stderr, "=== Part (ii): m=1..%d, max_n=%ld ===\n",
            M_MAX, max_n_ii);

    int found_ii = 0, not_found_ii = 0;

    #pragma omp parallel for schedule(dynamic, 10) \
            reduction(+:found_ii, not_found_ii)
    for (int m = 1; m <= M_MAX; m++) {
        int64_t m2 = (int64_t)m * m;
        if (m2 >= (int64_t)g_prime_count) { res_ii[m] = -2; continue; }
        int64_t n = find_min_n_4ii(m, max_n_ii, max_n_ii);
        res_ii[m] = n;
        if (n >= 0)       found_ii++;
        else if (n == -1) not_found_ii++;
    }

    /* ---- CSV output: Part (i) ---- */
    printf("# Sun Conjecture 4.1(i): (m+n) | (p_m^2 + p_n^2)\n");
    printf("# Sieve limit: %ld, Search bound n-index <= %ld\n",
           g_sieve_limit, (long)BOUND);
    printf("# Generated: %s", ctime(&t_start));
    printf("m,p_m,min_n,p_min_n,sum,sum_div,ratio_n_m\n");

    for (int m = 1; m <= M_MAX; m++) {
        if (m - 1 >= g_prime_count) break;
        int64_t n  = res_i[m];
        int64_t pm = g_prime[m - 1];
        if (n >= 0) {
            int64_t pn  = g_prime[n - 1];
            int64_t sum = pm * pm + pn * pn;
            printf("%d,%ld,%ld,%ld,%ld,%ld,%.6f\n",
                   m, pm, n, pn, sum, sum / (m + n), (double)n / m);
        } else {
            printf("%d,%ld,NOT_FOUND,-1,-1,-1,-1\n", m, pm);
        }
    }

    /* ---- Summary ---- */
    fprintf(stderr, "\n=== SUMMARY Part (i) ===\n");
    fprintf(stderr, "m range: 1..%d\n", M_MAX);
    fprintf(stderr, "Search bound: n-index <= %ld", (long)BOUND);
    if (BOUND > 0) fprintf(stderr, " (p_n <= %ld)", g_prime[BOUND - 1]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Found:     %d / %d (%.2f%%)\n",
            found_i, M_MAX, 100.0 * found_i / M_MAX);
    fprintf(stderr, "Not found: %d\n", not_found_i);
    if (max_m_i > 0)
        fprintf(stderr, "Max min_n: %ld at m=%ld (p_m=%ld)\n",
                max_n_i, max_m_i, g_prime[max_m_i - 1]);

    fprintf(stderr, "\n=== SUMMARY Part (ii) ===\n");
    fprintf(stderr, "max n s.t. n^2 < pi(sieve): %ld\n", max_n_ii);
    fprintf(stderr, "Found:     %d / %d\n", found_ii, M_MAX);
    fprintf(stderr, "Not found (within bound): %d\n", not_found_ii);

    fprintf(stderr, "\nTotal elapsed: %.0f s\n", difftime(time(NULL), t_start));

    /* ---- Heuristic obstruction analysis: NOT FOUND cases ---- */
    fprintf(stderr, "\n=== Heuristic partial obstruction analysis (NOT FOUND, Part i) ===\n");
    fprintf(stderr, "NOTE: heuristic — confirmed obstructions only.\n");
    fprintf(stderr, "%-8s %-12s %-12s\n", "m", "p_m", "blocked%");
    for (int m = 1; m <= M_MAX; m++) {
        if (res_i[m] >= 0) continue;
        int64_t pm = g_prime[m - 1];
        int blocked = 0, check_range = 10000;
        for (int n = 1; n <= check_range; n++)
            if (has_partial_obstruction((int64_t)m + n, pm)) blocked++;
        fprintf(stderr, "%-8d %-12ld %.1f%%\n",
                m, pm, 100.0 * blocked / check_range);
    }

    free(res_i); free(res_ii); free(g_prime); free(g_bsieve);
    return 0;
}
