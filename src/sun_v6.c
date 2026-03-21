/*
 * sun_v6.c  --  Computational verification of Sun's Conjecture 4.1(i)
 *
 * Conjecture (Z.-W. Sun, arXiv:1309.1679, 2014-09-29):
 *   For every m in Z+, there exists n in Z+ such that
 *   (m + n) | (p_m^2 + p_n^2),
 *   where p_k denotes the k-th prime.
 *
 * OEIS: A247975.  a(m) = least such n.
 *
 * Changes from sun_v1.c:
 *
 *   [V2-01] __int128 arithmetic in find_min_n.
 *     pm^2 and pnr^2 are now computed via __int128, then reduced mod s.
 *     This is safe for any BOUND: with BOUND > ~3e9, pnr can exceed
 *     ~3e9 and pnr^2 would overflow int64_t (silent UB).  __int128
 *     costs one MULQ instruction (3 cycles on Zen 4) and eliminates
 *     the overflow risk for all future use of this code.
 *
 *   [V2-02] Off-by-one fix in find_min_n.
 *     lim was set to g_prime_count-1 when BOUND >= g_prime_count,
 *     causing the last prime in the table to be silently skipped.
 *     Corrected to g_prime_count (accesses g_prime[g_prime_count-1],
 *     the last valid index).
 *
 *   [V2-03] Additive remainder filter in modular_obstruction_filter.
 *     The filter previously computed s % q and s % q^2 via hardware
 *     division (~26 cycles each) on every candidate s.  Since s
 *     increments by 1 at each step, the remainders can be maintained
 *     by simple increment-and-reset (1-2 cycles), eliminating all
 *     division instructions from the inner loop.
 *     The remainder state (rem_q, rem_q2) is initialised once per m
 *     call and passed into the filter.
 *     Correctness: identical to the division-based version; the
 *     remainder update is exact because s increases by exactly 1.
 *
 *   [V2-04] strtoll with range validation replaces atoll/atoi.
 *     atoll gives undefined behaviour on overflow; strtoll sets errno.
 *
 *   [V2-05] Sieve limit estimate uses factor 1.3 instead of 2.0.
 *     estimate_pn(BOUND) already has a 5% margin built in; factor 2.0
 *     was excessive (doubled the sieve RAM unnecessarily).  Factor 1.3
 *     provides a comfortable margin while reducing memory use.
 *
 *   [V2-06] Prime table allocated with posix_memalign (64-byte aligned).
 *     Ensures cache-line alignment for the hot array.
 *
 *   [V6-01] Removed orphaned BOUND_INT64_SAFE define.
 *     Superseded by V4-01 (__int128 + divq_mod).
 *   [V6-02] Updated stale sun_v2 references to sun_v6.
 *   [V6-03] Integer arithmetic for prime-table slot count.
 *   [V6-04] Explicit casts to silence -Wconversion/-Wsign-conversion.
 *   [V6-05] Explicit (double) cast in estimate_pn.
 *   [V6-06] Documented volatile rationale in divq_mod asm.
 *
 * Compile:
 *   gcc -O3 -march=znver4 -mtune=znver4 -funroll-loops -ffast-math \
 *       -fomit-frame-pointer -flto -fopenmp \
 *       sun_v6.c -o sun_v6 -lm
 *
 * Usage:
 *   ./sun_v6 [M_MAX [BOUND [SIEVE_LIMIT]]]
 *
 * Output:
 *   stdout : CSV (one row per m)
 *   stderr : progress, timing, summary, unresolved list
 *
 * Author:  Carlo Corti
 * Date:    2026
 * License: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <sys/mman.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ------------------------------------------------------------------ */
/*  Configuration defaults                                             */
/* ------------------------------------------------------------------ */
#define DEFAULT_M_MAX   120000
#define DEFAULT_BOUND   2000000000LL

/* [V3-07] Cap to prevent estimate_pn overflow for extreme BOUND */
#define BOUND_ESTIMATE_CAP  2000000000000LL  /* 2e12 */

/* ------------------------------------------------------------------ */
/*  Global state                                                       */
/* ------------------------------------------------------------------ */
static uint8_t  *g_bsieve      = NULL;
static int64_t  *g_prime       = NULL;
static int64_t   g_prime_count = 0;
static int64_t   g_sieve_limit = 0;

/* ------------------------------------------------------------------ */
/*  Bitset sieve (one bit per odd integer >= 3)                        */
/*  Mapping: odd k = 2*i + 3  <-->  bit i  (i = (k-3)/2)             */
/*  Bit value: 0 = prime, 1 = composite.                              */
/* ------------------------------------------------------------------ */
#define BSIEVE_IDX(k)    (((k) - 3) >> 1)
#define BSIEVE_BYTES(l)  ((((l) - 3) >> 1) / 8 + 2)

static inline int bsieve_is_prime(int64_t k)
{
    if (k < 2)       return 0;
    if (k == 2)      return 1;
    if (k % 2 == 0)  return 0;
    int64_t i = BSIEVE_IDX(k);
    return !((g_bsieve[i >> 3] >> (i & 7)) & 1);
}

static inline void bsieve_set_composite(int64_t k)
{
    int64_t i = BSIEVE_IDX(k);
    g_bsieve[i >> 3] |= (uint8_t)(1 << (i & 7));
}

static void build_sieve(int64_t limit)
{
    size_t bsz = (size_t)BSIEVE_BYTES(limit);
    /* [V5-02] 64-byte aligned allocation for cache-line consistency */
    if (posix_memalign((void **)&g_bsieve, 64, bsz) != 0) {
        fprintf(stderr, "FATAL: bitset allocation failed (%.2f GB)\n",
                (double)bsz / 1e9);
        exit(1);
    }
    memset(g_bsieve, 0, bsz);
    fprintf(stderr, "  Bitset sieve:  %.2f GB\n", (double)bsz / 1e9);
    /* [V4-03] Huge pages for the sieve bitset */
#ifdef MADV_HUGEPAGE
    madvise(g_bsieve, bsz, MADV_HUGEPAGE);
#endif

    for (int64_t i = 3; i * i <= limit; i += 2)
        if (bsieve_is_prime(i))
            for (int64_t j = i * i; j <= limit; j += 2 * i)
                bsieve_set_composite(j);

    /* Pass 1: count primes */
    int64_t count = 1;  /* prime 2 */
    for (int64_t k = 3; k <= limit; k += 2)
        if (bsieve_is_prime(k)) count++;

    /* [V2-06] 64-byte aligned allocation for the prime table */
    /* [V6-03] Integer arithmetic avoids int64_t→double conversion */
    size_t slots = (size_t)(count + count / 100 + 64);
    if (posix_memalign((void **)&g_prime, 64,
                       sizeof(int64_t) * slots) != 0) {
        fprintf(stderr, "FATAL: prime table allocation failed (%.2f GB)\n",
                (double)(slots * sizeof(int64_t)) / 1e9);
        exit(1);
    }
    fprintf(stderr, "  Prime table:   %ld primes, %.2f GB\n",
            count, (double)(slots * sizeof(int64_t)) / 1e9);

    /* [V3-04] Request huge pages to reduce TLB misses on large tables */
#ifdef MADV_HUGEPAGE
    madvise(g_prime, sizeof(int64_t) * slots, MADV_HUGEPAGE);
#endif

    /* Pass 2: fill prime table */
    g_prime[g_prime_count++] = 2;
    for (int64_t k = 3; k <= limit; k += 2)
        if (bsieve_is_prime(k))
            g_prime[g_prime_count++] = k;
}

/* ------------------------------------------------------------------ */
/*  [V3-01] Inline asm: 128-bit dividend divided by 64-bit divisor.   */
/*  Uses hardware DIVQ instead of __umodti3 software routine.          */
/*  Precondition: mathematical quotient (hi:lo) / m fits in 64 bits.  */
/* ------------------------------------------------------------------ */
static inline int64_t divq_mod(uint64_t hi, uint64_t lo, uint64_t m)
{
    /* [V5-03] DEBUG: verify quotient fits in 64 bits (hi < m required).
     * Guaranteed by problem structure: pn^2/s ≈ n*log^2(n) ≪ 2^64.
     * This assert catches future code changes that might violate it. */
#ifdef DEBUG
    if (hi >= m) {
        fprintf(stderr,
            "FATAL divq_mod: hi=%lu >= m=%lu (quotient overflow)\n",
            (unsigned long)hi, (unsigned long)m);
        abort();
    }
#endif
    uint64_t rem;
    /* [V6-06] volatile is technically redundant (no side effects beyond
     * declared outputs), but retained as a safety net to prevent the
     * compiler from reordering or eliminating this asm in future
     * refactors where the same divq_mod call might appear twice. */
    __asm__ volatile (
        "divq %2"
        : "=d" (rem), "+a" (lo)
        : "rm" (m), "d" (hi)
        : "cc"
    );
    return (int64_t)rem;
}

/* ------------------------------------------------------------------ */
/*  Modular obstruction filter                                         */
/*                                                                     */
/*  Returns 1 if (m+n) | (p_m^2 + p_n^2) is provably impossible      */
/*  for ALL n with m+n = s.  Returns 0 if s cannot be ruled out.      */
/*                                                                     */
/*  [V2-03] Additive remainder maintenance.                            */
/*  Instead of computing s % q via hardware division at each step,    */
/*  the caller maintains rem_q[i] and rem_q2[i]: the remainders of    */
/*  the current s modulo FILTER_Q[i] and FILTER_Q2[i] respectively.   */
/*  These are updated by the caller with a simple increment-and-reset. */
/* ------------------------------------------------------------------ */

#define N_FILTER_PRIMES 4
static const int64_t FILTER_Q[N_FILTER_PRIMES]  = {3, 7, 11, 19};
static const int64_t FILTER_Q2[N_FILTER_PRIMES] = {9, 49, 121, 361};

/* Returns 1 if v_q(s) is odd given precomputed remainders.
   rem_q  = s % q   (maintained additively by caller)
   rem_q2 = s % q^2 (maintained additively by caller)
   For v_q >= 2 (rare: rem_q==0 && rem_q2==0), falls back to division. */
static inline int vq_is_odd_fast(int64_t s, int qi,
                                  int rem_q, int rem_q2)
{
    if (rem_q  != 0) return 0;   /* q does not divide s */
    if (rem_q2 != 0) return 1;   /* v_q = 1 (odd) — most common case */
    /* v_q >= 2: count fully (rare, cost amortised) */
    int64_t t = s;
    int e = 0;
    int64_t q = FILTER_Q[qi];
    while (t % q == 0) { t /= q; e++; }
    return e & 1;
}

/*
 * Precompute bitmask: bit i set iff p_m == FILTER_Q[i].
 * For m > 8 (p_m > 19) always returns 0.  Call once per m.
 */
static inline uint8_t compute_pm_qmask(int64_t pm)
{
    uint8_t mask = 0;
    for (int i = 0; i < N_FILTER_PRIMES; i++)
        if (pm == FILTER_Q[i]) mask |= (uint8_t)(1 << i);
    return mask;
}

/*
 * rem[2*i]   = current s % FILTER_Q[i]
 * rem[2*i+1] = current s % FILTER_Q2[i]
 * Caller updates these after each increment of s.
 */
static inline int modular_obstruction_filter(int64_t s, int64_t pm,
                                              uint8_t pm_qmask,
                                              const int * restrict rem)
{
    for (int i = 0; i < N_FILTER_PRIMES; i++) {
        if (pm_qmask & (1 << i)) continue;
        if (vq_is_odd_fast(s, i, rem[2*i], rem[2*i+1]))
            return 1;
    }
    /* Retained: prime s ≡ 3 (mod 4) with s ≠ p_m */
    if (s <= g_sieve_limit && bsieve_is_prime(s)
            && (s % 4 == 3) && (s != pm))
        return 1;
    return 0;
}

/* ------------------------------------------------------------------ */
/*  Core search: find a(m) = least n in [1, bound] s.t.               */
/*  (m+n) | (p_m^2 + p_n^2).  Returns -1 if not found.               */
/*                                                                     */
/*  [V2-01] __int128 arithmetic.                                       */
/*  pm_sq = pm^2 computed once as __int128, then pm_r = pm_sq % s     */
/*  recomputed at each s (cheap: just a 128÷64 division, 1 DIVQ).     */
/*  pnr^2 computed as __int128 then reduced mod s.                    */
/*  This is safe for any BOUND without overflow risk.                  */
/*                                                                     */
/*  [V2-02] lim = g_prime_count (not g_prime_count-1).                */
/*  [V2-03] Additive remainder state rem[8] initialised here.         */
/* ------------------------------------------------------------------ */
static inline int64_t find_min_n(int m, int64_t bound)
{
    /* [V3-03] a(1) = 1 is known: 2 | p_1^2 + p_1^2 = 8 */
    if (m == 1) return 1;

    int64_t pm = g_prime[m - 1];

    /* [V3-01] Decompose pm^2 into hi:lo for DIVQ
     * [V6-04] Explicit uint64_t casts silence -Wsign-conversion;
     *         redundant mask on pm_lo removed (uint64_t cast truncates). */
    unsigned __int128 pm_sq128 = (unsigned __int128)(uint64_t)pm
                                 * (uint64_t)pm;
    uint64_t pm_hi = (uint64_t)(pm_sq128 >> 64);
    uint64_t pm_lo = (uint64_t)pm_sq128;

    int64_t lim      = (bound <= g_prime_count) ? bound : g_prime_count;
    uint8_t pm_qmask = compute_pm_qmask(pm);

    /* [V2-03] Initialise additive remainder state for s = m+1 */
    int64_t s0 = (int64_t)m + 1;
    int rem[2 * N_FILTER_PRIMES];
    for (int i = 0; i < N_FILTER_PRIMES; i++) {
        rem[2*i]   = (int)(s0 % FILTER_Q[i]);
        rem[2*i+1] = (int)(s0 % FILTER_Q2[i]);
    }

    for (int64_t s = s0; s <= (int64_t)m + lim; s++) {

        /* [V3-05] hint: filter rejects ~49%, so "pass" is slightly
           more common; false is the "rarer" branch for the outer if */
        if (__builtin_expect(
                !modular_obstruction_filter(s, pm, pm_qmask, rem), 1)) {
            int64_t n   = s - m;
            int64_t pn  = g_prime[n - 1];

            /* [V3-01] pm^2 % s via hardware DIVQ */
            int64_t r   = divq_mod(pm_hi, pm_lo, (uint64_t)s);

            /* [V4-01] pn^2 % s via __int128 + divq_mod.
             * Eliminates one hardware division (previously pn%s then pnr*pnr%s).
             * Also correct for any BOUND (pn can reach ~4.7e10 at BOUND=2e9,
             * so pn^2 ~ 2.2e21 overflows int64_t; __int128 handles this). */
            unsigned __int128 pn_sq128 = (unsigned __int128)pn * (unsigned __int128)pn;
            uint64_t pn_hi = (uint64_t)(pn_sq128 >> 64);
            uint64_t pn_lo = (uint64_t)(pn_sq128);
            int64_t pnr2   = (int64_t)divq_mod(pn_hi, pn_lo, (uint64_t)s);

            /* [V3-03] direct comparison replaces (r + pnr2) % s == 0 */
            if (r + pnr2 == s)
                return n;
        }

        /* [V2-03] Update remainders for s+1 */
        for (int i = 0; i < N_FILTER_PRIMES; i++) {
            if (++rem[2*i]   == (int)FILTER_Q[i])  rem[2*i]   = 0;
            if (++rem[2*i+1] == (int)FILTER_Q2[i]) rem[2*i+1] = 0;
        }
    }
    return -1;
}

/* ------------------------------------------------------------------ */
/*  Estimate p_n ~ n*(ln n + ln ln n - 1) with safety margin          */
/* ------------------------------------------------------------------ */
static int64_t estimate_pn(int64_t n)
{
    /* [V3-07] Cap to prevent int64_t overflow for extreme BOUND */
    if (n > BOUND_ESTIMATE_CAP) {
        fprintf(stderr,
            "WARNING: BOUND exceeds estimate cap %ld; "
            "provide SIEVE_LIMIT explicitly.\n",
            (long)BOUND_ESTIMATE_CAP);
        n = BOUND_ESTIMATE_CAP;
    }
    double ln_n  = log((double)n);
    double ln2_n = log(ln_n);
    /* [V6-05] Explicit (double)n documents intentional conversion */
    return (int64_t)((double)n * (ln_n + ln2_n - 1.0) * 1.05);
}

/* ------------------------------------------------------------------ */
/*  [V2-04] Safe argument parsing with strtoll                         */
/* ------------------------------------------------------------------ */
static int64_t parse_int64(const char *s, const char *name,
                           int64_t lo, int64_t hi)
{
    char *end;
    errno = 0;
    long long v = strtoll(s, &end, 10);
    if (errno || *end != '\0' || v < lo || v > hi) {
        fprintf(stderr, "FATAL: invalid value for %s: '%s' "
                "(must be %lld..%lld)\n", name, s, (long long)lo,
                (long long)hi);
        exit(1);
    }
    return (int64_t)v;
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int main(int argc, char *argv[])
{
    int     M_MAX = DEFAULT_M_MAX;
    int64_t BOUND = DEFAULT_BOUND;

    /* [V2-04] validated parsing */
    if (argc > 1) M_MAX         = (int)parse_int64(argv[1], "M_MAX",
                                                    1, 10000000);
    if (argc > 2) BOUND         = parse_int64(argv[2], "BOUND",
                                              1, (int64_t)4e18);
    if (argc > 3) g_sieve_limit = parse_int64(argv[3], "SIEVE_LIMIT",
                                              1, (int64_t)4e18);

    /* [V4-01] pnr^2 now uses __int128 + divq_mod → no BOUND limit for pnr^2 */

    /* [V2-05] Auto sieve limit: factor 1.3 instead of 2.0 */
    if (g_sieve_limit == 0)
        g_sieve_limit = estimate_pn(BOUND) * 13 / 10;

    /* ---- Header ---- */
    time_t t_start = time(NULL);
    fprintf(stderr, "sun_v6  --  Sun Conjecture 4.1(i), sequence A247975\n");
    fprintf(stderr, "  M_MAX       : %d\n", M_MAX);
    fprintf(stderr, "  BOUND (n)   : %ld\n", (long)BOUND);
    fprintf(stderr, "  SIEVE_LIMIT : %ld\n", (long)g_sieve_limit);
#ifdef _OPENMP
    fprintf(stderr, "  OpenMP      : %d threads\n", omp_get_max_threads());
#else
    fprintf(stderr, "  OpenMP      : disabled\n");
#endif
    fprintf(stderr, "\n");

    /* ---- Build sieve and prime table ---- */
    fprintf(stderr, "Building sieve up to %ld ...\n", (long)g_sieve_limit);
    build_sieve(g_sieve_limit);
    fprintf(stderr, "  Largest prime in table: %ld\n",
            g_prime[g_prime_count - 1]);
    if (g_prime_count >= 100000)
        fprintf(stderr,
                "  Spot check: p(1000)=%ld  p(10000)=%ld  p(100000)=%ld\n",
                g_prime[999], g_prime[9999], g_prime[99999]);
    fprintf(stderr, "  Sieve time: %.1f s\n\n",
            difftime(time(NULL), t_start));

    /* Clamp BOUND if sieve is too small */
    if (BOUND > g_prime_count) {
        fprintf(stderr,
            "WARNING: BOUND=%ld > pi(SIEVE_LIMIT)=%ld. "
            "Clamping to %ld.\n"
            "  Rerun with SIEVE_LIMIT >= %.3e to search n up to %ld.\n\n",
            (long)BOUND, (long)g_prime_count, (long)g_prime_count,
            (double)BOUND * log((double)BOUND), (long)BOUND);
        BOUND = g_prime_count;
    }

    /* ---- Allocate result array ---- */
    /* [V6-04] Explicit size_t cast silences -Wsign-conversion */
    int64_t *result = (int64_t *)calloc((size_t)(M_MAX + 1), sizeof(int64_t));
    if (!result) { fprintf(stderr, "FATAL: out of memory\n"); exit(1); }

    /* ---- Parallel search ---- */
    fprintf(stderr, "Searching m = 1..%d, n-bound = %ld ...\n",
            M_MAX, (long)BOUND);

    int found     = 0;
    int not_found = 0;

    #pragma omp parallel for schedule(dynamic, 10) \
            reduction(+:found, not_found)
    for (int m = 1; m <= M_MAX; m++) {
        if (m - 1 >= g_prime_count) continue;
        int64_t n = find_min_n(m, BOUND);
        result[m] = n;
        if (n >= 0) found++;
        else        not_found++;
    }

    /* ---- CSV output ---- */
    printf("# sun_v6 output: Sun Conjecture 4.1(i), OEIS A247975\n");
    printf("# M_MAX=%d  BOUND=%ld  SIEVE_LIMIT=%ld\n",
           M_MAX, (long)BOUND, (long)g_sieve_limit);
    printf("# %s", ctime(&t_start));
    printf("m,p_m,a_m,p_am,s,cramér_ratio\n");

    for (int m = 1; m <= M_MAX; m++) {
        if (m - 1 >= g_prime_count) break;
        int64_t pm = g_prime[m - 1];
        int64_t n  = result[m];
        if (n >= 0) {
            int64_t pn = g_prime[n - 1];
            int64_t s  = (int64_t)m + n;
            /* [V3-06] avoid div-by-zero: log(1) = 0 */
            double  cr = (m > 1)
                ? (double)n / ((double)m * log((double)m))
                : 0.0;
            printf("%d,%ld,%ld,%ld,%ld,%.4f\n", m, pm, n, pn, s, cr);
        } else {
            printf("%d,%ld,UNRESOLVED,-1,-1,-1\n", m, pm);
        }
    }

    /* ---- Summary ---- */
    double elapsed = difftime(time(NULL), t_start);
    fprintf(stderr, "\n=== SUMMARY ===\n");
    fprintf(stderr, "m range    : 1..%d\n",   M_MAX);
    fprintf(stderr, "n-bound    : %ld\n",      (long)BOUND);
    fprintf(stderr, "p_{BOUND}  : %ld\n",      g_prime[BOUND - 1]);
    fprintf(stderr, "Resolved   : %d / %d (%.4f%%)\n",
            found, M_MAX, 100.0 * found / M_MAX);
    fprintf(stderr, "Unresolved : %d\n",       not_found);
    fprintf(stderr, "Elapsed    : %.0f s (%.1f min)\n",
            elapsed, elapsed / 60.0);

    /* ---- Unresolved list ---- */
    if (not_found > 0) {
        fprintf(stderr, "\n=== UNRESOLVED CASES ===\n");
        fprintf(stderr, "%-8s  %-14s  %-10s\n", "m", "p_m", "C(m) >");
        for (int m = 1; m <= M_MAX; m++) {
            if (result[m] >= 0) continue;
            int64_t pm   = g_prime[m - 1];
            double  cr_lb = (double)BOUND
                            / ((double)m * log((double)m));
            fprintf(stderr, "%-8d  %-14ld  %-10.0f\n", m, pm, cr_lb);
        }
    }

    free(result);
    free(g_prime);
    free(g_bsieve);
    return 0;
}
