/* =============================================================================
 * verify_min_ext.c
 * Independent minimality verifier for Sun's Conjecture 4.1(i) / OEIS A247975
 * Extended version — covers Phase 1 AND Extended search (n up to 2×10^11).
 *
 * Condition verified: (m + n) | (p_m^2 + p_n^2)
 *
 * Usage:
 *   ./verify_min_ext <m> <p_m> <a_m> <p_a_m>
 *
 *   m     = index (positive integer)
 *   p_m   = p_m, the m-th prime (for cross-check)
 *   a_m   = claimed answer a(m) to verify
 *   p_a_m = p_{a(m)}, the a(m)-th prime (used as sieve-limit hint; an
 *            analytical upper bound is also computed as a safety margin)
 *
 * Output (one line, parseable by the Python orchestrator):
 *   VERIFIED   m first_n p_first elapsed_s
 *   MISMATCH   m first_n p_first elapsed_s   (found earlier solution)
 *   NOT_FOUND  m elapsed_s                    (no solution in [1, a_m])
 *
 * Algorithm:
 *   Fully independent segmented sieve of Eratosthenes (no code shared
 *   with the original Phase-6 C program).  Primes are generated in order;
 *   for each prime p_n at index n we test the divisibility condition.
 *   The first n that satisfies it is reported and compared with a_m.
 *
 * Arithmetic safety (Extended range, n up to 2×10^11):
 *   - pm  <= p_{120000} = 1,567,271     ->  pm^2  < 2^42      (uint64_t OK)
 *   - pn  up to p_{2e11} ~ 5.2e12       ->  pn^2  ~ 2.7e25    (needs __uint128_t)
 *   - s   = m + n < 120000 + 2e11 ~ 2e11 < 2^38               (uint64_t OK)
 *   - pm^2 + pn^2 < 2.46e12 + 2.7e25 ~ 2^85                   (safe in __uint128_t)
 *
 * Version history:
 *   v1 — initial Extended version
 *   v2 — fix: `(__uint128_t)p * (uint64_t)p` replaced by
 *        `(__uint128_t)p * (__uint128_t)p` to guarantee 128-bit multiply
 *        with no ambiguity in C promotion rules (both occurrences: p=2
 *        special case and main sieve loop).
 *
 * Changes from verify_min.c (Phase 1 only):
 *   [EXT-1] Divisibility check uses __uint128_t to avoid overflow for
 *           large pn (pn^2 overflows int64_t when pn > ~3e9).
 *           Both operands of the multiply are cast to __uint128_t.
 *   [EXT-2] MAX_SMALLS raised from 25,000 to 200,000.
 *           Rationale: sqrt(p_{2e11}) ~ sqrt(5.2e12) ~ 2.28e6;
 *           pi(2.28e6) ~ 166,000 primes needed as sieve factors.
 *   [EXT-3] sml[] widened from int32_t to int64_t to hold primes > 2^31
 *           (largest small prime can reach ~2.28e6; still fits int32_t,
 *            but int64_t is used for uniformity and future safety).
 *   [EXT-4] build_small() allocation uses int64_t for sq to avoid
 *           truncation when sq > INT_MAX (not expected here, but safe).
 *   [EXT-5] Header comment updated to reflect Extended arithmetic bounds.
 *
 * This program is a strict superset of verify_min.c: it produces identical
 * results for all Phase 1 cases (a_m < 2×10^9) and correct results for
 * Extended cases (a_m up to 2×10^11).
 *
 * Compile:
 *   gcc -O2 -o verify_min_ext verify_min_ext.c -lm
 *
 * Smoke test (absolute reference value):
 *   ./verify_min_ext 11924 127289 6119832581 151142662267
 *   Expected: VERIFIED 11924 6119832581 151142662267 <elapsed>
 *
 * Author:  generated for Carlo Corti's A247975 verification campaign
 * =============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

/* ---------- Sieve configuration ---------- */

/* One byte per odd number in the segment.
 * SEG_BYTES = 2^22 = 4,194,304  ->  segment covers 8,388,608 numbers.
 * Fits comfortably in L3 cache. */
#define SEG_BYTES  (1 << 22)

/* [EXT-2] Max number of small primes needed for Extended range.
 * sqrt(p_{2e11}) ~ sqrt(5.2e12) ~ 2,280,351
 * pi(2,280,351) ~ 166,000 primes needed as sieve factors.
 * 200,000 is a safe upper bound with margin. */
#define MAX_SMALLS 200000

/* ---------- Globals ---------- */
static unsigned char  seg[SEG_BYTES];        /* segment: 0=prime cand, 1=composite */
/* [EXT-3] sml[] widened to int64_t for uniformity (values still fit int32_t) */
static int64_t        sml[MAX_SMALLS];       /* small primes (odd only) */
static int64_t        sml_off[MAX_SMALLS];   /* next odd multiple per small prime */
static int            sml_cnt = 0;

/* ---------- Build small primes up to `sq` via basic sieve ---------- */
static void build_small(int64_t sq) {
    /* [EXT-4] sq is int64_t; allocation size must not overflow size_t.
     * For Extended: sq ~ 2.28e6, allocation ~ 2.28 MB — well within limits. */
    size_t n = (size_t)(sq + 1);
    unsigned char *s = (unsigned char *)calloc(n, 1);
    if (!s) { fprintf(stderr, "OOM in build_small(%lld)\n", (long long)sq); exit(1); }
    s[0] = s[1] = 1;
    for (int64_t i = 2; i * i <= sq; i++)
        if (!s[i])
            for (int64_t j = i * i; j < (int64_t)n; j += i)
                s[j] = 1;
    sml_cnt = 0;
    for (int64_t i = 3; i < (int64_t)n && sml_cnt < MAX_SMALLS; i += 2)
        if (!s[i]) sml[sml_cnt++] = i;
    free(s);
    /* Initialise each small prime's offset to p^2 (standard sieve start).
     * All p >= 3  ->  p^2 >= 9 >= 3 = first seg_lo, so no underflow. */
    for (int i = 0; i < sml_cnt; i++)
        sml_off[i] = sml[i] * sml[i];
}

/* ---------- Compute a safe upper bound for p_N ----------
 * Uses Robin's bound: p_n < n*(ln n + ln ln n + 2.0) for n >= 6.
 * Valid and generous for all n up to at least 2×10^11. */
static int64_t prime_upper_bound(int64_t N) {
    if (N < 6) return 20;
    double lnN = log((double)N);
    double bound = (double)N * (lnN + log(lnN) + 2.0);
    return (int64_t)bound + 1000000LL;
}

/* ---------- Main ---------- */
int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr,
            "Usage: %s <m> <p_m> <a_m> <p_a_m>\n"
            "  Verifies that a(m) is the minimum n with (m+n)|(p_m^2+p_n^2).\n"
            "  Covers both Phase 1 (a_m < 2e9) and Extended (a_m up to 2e11).\n",
            argv[0]);
        return 2;
    }

    int64_t m   = atoll(argv[1]);
    int64_t pm  = atoll(argv[2]);
    int64_t a_m = atoll(argv[3]);
    int64_t pam = atoll(argv[4]);   /* p_{a(m)} from CSV, used as sieve-limit hint */

    /* Sieve limit: take the maximum of:
     *   (a) pam + 1,000,000  (exact hint from the CSV)
     *   (b) analytical upper bound (safety margin if pam is slightly off) */
    int64_t sieve_lim = pam + 1000000LL;
    int64_t analytical = prime_upper_bound(a_m);
    if (analytical > sieve_lim) sieve_lim = analytical;

    fprintf(stderr,
        "[verify_min_ext] m=%lld  p_m=%lld  a_m=%lld  p_a_m=%lld  sieve_lim=%lld\n",
        (long long)m, (long long)pm, (long long)a_m,
        (long long)pam, (long long)sieve_lim);
    fflush(stderr);

    /* Build small primes up to ceil(sqrt(sieve_lim)) */
    int64_t sq = (int64_t)sqrt((double)sieve_lim) + 2;
    build_small(sq);
    fprintf(stderr, "[verify_min_ext] %d small primes up to sqrt(sieve_lim)~%lld\n",
            sml_cnt, (long long)sq);
    fflush(stderr);

    /* Precompute pm^2.
     * pm <= p_{120000} = 1,567,271 < 2^21  ->  pm^2 < 2^42  (uint64_t OK) */
    uint64_t pm2 = (uint64_t)pm * (uint64_t)pm;

    /* ---- Progress timing ---- */
    struct timespec ts0, ts_now;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    double last_report = 0.0;

    int64_t first_n = -1;   /* index of first solution found, -1=not yet */
    int64_t first_p = -1;   /* corresponding prime */
    int64_t n       = 0;    /* running prime count */

    /* ------ Check p=2 (n=1) separately ------ */
    {
        n = 1;
        /* [EXT-1] Use __uint128_t for the divisibility check.
         * s = m + n; we test: (pm^2 + pn^2) % s == 0
         * For p=2: pn^2 = 4, no overflow risk, but use the unified formula
         * for consistency. */
        uint64_t s = (uint64_t)(m + n);
        __uint128_t sum = (__uint128_t)pm2 + (__uint128_t)2 * (__uint128_t)2;
        if (sum % s == 0) {
            first_n = 1; first_p = 2;
        }
        if (n >= a_m && first_n < 0) {
            first_n = 0;   /* sentinel: searched all, none found */
        }
    }

    /* ------ Segmented sieve over odd numbers >= 3 ------ */
    int64_t seg_lo = 3;
    while (first_n < 0 && seg_lo <= sieve_lim) {
        /* Segment covers odds: seg_lo, seg_lo+2, ..., seg_lo + 2*(SEG_BYTES-1) */
        int64_t seg_hi = seg_lo + 2LL * (SEG_BYTES - 1);
        if (seg_hi > sieve_lim) seg_hi = sieve_lim;

        /* Clear: all candidates are prime until marked composite */
        memset(seg, 0, SEG_BYTES);

        /* Sieve: for each small prime p, mark its odd multiples in [seg_lo, seg_hi] */
        for (int i = 0; i < sml_cnt; i++) {
            int64_t p = sml[i];
            int64_t j = sml_off[i];

            /* Advance offset to first odd multiple of p >= seg_lo if needed */
            if (j < seg_lo) {
                int64_t r = seg_lo % p;
                j = (r == 0) ? seg_lo : seg_lo + (p - r);
                if (j % 2 == 0) j += p;   /* ensure odd (p is odd) */
            }

            /* Mark odd multiples of p within this segment */
            for (; j <= seg_hi; j += 2 * p)
                seg[(j - seg_lo) >> 1] = 1;

            sml_off[i] = j;   /* first unchecked odd multiple for the next segment */
        }

        /* Enumerate primes in this segment and check the divisibility condition */
        int64_t cnt = (seg_hi - seg_lo) / 2 + 1;   /* number of odds in segment */
        for (int64_t k = 0; k < cnt; k++) {
            if (!seg[k]) {
                int64_t p = seg_lo + 2 * k;
                n++;

                /* [EXT-1] Check: (m + n) | (pm^2 + pn^2) using __uint128_t.
                 *
                 * Bounds:
                 *   pm^2  < 2^42          (uint64_t, precomputed above)
                 *   pn    up to ~5.2e12   ->  pn^2 ~ 2.7e25 (overflows int64_t!)
                 *   s     = m + n < 2^38  (uint64_t OK)
                 *   sum   = pm^2 + pn^2 < 2^85  (safe in __uint128_t)
                 *
                 * Note: s < 2^38, so (uint64_t)s is exact for the % operator. */
                uint64_t s = (uint64_t)(m + n);
                __uint128_t sum = (__uint128_t)pm2 + (__uint128_t)p * (__uint128_t)p;
                if (sum % s == 0) {
                    first_n = n;
                    first_p = p;
                    goto done;
                }

                if (n >= a_m) {
                    first_n = 0;   /* searched all n in [1, a_m], none found */
                    goto done;
                }
            }
        }

        /* Periodic progress report to stderr (every 30 seconds) */
        clock_gettime(CLOCK_MONOTONIC, &ts_now);
        double elapsed = (ts_now.tv_sec - ts0.tv_sec)
                       + 1e-9 * (ts_now.tv_nsec - ts0.tv_nsec);
        if (elapsed - last_report >= 30.0) {
            double pct = 100.0 * (double)(seg_hi) / (double)sieve_lim;
            fprintf(stderr,
                "[verify_min_ext] m=%lld  n=%lld / %lld  seg_hi=%lld (%.1f%%)  %.0fs\n",
                (long long)m, (long long)n, (long long)a_m,
                (long long)seg_hi, pct, elapsed);
            fflush(stderr);
            last_report = elapsed;
        }

        seg_lo += 2LL * SEG_BYTES;
    }

done:;
    /* Final timing */
    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    double elapsed = (ts_now.tv_sec - ts0.tv_sec)
                   + 1e-9 * (ts_now.tv_nsec - ts0.tv_nsec);

    /* Print machine-readable result line */
    if (first_n == a_m) {
        printf("VERIFIED %lld %lld %lld %.2f\n",
               (long long)m, (long long)first_n,
               (long long)first_p, elapsed);
    } else if (first_n > 0 && first_n < a_m) {
        /* Found a SMALLER solution -> claimed a_m is not minimal */
        printf("MISMATCH %lld %lld %lld %.2f\n",
               (long long)m, (long long)first_n,
               (long long)first_p, elapsed);
    } else {
        /* first_n == 0: exhausted [1, a_m] with no solution */
        printf("NOT_FOUND %lld %.2f\n", (long long)m, elapsed);
    }

    fflush(stdout);
    return (first_n == a_m) ? 0 : 1;
}
