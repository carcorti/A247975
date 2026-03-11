/*
 * sun41_v12c.c  --  Targeted search for unresolved cases of Sun's Conjecture 4.1(i)
 *
 * Conjecture (Z.-W. Sun, 2014-09-29, arXiv:1309.1679):
 *   For every m in Z+, there exists n in Z+ such that (m+n) | (p_m^2 + p_n^2),
 *   where p_k denotes the k-th prime.
 *
 * OEIS reference: A247975.
 *
 * PURPOSE:
 *   This program searches exclusively for the 49 values of m <= 100000 that
 *   were unresolved in Phase 6 (sun41_v6b.c, bound n <= 2e9).  Instead of
 *   building a monolithic prime table in RAM (which required ~20 GB for the
 *   Phase 6 run), it uses a segmented sieve that generates primes on the fly
 *   in blocks, reducing peak RAM to O(block_size) ~ a few hundred MB.
 *   This makes bounds n <= 5e9 (and beyond) feasible on 32 GB machines.
 *
 * CHANGES v12c (fase 3: n = 4e9 .. 5e9):
 *   - PHASE6_BOUND updated to 4000000000 (fase-2 bound).
 *   - UNRESOLVED_M reduced from 38 to 34 cases (4 resolved in fase 2).
 *     Resolved in fase 2: 49889, 61007, 90161, 98008.
 *   - Binary renamed sun41_v12c.
 *
 * CHANGES v12b (fase 2: n = 3e9 .. 4e9):
 *   - PHASE6_BOUND updated to 3000000000 (fase-1 bound).
 *   - UNRESOLVED_M reduced from 49 to 37 cases (11 resolved in fase 1,
 *     plus m=37253 which was a duplicate entry — see note below).
 *     Resolved in fase 1: 18241, 20686, 36726, 67273, 68744, 71049,
 *                         84521, 92335, 93790, 93890, 94768.
 *   - Binary renamed sun41_v12b.
 *
 * CHANGES v12 (progress report polish):
 *   - REPORT_INTERVAL macro moved from inside search_case() to file scope,
 *     following standard C style for preprocessor constants.
 *   - Progress percentage clamped to [0, 100] to prevent ETA going negative
 *     if g_segs_done ever exceeds the (estimated) g_segs_total.
 *   - Version strings and binary name updated throughout.
 *
 * CHANGES v11 (targeted segmented sieve -- based on v6b):
 *   - Global prime table g_prime[] removed.  p_n values are generated
 *     on the fly via a segmented Sieve of Eratosthenes (block size
 *     configurable, default 20000000 odd slots ~ 2.5 MB per thread).
 *   - Part (ii) removed (not needed for the 49 unresolved cases).
 *   - M_MAX loop replaced by a hardcoded list of the 49 unresolved m-values.
 *   - p_m is looked up from a small auxiliary prime table covering only
 *     the range [2, SMALL_SIEVE_LIMIT=1400000], requiring < 1 MB.
 *   - The sieve base primes (for crossing off composites in each segment)
 *     need only go up to sqrt(p_{BOUND}); for BOUND=5e9, sqrt(p_n) ~ 3.5e5,
 *     which is well within the small prime table.
 *   - Seek strategy: to locate the segment containing p_{start_n}
 *     (the first prime index beyond the Phase-6 bound n=2e9), the Cipolla
 *     asymptotic formula estimates p_{start_n} with error < 0.01%
 *     (~5e6 for n=2e9).  We back off by 4 blocks (~1.6e8 integers) from
 *     the estimate, count primes exactly in that window, then run the
 *     main search from the correct offset.  Seek cost: ~4 block scans
 *     (< 1 s) vs ~1100 sequential blocks (~ 60 s) without the PNT skip.
 *   - Modular obstruction filter (v4) and pm_qmask optimisation (v5)
 *     retained unchanged.
 *   - OpenMP parallelisation: outer loop over the 49 cases, dynamic
 *     scheduling.  Each thread owns its own segment buffer.
 *
 * RAM budget (BOUND=5e9, BLOCK_SIZE=20000000, 16 threads):
 *   Small prime table (up to 1.4e6)           :  < 1 MB
 *   Base sieve bitset (up to 1.4e6)           :  < 1 MB
 *   Per-thread segment buffer (20M bits)      : ~2.5 MB x 16 = 40 MB
 *   Temp counting buffer (one per thread)     : ~2.5 MB x 16 = 40 MB
 *   Total                                     : < 100 MB
 *
 * Compile (Ryzen 9 7940HS, GCC):
 *   gcc -O3 -march=znver4 -mtune=znver4 -funroll-loops -ffast-math \
 *       -fomit-frame-pointer -flto -fopenmp \
 *       sun41_v12c.c -o sun41_v12c -lm
 *
 * Usage:
 *   ./sun41_v12c [BOUND [BLOCK_SIZE]]
 *   Defaults: BOUND=5000000000  BLOCK_SIZE=20000000
 *
 *   BOUND      : upper limit on n-index (search n = 1 .. BOUND)
 *   BLOCK_SIZE : number of odd-integer slots per sieve segment
 *
 * Output:
 *   stdout : CSV  m,p_m,min_n,status
 *   stderr : progress and summary
 *
 * Note on resolved output: the CSV records the n-index of the solution.
 *   To recover p_{min_n}, run any prime generator for that index, or use
 *   the estimate p_n ~ n*(ln(n)+ln(ln(n))-1) as a starting point.
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

#ifdef _OPENMP
#include <omp.h>
#endif

/* ================================================================== */
/*  Configuration defaults                                             */
/* ================================================================== */
#define DEFAULT_BOUND       5000000000LL
#define DEFAULT_BLOCK_SIZE    20000000LL   /* odd-integer slots per segment */

/* Small prime table covers all p_m (m up to 99105, p_m up to 1286953) and
 * all sieve base primes up to sqrt(p_{5e9}) ~ 3.5e5.                       */
#define SMALL_SIEVE_LIMIT   1400000L

/* Phase-2 search bound (n-index). v12c continues from here. */
#define PHASE6_BOUND  4000000000LL

/* ================================================================== */
/*  The 34 unresolved m-values after Phase 2 (n-bound 4e9)           */
/*  Removed (resolved in fase 2): 49889, 61007, 90161, 98008.        */
/* ================================================================== */
static const int UNRESOLVED_M[] = {
    11924, 19234, 22186, 22802, 23995, 33566, 36303,
    37249, 37253, 37949, 43297, 48834, 51816, 53963,
    55468, 55874, 66257, 67203, 67447, 69546, 72568,
    76761, 76868, 81986, 82096, 82300, 86557, 90005,
    94429, 97251, 97351, 98258, 98379, 99105
};
#define N_UNRESOLVED ((int)(sizeof(UNRESOLVED_M)/sizeof(UNRESOLVED_M[0])))

/* ================================================================== */
/*  Small base prime table                                             */
/* ================================================================== */
static int64_t *g_small_prime = NULL;
static int64_t  g_small_count = 0;
static uint8_t *g_base_bsieve = NULL;

#define BS_IDX(k)    (((k)-3)>>1)
#define BS_BYTES(l)  ((((l)-3)>>1)/8+2)

static inline int bs_is_prime(long k)
{
    if (k < 2)    return 0;
    if (k == 2)   return 1;
    if (k%2 == 0) return 0;
    long i = BS_IDX(k);
    return !((g_base_bsieve[i>>3] >> (i&7)) & 1);
}

static void build_small_sieve(long limit)
{
    size_t bsz = BS_BYTES(limit);
    g_base_bsieve = (uint8_t *)calloc(bsz, 1);
    if (!g_base_bsieve) { fprintf(stderr,"OOM: base bsieve\n"); exit(1); }

    for (long i = 3; i*i <= limit; i += 2)
        if (bs_is_prime(i))
            for (long j = i*i; j <= limit; j += 2*i) {
                long idx = BS_IDX(j);
                g_base_bsieve[idx>>3] |= (uint8_t)(1<<(idx&7));
            }

    int64_t cnt = 1;
    for (long k = 3; k <= limit; k += 2)
        if (bs_is_prime(k)) cnt++;

    g_small_prime = (int64_t *)malloc(sizeof(int64_t)*((size_t)cnt+64));
    if (!g_small_prime) { fprintf(stderr,"OOM: small prime table\n"); exit(1); }

    g_small_prime[g_small_count++] = 2;
    for (long k = 3; k <= limit; k += 2)
        if (bs_is_prime(k))
            g_small_prime[g_small_count++] = k;

    fprintf(stderr,"  Small prime table: %ld primes up to %ld (%.2f MB)\n",
            g_small_count, limit,
            (double)(bsz + (size_t)cnt*sizeof(int64_t))/1e6);
}

/* ================================================================== */
/*  Cipolla estimate for p_n (1-based)                                */
/* ================================================================== */
static int64_t pnt_pn(int64_t n)
{
    if (n <= 1) return 2;
    double dn  = (double)n;
    double ln  = log(dn);
    double lln = log(ln);
    double est = dn * (ln + lln - 1.0
                       + (lln - 2.0)/ln
                       - (lln*lln - 6.0*lln + 11.0)/(2.0*ln*ln));
    if (est < 2.0) est = 2.0;
    int64_t r = (int64_t)est;
    if (r % 2 == 0) r--;
    if (r < 3) r = 3;
    return r;
}

/* ================================================================== */
/*  Segmented sieve: per-thread context                               */
/* ================================================================== */
typedef struct {
    uint8_t  *buf;
    int64_t   seg_low;
    int64_t   block_size;
    int64_t  *bp_val;
    int64_t  *bp_off;
    int64_t   n_bp;
} seg_ctx_t;

/*
 * Create a segment context.  Segments cover odd numbers only.
 * start_odd: first odd number >= 3 of the first segment.
 * sieve_hi : upper bound on numbers we will ever sieve (for sqrt check).
 */
static seg_ctx_t *seg_create(int64_t start_odd, int64_t block_size,
                              int64_t sieve_hi)
{
    seg_ctx_t *ctx = (seg_ctx_t *)calloc(1, sizeof(seg_ctx_t));
    ctx->block_size = block_size;
    ctx->buf = (uint8_t *)malloc((size_t)(block_size/8 + 2));
    if (!ctx->buf) { fprintf(stderr,"OOM: seg buf\n"); exit(1); }

    int64_t sq = (int64_t)sqrt((double)sieve_hi) + 2;
    int64_t n_bp = 0;
    for (int64_t i = 1; i < g_small_count && g_small_prime[i] <= sq; i++)
        n_bp++;
    ctx->n_bp   = n_bp;
    ctx->bp_val = (int64_t *)malloc(sizeof(int64_t)*(size_t)(n_bp+1));
    ctx->bp_off = (int64_t *)malloc(sizeof(int64_t)*(size_t)(n_bp+1));
    if (!ctx->bp_val || !ctx->bp_off) {
        fprintf(stderr,"OOM: bp arrays\n"); exit(1);
    }

    for (int64_t i = 0; i < n_bp; i++) {
        int64_t p = g_small_prime[i+1];
        ctx->bp_val[i] = p;
        /* First odd multiple of p that is >= start_odd and >= p^2 */
        int64_t first = p * p;
        if (first < start_odd) {
            int64_t r = start_odd % p;
            first = (r == 0) ? start_odd : start_odd + (p - r);
            if (first % 2 == 0) first += p;
        }
        ctx->bp_off[i] = (first - start_odd) / 2;
    }
    ctx->seg_low = start_odd;

    /* Fill first segment */
    int64_t slots = block_size;
    memset(ctx->buf, 0, (size_t)(slots/8+2));
    for (int64_t i = 0; i < n_bp; i++) {
        int64_t p   = ctx->bp_val[i];
        int64_t off = ctx->bp_off[i];
        for (; off < slots; off += p)
            ctx->buf[off>>3] |= (uint8_t)(1<<(off&7));
        ctx->bp_off[i] = off - slots;
    }
    return ctx;
}

static void seg_next(seg_ctx_t *ctx)
{
    ctx->seg_low += 2 * ctx->block_size;
    int64_t slots = ctx->block_size;
    memset(ctx->buf, 0, (size_t)(slots/8+2));
    for (int64_t i = 0; i < ctx->n_bp; i++) {
        int64_t p   = ctx->bp_val[i];
        int64_t off = ctx->bp_off[i];
        for (; off < slots; off += p)
            ctx->buf[off>>3] |= (uint8_t)(1<<(off&7));
        ctx->bp_off[i] = off - slots;
    }
}

static inline int seg_prime(const seg_ctx_t *ctx, int64_t j)
{
    return !((ctx->buf[j>>3] >> (j&7)) & 1);
}

static void seg_free(seg_ctx_t *ctx)
{
    if (!ctx) return;
    free(ctx->buf);
    free(ctx->bp_val);
    free(ctx->bp_off);
    free(ctx);
}

/* ================================================================== */
/*  Count primes in [lo_odd, hi_odd] (both odd, lo >= 3)             */
/*  Uses a caller-provided buffer of block_size/8+2 bytes.            */
/* ================================================================== */
static int64_t count_primes_odd(int64_t lo, int64_t hi, int64_t block_size,
                                 uint8_t *buf, int64_t sieve_hi)
{
    if (lo > hi) return 0;
    int64_t sq = (int64_t)sqrt((double)sieve_hi) + 2;
    int64_t n_bp = 0;
    for (int64_t i = 1; i < g_small_count && g_small_prime[i] <= sq; i++)
        n_bp++;

    int64_t count = 0;
    for (int64_t seg = lo; seg <= hi; seg += 2*block_size) {
        int64_t seg_end = seg + 2*block_size - 2;
        if (seg_end > hi) seg_end = hi;
        int64_t slots = (seg_end - seg)/2 + 1;

        memset(buf, 0, (size_t)(slots/8+2));

        for (int64_t i = 0; i < n_bp; i++) {
            int64_t p = g_small_prime[i+1];
            int64_t r = seg % p;
            int64_t start = (r == 0) ? seg : seg + (p - r);
            if (start % 2 == 0) start += p;
            if (start < p*p) {
                start = p*p;
                if (start % 2 == 0) start += p;
            }
            for (int64_t k = start; k <= seg_end; k += 2*p) {
                int64_t j = (k - seg)/2;
                if (j < slots) buf[j>>3] |= (uint8_t)(1<<(j&7));
            }
        }

        for (int64_t j = 0; j < slots; j++)
            if (!((buf[j>>3] >> (j&7)) & 1)) count++;
    }
    return count;
}

/* ================================================================== */
/*  Modular obstruction filter (unchanged from v4/v5)                 */
/* ================================================================== */
#define N_FILTER_PRIMES 4
static const int64_t FILTER_Q[N_FILTER_PRIMES]  = {3, 7, 11, 19};
static const int64_t FILTER_Q2[N_FILTER_PRIMES] = {9, 49, 121, 361};

static inline int vq_is_odd(int64_t s, int qi)
{
    int64_t q  = FILTER_Q[qi];
    int64_t q2 = FILTER_Q2[qi];
    if (s % q  != 0) return 0;
    if (s % q2 != 0) return 1;
    int e = 0;
    while (s % q == 0) { s /= q; e++; }
    return e & 1;
}

static inline uint8_t compute_pm_qmask(int64_t pm)
{
    uint8_t mask = 0;
    for (int i = 0; i < N_FILTER_PRIMES; i++)
        if (pm == FILTER_Q[i]) mask |= (uint8_t)(1<<i);
    return mask;
}

static inline int modular_obstruction_filter(int64_t s, uint8_t pm_qmask)
{
    for (int i = 0; i < N_FILTER_PRIMES; i++) {
        if (pm_qmask & (1<<i)) continue;
        if (vq_is_odd(s, i)) return 1;
    }
    return 0;
}

/* ================================================================== */
/*  Progress reporting interval (segments between status lines)       */
/* ================================================================== */
#define REPORT_INTERVAL 500

/* ================================================================== */
/*  Global progress state (updated atomically by worker threads)      */
/* ================================================================== */
static volatile int64_t g_segs_done   = 0;
static          int64_t g_segs_total  = 0;
static          time_t  g_t_parallel  = 0;

/* ================================================================== */
/*  find_nth_prime: compute p_n exactly (single-threaded, called once)*/
/*                                                                     */
/*  Uses the same PNT-skip strategy as the old per-thread seek, but   */
/*  executed once in main() before the parallel region.  Result is    */
/*  passed to all search_case() invocations, eliminating 49 redundant */
/*  full-range sieve scans (~58 s each on a 16-thread machine).       */
/*                                                                     */
/*  target_n : 1-based prime index to locate.                         */
/*  Returns  : p_{target_n}, or -1 on failure.                        */
/* ================================================================== */
static int64_t find_nth_prime(int64_t target_n, int64_t block_size)
{
    if (target_n <= g_small_count) return g_small_prime[target_n - 1];

    int64_t sieve_hi = (int64_t)((double)pnt_pn(target_n) * 1.05) + 1000000LL;
    int64_t p_est    = pnt_pn(target_n);
    int64_t win_lo   = p_est - 8 * block_size * 2;
    if (win_lo < 3) win_lo = 3;
    if (win_lo % 2 == 0) win_lo++;

    /* Count primes < win_lo */
    int64_t pi_before = 1;  /* prime 2 */
    int64_t small_top = (win_lo - 2 < SMALL_SIEVE_LIMIT)
                        ? (win_lo - 2) : SMALL_SIEVE_LIMIT;
    for (int64_t i = 1; i < g_small_count; i++) {
        if (g_small_prime[i] > small_top) break;
        pi_before++;
    }
    if (win_lo - 2 > SMALL_SIEVE_LIMIT + 2) {
        int64_t lo2 = SMALL_SIEVE_LIMIT + 2;
        if (lo2 % 2 == 0) lo2++;
        int64_t hi2 = win_lo - 2;
        if (hi2 % 2 == 0) hi2--;
        if (lo2 <= hi2) {
            uint8_t *cbuf = (uint8_t *)malloc((size_t)(block_size/8+2));
            if (!cbuf) { fprintf(stderr,"OOM: cbuf in find_nth_prime\n"); exit(1); }
            pi_before += count_primes_odd(lo2, hi2, block_size, cbuf, sieve_hi);
            free(cbuf);
        }
    }

    /* Walk forward from win_lo, counting until we hit target_n */
    seg_ctx_t *ctx   = seg_create(win_lo, block_size, sieve_hi);
    int64_t n_cur    = pi_before;
    int64_t target_p = -1;

    while (ctx->seg_low <= sieve_hi) {
        int64_t slots     = block_size;
        int64_t last_slot = slots;
        if (ctx->seg_low + 2*(slots-1) > sieve_hi)
            last_slot = (sieve_hi - ctx->seg_low)/2 + 1;
        for (int64_t j = 0; j < last_slot; j++) {
            if (!seg_prime(ctx, j)) continue;
            n_cur++;
            if (n_cur == target_n) {
                target_p = ctx->seg_low + 2*j;
                seg_free(ctx);
                return target_p;
            }
        }
        seg_next(ctx);
    }
    seg_free(ctx);
    return -1;  /* should not happen if sieve_hi is large enough */
}

/* ================================================================== */
/*  Per-case search function                                          */
/*                                                                     */
/*  start_prime: exact value of p_{start_n}, pre-computed in main().  */
/*  The segment is created directly at start_prime; no seek overhead. */
/* ================================================================== */
static int64_t search_case(int m, int64_t pm,
                            int64_t start_n, int64_t start_prime,
                            int64_t bound_n, int64_t block_size)
{
    int64_t pm2      = pm * pm;
    uint8_t pm_qmask = compute_pm_qmask(pm);
    int64_t sieve_hi = (int64_t)((double)pnt_pn(bound_n) * 1.05) + 1000000LL;

    /* Start segment exactly at p_{start_n}; n_cur begins at start_n - 1
     * so the first prime encountered increments it to start_n.           */
    seg_ctx_t *ctx = seg_create(start_prime, block_size, sieve_hi);
    int64_t n_cur  = start_n - 1;
    int64_t result = -1;

    while (ctx->seg_low <= sieve_hi) {
        int64_t slots     = block_size;
        int64_t last_slot = slots;
        if (ctx->seg_low + 2*(slots-1) > sieve_hi)
            last_slot = (sieve_hi - ctx->seg_low)/2 + 1;

        for (int64_t j = 0; j < last_slot; j++) {
            if (!seg_prime(ctx, j)) continue;

            n_cur++;
            if (n_cur > bound_n) goto done;

            int64_t pn = ctx->seg_low + 2*j;
            int64_t s  = (int64_t)m + n_cur;

            if (modular_obstruction_filter(s, pm_qmask)) continue;

            int64_t r   = pm2 % s;
            int64_t pnr = pn  % s;
            /* pnr < s can reach ~5e9 for BOUND=5e9; pnr^2 would overflow
             * int64_t (max ~9.2e18).  Use __int128 for the squaring step.
             * This is the same guard that v6b handled implicitly via its
             * tighter bound (s <= 2e9+1e5, pnr^2 < 4e18 < INT64_MAX).     */
            int64_t pnr_sq_mod = (int64_t)(((__int128)pnr * pnr) % s);
            if ((r + pnr_sq_mod) % s == 0) {
                result = n_cur;
                goto done;
            }
        }
        seg_next(ctx);

        /* Progress report: every REPORT_INTERVAL segments, one thread
         * prints a status line showing overall % done and ETA.
         * The atomic increment is cheap (~1 ns); the fprintf only fires
         * when the threshold is crossed, so there is no lock contention. */
        {
            int64_t done;
            #pragma omp atomic capture
            done = ++g_segs_done;

            if (done % REPORT_INTERVAL == 0 && g_segs_total > 0) {
                double pct = 100.0 * (double)done / (double)g_segs_total;
                if (pct > 100.0) pct = 100.0;   /* clamp: estimate may be low */
                double elapsed = difftime(time(NULL), g_t_parallel);
                double eta_s   = (pct > 0.1 && pct < 100.0)
                                 ? elapsed * (100.0 - pct) / pct : 0.0;
                #pragma omp critical
                fprintf(stderr,
                    "[progress] %5.1f%%  elapsed %4.0f min  ETA %4.0f min"
                    "  (seg %lld / %lld)\n",
                    pct, elapsed/60.0, eta_s/60.0,
                    (long long)done, (long long)g_segs_total);
            }
        }
    }

done:
    seg_free(ctx);
    return result;
}

/* ================================================================== */
/*  Main                                                               */
/* ================================================================== */
int main(int argc, char *argv[])
{
    int64_t BOUND      = DEFAULT_BOUND;
    int64_t BLOCK_SIZE = DEFAULT_BLOCK_SIZE;

    if (argc > 1) BOUND      = atoll(argv[1]);
    if (argc > 2) BLOCK_SIZE = atoll(argv[2]);

    int64_t START_N = PHASE6_BOUND + 1LL;

    time_t t0 = time(NULL);
    fprintf(stderr,
        "sun41_v12c -- targeted segmented-sieve search (A247975)\n");
    fprintf(stderr, "  Unresolved cases : %d (after fase 2, n-bound 4e9)\n",
            N_UNRESOLVED);
    fprintf(stderr, "  Search range     : n-index %lld .. %lld\n",
            (long long)START_N, (long long)BOUND);
    fprintf(stderr, "  Segment size     : %lld odd slots (%.2f MB/thread)\n",
            (long long)BLOCK_SIZE, (double)BLOCK_SIZE/8e6);
#ifdef _OPENMP
    fprintf(stderr, "  OpenMP threads   : %d\n", omp_get_max_threads());
#endif
    fprintf(stderr, "\n");

    fprintf(stderr, "Building small prime sieve up to %ld ...\n",
            SMALL_SIEVE_LIMIT);
    build_small_sieve(SMALL_SIEVE_LIMIT);

    int max_m = UNRESOLVED_M[N_UNRESOLVED - 1];
    if ((int64_t)max_m > g_small_count) {
        fprintf(stderr,"ERROR: small sieve has %ld primes, need index %d\n",
                g_small_count, max_m);
        exit(1);
    }
    fprintf(stderr, "  p_%d = %ld (largest p_m)\n\n",
            max_m, (long)g_small_prime[max_m - 1]);

    int64_t p_bound_est = pnt_pn(BOUND);
    fprintf(stderr, "  Estimated p_{BOUND} ~ %.4e\n\n", (double)p_bound_est);

    /* Pre-compute p_{START_N} exactly, once, before the parallel region.
     * This avoids 49 redundant full-range sieve scans (one per thread).
     * Cost: ~60 s serial, but paid only once.                            */
    fprintf(stderr, "Pre-computing p_{START_N} (n = %lld) ...\n",
            (long long)START_N);
    int64_t p_start_exact = find_nth_prime(START_N, BLOCK_SIZE);
    if (p_start_exact < 0) {
        fprintf(stderr, "FATAL: could not locate p_{START_N}.\n");
        exit(1);
    }
    fprintf(stderr, "  p_{%lld} = %lld\n\n",
            (long long)START_N, (long long)p_start_exact);

    int64_t result_n[N_UNRESOLVED];
    for (int i = 0; i < N_UNRESOLVED; i++) result_n[i] = -1;

    int completed = 0;

    /* Estimate total segments: each segment covers BLOCK_SIZE odd slots,
     * containing approximately BLOCK_SIZE / ln(p_{START_N}) primes.
     * We need n_range primes per case, so segs_per_case ~ n_range * ln(p) / BLOCK_SIZE.
     * p_{START_N} ~ 4.7e10 for START_N=2e9 -> ln ~ 24.6; use p_start_exact directly. */
    int64_t n_range      = BOUND - START_N + 1;
    double ln_p_start    = log((double)p_start_exact);
    int64_t segs_per_case = (int64_t)((double)n_range * ln_p_start / (double)BLOCK_SIZE) + 10;
    g_segs_total = segs_per_case * N_UNRESOLVED;
    g_t_parallel = time(NULL);

    fprintf(stderr,
        "Starting parallel search (%d cases, ~%lld segments total) ...\n\n",
        N_UNRESOLVED, (long long)g_segs_total);

    /* ---- Parallel search ---- */
    /* Note: no reduction clause; completed is updated atomically so that
     * the progress log inside the critical section always shows the true
     * global count, not a per-thread partial value.                       */
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < N_UNRESOLVED; i++) {
        int     m  = UNRESOLVED_M[i];
        int64_t pm = g_small_prime[m - 1];

        int64_t n = search_case(m, pm, START_N, p_start_exact, BOUND, BLOCK_SIZE);
        result_n[i] = n;

        int local_completed;
        #pragma omp atomic capture
        local_completed = ++completed;

        #pragma omp critical
        {
            if (n >= 0)
                fprintf(stderr,
                    "  FOUND  m=%-7d  p_m=%-10ld  n=%lld  [%d/%d]\n",
                    m, (long)pm, (long long)n, local_completed, N_UNRESOLVED);
            else
                fprintf(stderr,
                    "  open   m=%-7d  p_m=%-10ld  n > %lld  [%d/%d]\n",
                    m, (long)pm, (long long)BOUND, local_completed, N_UNRESOLVED);
            fflush(stderr);
        }
    }

    int n_res = 0, n_open = 0;
    for (int i = 0; i < N_UNRESOLVED; i++) {
        if (result_n[i] >= 0) n_res++; else n_open++;
    }

    fprintf(stderr, "\n=== SUMMARY ===\n");
    fprintf(stderr, "Search range   : n = %lld .. %lld\n",
            (long long)START_N, (long long)BOUND);
    fprintf(stderr, "Newly resolved : %d / %d\n", n_res, N_UNRESOLVED);
    fprintf(stderr, "Still open     : %d\n", n_open);
    fprintf(stderr, "Elapsed        : %.0f s\n", difftime(time(NULL), t0));

    /* ---- CSV output ---- */
    printf("# sun41_v12c: targeted search for A247975 unresolved cases\n");
    printf("# Phase-6 bound : n <= 2000000000\n");
    printf("# This run      : n <= %lld\n", (long long)BOUND);
    printf("# Generated     : %s", ctime(&t0));
    printf("m,p_m,min_n,status\n");

    for (int i = 0; i < N_UNRESOLVED; i++) {
        int     m  = UNRESOLVED_M[i];
        int64_t pm = g_small_prime[m - 1];
        int64_t n  = result_n[i];
        if (n >= 0)
            printf("%d,%ld,%lld,RESOLVED\n",   m, (long)pm, (long long)n);
        else
            printf("%d,%ld,NOT_FOUND,n_leq_%lld\n",
                   m, (long)pm, (long long)BOUND);
    }

    free(g_small_prime);
    free(g_base_bsieve);
    return 0;
}
