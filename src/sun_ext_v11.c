/*
 * sun_ext_v11.c  --  Extended search for unresolved cases of Sun's Conjecture 4.1(i)
 *
 * Conjecture (Z.-W. Sun, 2014-09-29, arXiv:1309.1679):
 *   For every m in Z+, there exists n in Z+ such that (m+n) | (p_m^2 + p_n^2),
 *   where p_k denotes the k-th prime.
 *
 * OEIS reference: A247975.
 *
 * PURPOSE:
 *   Searches for the cases left unresolved by sun_v6 (Phase 1, n <= 2e9).
 *   Reads the list of unresolved m-values from a file (one per line).
 *   Uses a segmented sieve of Eratosthenes with per-thread buffers, keeping
 *   peak RAM to ~100 MB regardless of the search bound.
 *
 * ARITHMETIC SAFETY:
 *   pm <= p_{120000} = 1567271 < 2^21  ->  pm^2 < 2^42  (int64_t OK)
 *   pn can exceed 47e9 for n > 2e9     ->  pn itself fits int64_t (< 2^63)
 *   s = m + n < 120000 + 200e9 < 2^38  ->  pnr = pn % s < 2^38
 *   pnr^2 < 2^76: __int128 used for the squaring step (same as sun41_v12g)
 *
 * ARCHITECTURE (mirrors sun41_v12g.c):
 *   - Small prime table up to SMALL_SIEVE_LIMIT covering all p_m (m<=120000)
 *     and all sieve base primes up to sqrt(p_{BOUND}).
 *   - p_{START_N} computed once, serially, before the parallel region.
 *   - Parallel loop over unresolved cases (OpenMP, dynamic scheduling).
 *   - Each thread owns its own seg_ctx_t with private sieve buffer.
 *   - Modular obstruction filter {3,7,11,19} (identical to sun_v6).
 *
 * CHANGES v11 (from v10):
 *   Base: sun_ext_v10.c — word-scanning verified by multi-LLM review.
 *   None of these changes touch the [W1] word-scanning loop or any
 *   mathematical logic. All are code-quality / robustness fixes from
 *   the adjudication of ChatGPT, DeepSeek, and Gemini reviews.
 *
 *   [R1]  Fix strict aliasing: replace pointer cast
 *           const uint64_t *buf64 = (const uint64_t *)ctx->buf
 *         with per-word memcpy(&word, ctx->buf + w*8, 8).
 *         Standard-compliant (C11 §6.5 ¶7); compiles to identical
 *         machine code (single MOV) on GCC/Clang with -O2+.
 *   [R2]  read_m_list(): validate trailing characters after strtol.
 *         Rejects malformed lines like "123abc" (previously accepted
 *         silently as m=123).
 *   [R3]  Remove volatile from g_segs_done. Redundant with
 *         #pragma omp atomic capture (documented in [E5]).
 *   [R4]  Replace time(NULL)/difftime() with clock_gettime(CLOCK_MONOTONIC)
 *         for monotonic, sub-second timing in progress reports and summary.
 *   [R5]  sieve_hi computed once in main() and passed to search_case(),
 *         eliminating redundant pnt_pn(bound_n) calls per case.
 *   [R6]  Header COMPILE section updated with PGO instructions and
 *         recommended OMP_PROC_BIND/OMP_PLACES environment variables.
 *   [R7]  Version strings updated to v11.
 *
 * CHANGES v10 (from v9):
 *   Base: sun_ext_v9.c — verified correct against OEIS b-file.
 *
 *   [W1]  Word-scanning in search_case() inner loop.
 *         Replaces slot-by-slot seg_prime(ctx, j) with 64-bit word reads:
 *         - Read uint64_t from sieve buffer, invert (~word: 1=prime).
 *         - If word==0, skip 64 composite slots in one branch.
 *         - Enumerate set bits via __builtin_ctzll / word &= word-1.
 *         CRITICAL INVARIANT PRESERVED: n_cur++ and s++ occur exactly
 *         once per set bit (= per prime found), identical to v9 slot-by-slot.
 *         No reordering, no batching, no approximation of the count.
 *         The word-scanning skips ONLY all-composite words where no
 *         increment would have occurred anyway.
 *   [W2]  seg_create: memset padding bytes [buf_bytes, bufsz) to zero once.
 *         Ensures the last uint64_t read by word-scanning contains zero
 *         bits beyond the valid slot range. The sieving loop never writes
 *         to padding bytes (off < slots guarantees byte index < buf_bytes),
 *         so the zeros are permanent after the initial clear.
 *   [W3]  Compile line updated with -fprefetch-loop-arrays -fno-plt
 *         and -DUSE_HUGEPAGE for optimal Zen 4 performance.
 *   [W4]  Version strings updated to v10.
 *
 *   ⚠ REGRESSION TEST REQUIRED: OBS3-P1 (March 2026 lesson) showed that
 *   any restructuring of the inner loop MUST be verified against known
 *   reference values before production use.
 *
 * CHANGES v9 (from v3):
 *   Base: sun_ext_v3.c — verified correct against OEIS b-file.
 *   Loop structure UNCHANGED: slot-by-slot (seg_prime), proven safe.
 *   OBS3-P1 (word scanning) and OBS3-M1 (global g_n_bp) NOT included:
 *   both produced wrong results confirmed by SymPy verification.
 *
 *   Additions from v4-v6 reviews (all safe, do not touch inner loop):
 *   [OBS3-C1] needed bytes formula: ((block_size+7)/8) — correct ceiling.
 *   [OBS3-P2] Remove redundant first memset in seg_create.
 *   [E1]  atoll -> strtoll + range validation for CLI args (parse_int64).
 *   [E2]  Version strings updated to v9.
 *   [E3]  Guard: reject even p_{START_N} (sieve requires odd start).
 *   [E4]  Progress variable 'done' renamed 'seg_done' (label conflict).
 *   [E5]  volatile on g_segs_done documented (redundant with omp atomic).
 *   [E6]  memset in count_primes_odd: ((slots+7)/8) formula.
 *   [E7]  Cramér ratio = n/(m*ln(m)), matching sun_v6 [V3-06].
 *         Column name: cramér_ratio (directly comparable with sun_v6).
 *
 * CHANGES v3 (from v2):
 *   [OBS2-C1] aligned_alloc size rounded up to multiple of 64 (C11 UB fix).
 *   [OBS2-C2] seg_ctx_t stores buf_bytes (exact needed bytes); seg_next
 *             memset uses buf_bytes instead of padded allocation size.
 *   [OBS2-P1] divq_mod128(): inline asm DIVQ replaces __umodti3 call for
 *             __uint128_t % uint64_t. Precondition (high < s) is guaranteed
 *             algebraically; DEBUG assertion verifies at runtime.
 *   [OBS2-P2] Optional madvise(MADV_HUGEPAGE) for sieve buffers,
 *             enabled with -DUSE_HUGEPAGE compile flag.
 *   [OBS2-P3] REPORT_INTERVAL 500 -> 1000 (halves progress log frequency).
 *   [OBS2-P4] __builtin_expect(mod == 0, 0) on divisibility branch.
 *   [OBS2-P5] Incremental s variable in inner loop (s++ per prime found).
 *
 * CHANGES v2 (from v1):
 *   [OBS-C1] Divisibility check: replace (r+pnr_sq_mod)%s==0 with
 *            conditional subtract (no division, same result, cleaner).
 *   [OBS-C2] search_case() now returns a result_t {n, pn} struct,
 *            saving p_{a_m} at the moment of discovery. Eliminates
 *            find_nth_prime() calls in the post-processing CSV loop.
 *   [OBS-C5] read_m_list(): replaced atoi() with strtol() for safer
 *            input parsing.
 *   [OBS-A1] Divisibility test simplified to single __int128 division:
 *            ((__int128)pm2 + (__int128)pn*pn) % s == 0. Reduces from
 *            ~3 divisions to 1 per candidate that passes the filter.
 *   [OBS-A5] DEFAULT_BLOCK_SIZE reduced from 20,000,000 to 4,000,000
 *            (512 KB/thread × 16 = 8 MB, fits in L3 cache of 7940HS).
 *   [OBS-A6] sieve buffers allocated with aligned_alloc(64,...) for
 *            cache-line alignment.
 *
 * KEY DIFFERENCES FROM sun41_v12g.c:
 *   - Reads m-list from file (not hardcoded).
 *   - SMALL_SIEVE_LIMIT = 2500000 (covers sqrt(p_{200e9}) ~ 2.28e6).
 *   - START_N default = 2000000001 (sun_v6 Phase 1 bound + 1).
 *   - Output CSV format matches sun_v6:
 *       m, p_m, a_m, p_{a_m}, s=m+a_m, cramér_ratio
 *     Unresolved: m, p_m, UNRESOLVED, -1, -1, -1
 *   - p_{a_m} saved directly during search (no post-processing sieve).
 *
 * USAGE:
 *   ./sun_ext_v11 <input_file> [START_N [BOUND [BLOCK_SIZE]]]
 *
 *   input_file : path to file with unresolved m-values (one per line;
 *                lines starting with '#' are ignored)
 *   START_N    : first n-index to search (default: 2000000001)
 *   BOUND      : last  n-index to search (default: 200000000000)
 *   BLOCK_SIZE : odd-integer slots per sieve segment (default: 4000000)
 *
 *   Recommended environment variables for thread affinity on Zen 4:
 *     export OMP_PROC_BIND=close
 *     export OMP_PLACES=cores
 *
 * OUTPUT:
 *   stdout -> CSV (header + one row per m)
 *   stderr -> progress, timing, summary, list of still-unresolved cases
 *
 * COMPILE:
 *   gcc -O3 -march=znver4 -mtune=znver4 -funroll-loops -ffast-math \
 *       -fomit-frame-pointer -flto -fopenmp \
 *       -fprefetch-loop-arrays -fno-plt -DUSE_HUGEPAGE \
 *       sun_ext_v11.c -o sun_ext_v11 -lm
 *
 *   For additional 5-10% improvement, use Profile-Guided Optimization:
 *     Step 1: gcc [flags above] -fprofile-generate sun_ext_v11.c -o sun_ext_v11_pg -lm
 *     Step 2: ./sun_ext_v11_pg <input> START_N <short_bound>  (representative run)
 *     Step 3: gcc [flags above] -fprofile-use sun_ext_v11.c -o sun_ext_v11 -lm
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

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_HUGEPAGE
#include <sys/mman.h>
#endif

/* ================================================================== */
/*  Configuration defaults                                             */
/* ================================================================== */

#define PHASE1_BOUND       2000000000LL       /* sun_v6 Phase 1 n-bound */

#define DEFAULT_START_N    (PHASE1_BOUND + 1LL)
#define DEFAULT_BOUND      200000000000LL     /* 200e9 */
#define DEFAULT_BLOCK_SIZE   4000000LL        /* odd slots per segment; 512 KB/thread [OBS-A5] */

/* Small prime table limit.
 * Must cover:
 *   (a) p_{120000} = 1567271  (all p_m for m <= 120000)
 *   (b) sqrt(p_{200e9}): p_{200e9} ~ 5.2e12, sqrt ~ 2.28e6
 *       -> 2500000 provides a safe margin.                             */
#define SMALL_SIEVE_LIMIT  2500000L

/* Maximum number of m-values accepted from the input file */
#define MAX_CASES 500

/* ================================================================== */
/*  Small base prime table                                             */
/* ================================================================== */
static int64_t *g_small_prime = NULL;
static int64_t  g_small_count = 0;
static uint8_t *g_base_bsieve = NULL;

/* Bit-sieve macros (odd numbers only, starting from 3) */
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
    size_t bsz = (size_t)BS_BYTES(limit);
    g_base_bsieve = (uint8_t *)calloc(bsz, 1);
    if (!g_base_bsieve) { fprintf(stderr,"OOM: base bsieve\n"); exit(1); }

    for (long i = 3; i*i <= limit; i += 2)
        if (bs_is_prime(i))
            for (long j = i*i; j <= limit; j += 2*i) {
                long idx = BS_IDX(j);
                g_base_bsieve[idx>>3] |= (uint8_t)(1<<(idx&7));
            }

    int64_t cnt = 1; /* prime 2 */
    for (long k = 3; k <= limit; k += 2)
        if (bs_is_prime(k)) cnt++;

    g_small_prime = (int64_t *)malloc(sizeof(int64_t)*((size_t)cnt+64));
    if (!g_small_prime) { fprintf(stderr,"OOM: small prime table\n"); exit(1); }

    g_small_prime[g_small_count++] = 2;
    for (long k = 3; k <= limit; k += 2)
        if (bs_is_prime(k))
            g_small_prime[g_small_count++] = k;

    fprintf(stderr,"  Small prime table: %lld primes up to %ld (%.2f MB)\n",
            (long long)g_small_count, limit,
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
    size_t    buf_bytes;  /* [OBS2-C2] exact bytes needed by sieve (= block_size/8+1, rounded up) */
    int64_t   seg_low;
    int64_t   block_size;
    int64_t  *bp_val;
    int64_t  *bp_off;
    int64_t   n_bp;
} seg_ctx_t;

static seg_ctx_t *seg_create(int64_t start_odd, int64_t block_size,
                              int64_t sieve_hi)
{
    seg_ctx_t *ctx = (seg_ctx_t *)calloc(1, sizeof(seg_ctx_t));
    if (!ctx) { fprintf(stderr,"OOM: seg_ctx\n"); exit(1); }
    ctx->block_size = block_size;
    /* [OBS3-C1] Correct ceiling formula for byte count.
     * [OBS2-C1] Round up to multiple of 64 for aligned_alloc.           */
    size_t needed = ((size_t)block_size + 7) / 8;
    size_t bufsz  = (needed + 63) & ~(size_t)63;
    ctx->buf_bytes = needed;
    ctx->buf = (uint8_t *)aligned_alloc(64, bufsz);
    if (!ctx->buf) { fprintf(stderr,"OOM: seg buf\n"); exit(1); }
    /* [W2] Zero the padding bytes [buf_bytes, bufsz) once.
     * The sieving loop only writes to bytes [0, buf_bytes-1] (since
     * off < slots guarantees byte index < needed), so the padding stays
     * zero permanently.  This ensures the last uint64_t read by the
     * word-scanning loop [W1] contains zero bits beyond valid slots.    */
    if (bufsz > needed)
        memset(ctx->buf + needed, 0, bufsz - needed);
    /* [OBS3-P2] No initial full memset: the pre-fill memset below
     * covers [0, buf_bytes), and [W2] above covers the padding.        */
#ifdef USE_HUGEPAGE
    madvise(ctx->buf, bufsz, MADV_HUGEPAGE);  /* [OBS2-P2] optional THP hint */
#endif

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
    memset(ctx->buf, 0, ctx->buf_bytes);  /* [OBS2-C2] clear only used bytes */
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
    memset(ctx->buf, 0, ctx->buf_bytes);  /* [OBS2-C2] clear only used bytes */
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

        memset(buf, 0, ((size_t)slots + 7) / 8);  /* [E6] correct ceiling */

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
/*  find_nth_prime: compute p_n exactly (single-threaded)             */
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
    int64_t pi_before = 1; /* prime 2 */
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
            if (!cbuf) {
                fprintf(stderr,"OOM: cbuf in find_nth_prime\n"); exit(1);
            }
            pi_before += count_primes_odd(lo2, hi2, block_size, cbuf, sieve_hi);
            free(cbuf);
        }
    }

    /* Walk forward from win_lo counting until we hit target_n */
    seg_ctx_t *ctx = seg_create(win_lo, block_size, sieve_hi);
    int64_t n_cur  = pi_before;

    while (ctx->seg_low <= sieve_hi) {
        int64_t slots     = block_size;
        int64_t last_slot = slots;
        if (ctx->seg_low + 2*(slots-1) > sieve_hi)
            last_slot = (sieve_hi - ctx->seg_low)/2 + 1;
        for (int64_t j = 0; j < last_slot; j++) {
            if (!seg_prime(ctx, j)) continue;
            n_cur++;
            if (n_cur == target_n) {
                int64_t result = ctx->seg_low + 2*j;
                seg_free(ctx);
                return result;
            }
        }
        seg_next(ctx);
    }
    seg_free(ctx);
    return -1;
}

/* ================================================================== */
/*  Modular obstruction filter {3, 7, 11, 19} (identical to sun_v6)  */
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
    int64_t t = s;
    while (t % q == 0) { t /= q; e++; }
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
/*  Result type: n-index and prime value, returned by search_case     */
/*  [OBS-C2] Saves p_{a_m} at discovery; eliminates post-processing  */
/*           find_nth_prime() calls in the CSV output loop.           */
/* ================================================================== */
typedef struct {
    int64_t n;   /* n-index of solution, or -1 if not found */
    int64_t pn;  /* p_n at solution; valid only when n >= 0  */
} result_t;

/* ================================================================== */
/*  Progress reporting                                                 */
/* ================================================================== */
#define REPORT_INTERVAL 1000  /* [OBS2-P3] halved log frequency */

/* [R3] volatile removed: #pragma omp atomic capture already guarantees
 * visibility across threads (OpenMP 4.0+ §2.12.6). Previously kept for
 * documentation clarity [E5], but all three external reviewers flagged
 * it as potentially inhibiting optimization.                             */
static int64_t g_segs_done  = 0;
static int64_t g_segs_total = 0;
/* [R4] Monotonic clock for progress and ETA (replaces time(NULL)). */
static struct timespec g_ts_parallel;

/* ================================================================== */
/*  divq_mod128: compute (pm2 + pn*pn) % s using hardware DIVQ       */
/*                                                                     */
/*  [OBS2-P1] GCC emits __umodti3 (software, ~50 cycles) for         */
/*  __uint128_t % uint64_t. Since sum < 2^87 and s < 2^38, the high  */
/*  64-bit word h = sum >> 64 < 2^23 < s, so hardware DIVQ is valid  */
/*  (DIVQ requires h < divisor to avoid #DE fault). This replaces     */
/*  the software helper with a single ~9-19 cycle DIVQ instruction.  */
/*                                                                     */
/*  Precondition (guaranteed algebraically):                          */
/*    pm2 < 2^42, pn < 2^43 -> sum < 2^87 -> h < 2^23               */
/*    s = m + n > 120000 always -> h < 2^23 < s  QED                 */
/*                                                                     */
/*  DEBUG build verifies precondition at runtime via assertion.       */
/* ================================================================== */
static inline int divq_mod128(uint64_t pm2, uint64_t pn, uint64_t s)
{
    __uint128_t sum  = (__uint128_t)pm2 + (__uint128_t)pn * pn;
    uint64_t    high = (uint64_t)(sum >> 64);
    uint64_t    low  = (uint64_t)sum;

#ifdef DEBUG
    if (__builtin_expect(high >= s, 0)) {
        fprintf(stderr,
            "FATAL divq_mod128: precondition violated: high=%llu >= s=%llu\n",
            (unsigned long long)high, (unsigned long long)s);
        abort();
    }
#endif

    uint64_t rem;
    __asm__("divq %[d]"
            : "=d"(rem), "+a"(low)
            : "d"(high), [d]"rm"(s));
    return rem == 0;
}

/* ================================================================== */
/*  Per-case search                                                    */
/*                                                                     */
/*  Returns result_t {n, pn}: n-index and prime value of the first    */
/*  solution in [start_n, bound_n], or {-1, -1} if none found.       */
/*                                                                     */
/*  [OBS-A1]   Divisibility test: single hardware DIVQ via divq_mod128().    */
/*  [OBS2-P1]  divq_mod128() replaces __umodti3 software call.               */
/*  [OBS-C1]   No overflow risk: sum < 2*s handled inside divq_mod128.       */
/*  [OBS2-P5]  s maintained as running variable (s++ per prime found).       */
/*  [W1]       Word-scanning: process 64 sieve slots per uint64_t read.      */
/*             INVARIANT: n_cur++ and s++ occur exactly once per set bit      */
/*             (= per prime found), identical to the v9 slot-by-slot loop.    */
/* ================================================================== */
static result_t search_case(int m, int64_t pm,
                             int64_t start_n, int64_t start_prime,
                             int64_t bound_n, int64_t block_size,
                             int64_t sieve_hi)       /* [R5] computed once in main */
{
    uint64_t pm2      = (uint64_t)(pm * pm);
    uint8_t  pm_qmask = compute_pm_qmask(pm);
    result_t res      = {-1, -1};

    seg_ctx_t *ctx = seg_create(start_prime, block_size, sieve_hi);
    int64_t n_cur  = start_n - 1;
    int64_t s      = (int64_t)m + n_cur;   /* [OBS2-P5] incremental s */

    /* [W1] Number of full 64-bit words covering block_size slots.
     * Padding bytes [buf_bytes, bufsz) are guaranteed zero by [W2],
     * so reading the last word is safe even if last_slot < nwords*64. */
    const int64_t nwords_full = block_size / 64;
    const int     tail_bits   = (int)(block_size % 64);  /* 0..63 */

    while (ctx->seg_low <= sieve_hi) {
        int64_t last_slot = block_size;
        if (ctx->seg_low + 2*(block_size-1) > sieve_hi)
            last_slot = (sieve_hi - ctx->seg_low)/2 + 1;

        /* [W1] Word-scanning loop.
         * For each 64-bit word from the sieve buffer:
         *   - Invert: in the sieve, bit=0 means prime, bit=1 means composite.
         *     After ~word, bit=1 means prime.
         *   - If word==0 (all 64 slots are composite), skip entirely.
         *   - Otherwise, extract set bits one by one with ctzll + clear.
         *   - Each set bit corresponds to exactly one prime → n_cur++, s++.
         *
         * CORRECTNESS ARGUMENT:
         *   Let P = number of primes in the segment (set bits after inversion).
         *   The v9 slot-by-slot loop increments n_cur exactly P times.
         *   This word-scanning loop also increments n_cur exactly P times:
         *   each set bit is visited exactly once by the ctzll/clear cycle,
         *   and the outer word loop covers all slots [0, last_slot).
         *   The per-prime logic (filter + divq_mod128) is identical.         */
        const int64_t nw = (last_slot < block_size)
                           ? (last_slot + 63) / 64
                           : nwords_full + (tail_bits > 0 ? 1 : 0);

        for (int64_t w = 0; w < nw; w++) {
            /* [R1] memcpy instead of *(uint64_t*) pointer cast: avoids
             * strict aliasing UB (C11 §6.5 ¶7) when reading uint8_t[]
             * buffer as uint64_t. Compiles to a single MOV with -O2+. */
            uint64_t raw;
            memcpy(&raw, ctx->buf + w * 8, sizeof(raw));
            uint64_t word = ~raw;        /* invert: 1 = prime */

            /* Mask off bits beyond last_slot in the final word.
             * [W2] guarantees padding bytes are zero, so after inversion
             * they become all-ones → must be masked out explicitly.        */
            int64_t base_slot = w * 64;
            int64_t remaining = last_slot - base_slot;
            if (remaining < 64) {
                /* remaining is in [1, 63]; shift is safe (no UB) */
                word &= (1ULL << remaining) - 1;
            }

            if (word == 0) continue;   /* 64 composites → skip */

            /* Enumerate primes (set bits) in this word */
            do {
                int bit = __builtin_ctzll(word);
                word &= word - 1;      /* clear lowest set bit */

                n_cur++;
                s++;                    /* [OBS2-P5] s = m + n_cur */
                if (__builtin_expect(n_cur > bound_n, 0)) goto done;

                int64_t pn = ctx->seg_low + 2 * (base_slot + bit);

                if (modular_obstruction_filter(s, pm_qmask)) continue;

                /* [OBS2-P1] Hardware DIVQ: ~9-19 cycles vs ~50+.
                 * [OBS2-P4] __builtin_expect: solution is extremely rare. */
                if (__builtin_expect(divq_mod128(pm2, (uint64_t)pn,
                                                 (uint64_t)s), 0)) {
                    res.n  = n_cur;
                    res.pn = pn;
                    goto done;
                }
            } while (word);
        }
        seg_next(ctx);

        /* Progress report (one status line every REPORT_INTERVAL segments) */
        {
        /* [E4] Renamed 'done' -> 'seg_done' to avoid conflict with label done: */
            int64_t seg_done;
            #pragma omp atomic capture
            seg_done = ++g_segs_done;

            if (seg_done % REPORT_INTERVAL == 0 && g_segs_total > 0) {
                double pct = 100.0 * (double)seg_done / (double)g_segs_total;
                if (pct > 100.0) pct = 100.0;
                /* [R4] Monotonic clock, sub-second precision. */
                struct timespec ts_now;
                clock_gettime(CLOCK_MONOTONIC, &ts_now);
                double elapsed = (double)(ts_now.tv_sec - g_ts_parallel.tv_sec)
                               + 1e-9 * (double)(ts_now.tv_nsec - g_ts_parallel.tv_nsec);
                double eta_s   = (pct > 0.1 && pct < 100.0)
                                 ? elapsed * (100.0 - pct) / pct : 0.0;
                #pragma omp critical
                fprintf(stderr,
                    "[progress] %5.1f%%  elapsed %4.0f min  ETA %4.0f min"
                    "  (seg %lld / %lld)\n",
                    pct, elapsed/60.0, eta_s/60.0,
                    (long long)seg_done, (long long)g_segs_total);
            }
        }
    }

done:
    seg_free(ctx);
    return res;
}

/* ================================================================== */
/*  Read m-list from file                                             */
/* ================================================================== */
static int read_m_list(const char *path, int *m_list, int max_cases)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "ERROR: cannot open input file '%s'\n", path);
        exit(1);
    }
    int count = 0;
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        /* [OBS-C5] Use strtol for safe parsing with error detection */
        char *end;
        long val = strtol(line, &end, 10);
        if (end == line || val <= 0 || val > 200000) continue;
        /* [R2] Reject lines with trailing non-whitespace (e.g. "123abc"). */
        while (*end == ' ' || *end == '\t') end++;
        if (*end != '\n' && *end != '\r' && *end != '\0') continue;
        int m = (int)val;
        if (count >= max_cases) {
            fprintf(stderr,
                "ERROR: too many cases in input (max %d)\n", max_cases);
            fclose(f);
            exit(1);
        }
        m_list[count++] = m;
    }
    fclose(f);
    return count;
}

/* ================================================================== */
/*  [E1] Safe CLI argument parser (matches sun_v6 V2-04)             */
/* ================================================================== */
static int64_t parse_int64(const char *s, const char *name,
                            int64_t lo, int64_t hi)
{
    char *end;
    errno = 0;
    long long v = strtoll(s, &end, 10);
    if (errno || end == s || *end != '\0' || v < lo || v > hi) {
        fprintf(stderr, "ERROR: invalid value for %s: '%s' (must be %lld..%lld)\n",
                name, s, (long long)lo, (long long)hi);
        exit(1);
    }
    return (int64_t)v;
}

/* ================================================================== */
/*  Main                                                               */
/* ================================================================== */
int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr,
            "sun_ext_v11  --  Extended search for A247975 unresolved cases\n\n"
            "Usage: %s <input_file> [START_N [BOUND [BLOCK_SIZE]]]\n"
            "  input_file : file with unresolved m-values (one per line)\n"
            "  START_N    : first n-index to search (default: %lld)\n"
            "  BOUND      : last  n-index to search (default: %lld)\n"
            "  BLOCK_SIZE : odd slots per segment    (default: %lld)\n",
            argv[0],
            (long long)DEFAULT_START_N,
            (long long)DEFAULT_BOUND,
            (long long)DEFAULT_BLOCK_SIZE);
        return 2;
    }

    const char *input_file = argv[1];
    int64_t START_N    = (argc > 2)
        ? parse_int64(argv[2], "START_N",    2, (int64_t)4e18)
        : DEFAULT_START_N;
    int64_t BOUND      = (argc > 3)
        ? parse_int64(argv[3], "BOUND",      START_N, (int64_t)4e18)
        : DEFAULT_BOUND;
    int64_t BLOCK_SIZE = (argc > 4)
        ? parse_int64(argv[4], "BLOCK_SIZE", 10000, 100000000LL)
        : DEFAULT_BLOCK_SIZE;

    /* [R4] Monotonic wall-clock for elapsed time. */
    struct timespec ts0;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    time_t t0_wall = time(NULL);  /* wall-clock for CSV timestamp only */

    /* Read unresolved m-list */
    static int m_list[MAX_CASES];
    int n_cases = read_m_list(input_file, m_list, MAX_CASES);

    if (n_cases == 0) {
        fprintf(stderr, "ERROR: no valid m-values found in '%s'\n", input_file);
        return 1;
    }

    fprintf(stderr,
        "sun_ext_v11  --  Extended search for A247975 unresolved cases\n");
    fprintf(stderr, "  Input file       : %s\n", input_file);
    fprintf(stderr, "  Unresolved cases : %d\n", n_cases);
    fprintf(stderr, "  Search range     : n-index %lld .. %lld\n",
            (long long)START_N, (long long)BOUND);
    fprintf(stderr, "  Segment size     : %lld odd slots (%.2f MB/thread)\n",
            (long long)BLOCK_SIZE, (double)BLOCK_SIZE/8e6);
#ifdef _OPENMP
    fprintf(stderr, "  OpenMP threads   : %d\n", omp_get_max_threads());
#endif
    fprintf(stderr, "\n");

    /* Build small prime sieve */
    fprintf(stderr, "Building small prime sieve up to %ld ...\n",
            SMALL_SIEVE_LIMIT);
    build_small_sieve(SMALL_SIEVE_LIMIT);

    /* Verify small sieve covers all p_m */
    int max_m = 0;
    for (int i = 0; i < n_cases; i++)
        if (m_list[i] > max_m) max_m = m_list[i];

    if ((int64_t)max_m > g_small_count) {
        fprintf(stderr,
            "ERROR: small sieve has %lld primes, need index %d\n"
            "  Increase SMALL_SIEVE_LIMIT and recompile.\n",
            (long long)g_small_count, max_m);
        exit(1);
    }
    fprintf(stderr, "  p_%d = %lld (largest p_m needed)\n\n",
            max_m, (long long)g_small_prime[max_m - 1]);

    fprintf(stderr, "  Estimated p_{BOUND} ~ %.4e\n\n",
            (double)pnt_pn(BOUND));

    /* Pre-compute p_{START_N} exactly, once, before parallel region */
    fprintf(stderr, "Pre-computing p_{START_N} (n = %lld) ...\n",
            (long long)START_N);
    int64_t p_start_exact = find_nth_prime(START_N, BLOCK_SIZE);
    if (p_start_exact < 0) {
        fprintf(stderr, "FATAL: could not locate p_{START_N}.\n");
        exit(1);
    }
    /* [E3] Guard: segmented sieve requires an odd start prime >= 3. */
    if (p_start_exact < 3 || p_start_exact % 2 == 0) {
        fprintf(stderr,
            "FATAL: p_{START_N} = %lld is not an odd prime >= 3.\n"
            "  START_N must be >= 2.\n", (long long)p_start_exact);
        exit(1);
    }
    fprintf(stderr, "  p_{%lld} = %lld\n\n",
            (long long)START_N, (long long)p_start_exact);

    /* Allocate result arrays */
    result_t *results = (result_t *)malloc(sizeof(result_t)*(size_t)n_cases);
    if (!results) { fprintf(stderr,"OOM: results\n"); exit(1); }
    for (int i = 0; i < n_cases; i++) { results[i].n = -1; results[i].pn = -1; }

    /* [R5] Compute sieve_hi once (identical for all cases since bound_n
     * is the same). Passed to search_case() as a parameter.             */
    int64_t sieve_hi = (int64_t)((double)pnt_pn(BOUND) * 1.05) + 1000000LL;
    fprintf(stderr, "  sieve_hi         = %lld (PNT × 1.05 + 10^6)\n\n",
            (long long)sieve_hi);

    /* Estimate total segments for progress reporting */
    int64_t n_range       = BOUND - START_N + 1;
    double  ln_p_start    = log((double)p_start_exact);
    int64_t segs_per_case = (int64_t)((double)n_range * ln_p_start
                                       / (double)BLOCK_SIZE) + 10;
    g_segs_total = segs_per_case * n_cases;
    clock_gettime(CLOCK_MONOTONIC, &g_ts_parallel);  /* [R4] */

    int completed = 0;

    fprintf(stderr,
        "Starting parallel search (%d cases, ~%lld segments total) ...\n\n",
        n_cases, (long long)g_segs_total);

    /* ---- Parallel search ---- */
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < n_cases; i++) {
        int     m  = m_list[i];
        int64_t pm = g_small_prime[m - 1];

        result_t r = search_case(m, pm, START_N, p_start_exact,
                                  BOUND, BLOCK_SIZE, sieve_hi);
        results[i] = r;

        int local_completed;
        #pragma omp atomic capture
        local_completed = ++completed;

        #pragma omp critical
        {
            if (r.n >= 0)
                fprintf(stderr,
                    "  FOUND  m=%-7d  p_m=%-12lld  n=%lld  [%d/%d]\n",
                    m, (long long)pm, (long long)r.n,
                    local_completed, n_cases);
            else
                fprintf(stderr,
                    "  open   m=%-7d  p_m=%-12lld  n > %lld  [%d/%d]\n",
                    m, (long long)pm, (long long)BOUND,
                    local_completed, n_cases);
            fflush(stderr);
        }
    }

    /* ---- Summary ---- */
    int n_resolved = 0, n_open = 0;
    for (int i = 0; i < n_cases; i++) {
        if (results[i].n >= 0) n_resolved++; else n_open++;
    }

    fprintf(stderr, "\n=== SUMMARY ===\n");
    fprintf(stderr, "Input file     : %s\n", input_file);
    fprintf(stderr, "Search range   : n = %lld .. %lld\n",
            (long long)START_N, (long long)BOUND);
    fprintf(stderr, "Cases input    : %d\n", n_cases);
    fprintf(stderr, "Newly resolved : %d\n", n_resolved);
    fprintf(stderr, "Still open     : %d\n", n_open);
    /* [R4] Monotonic elapsed time. */
    {
        struct timespec ts_end;
        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        double elapsed = (double)(ts_end.tv_sec - ts0.tv_sec)
                       + 1e-9 * (double)(ts_end.tv_nsec - ts0.tv_nsec);
        fprintf(stderr, "Elapsed        : %.1f s (%.1f min)\n",
                elapsed, elapsed / 60.0);
    }

    if (n_open > 0) {
        fprintf(stderr, "\n=== UNRESOLVED CASES ===\n");
        fprintf(stderr, "%-10s  %-14s\n", "m", "p_m");
        for (int i = 0; i < n_cases; i++) {
            if (results[i].n < 0) {
                int     m  = m_list[i];
                int64_t pm = g_small_prime[m - 1];
                fprintf(stderr, "%-10d  %-14lld\n", m, (long long)pm);
            }
        }
    }

    /* ---- CSV output (format matches sun_v6) ---- */
    /* [OBS-C2] p_{a_m} taken from results[i].pn — no find_nth_prime() calls */
    /* [E7] cramér_ratio = n/(m*ln(m)), identical to sun_v6 formula.
     *      Directly comparable and concatenable with Phase 1 CSV output. */
    printf("# sun_ext_v11: extended search for A247975 unresolved cases\n");
    printf("# Phase 1 bound  : n <= %lld\n", (long long)(START_N - 1));
    printf("# This run bound : n <= %lld\n", (long long)BOUND);
    printf("# Generated      : %s", ctime(&t0_wall));
    printf("m, p_m, a_m, p_{a_m}, s=m+a_m, cramér_ratio\n");

    for (int i = 0; i < n_cases; i++) {
        int     m   = m_list[i];
        int64_t pm  = g_small_prime[m - 1];
        int64_t n   = results[i].n;
        int64_t pan = results[i].pn;

        if (n >= 0) {
            int64_t s      = (int64_t)m + n;
            /* [E7] n/(m*ln(m)), matching sun_v6 [V3-06]. Guard m>1. */
            double  cramer = (m > 1)
                             ? (double)n / ((double)m * log((double)m))
                             : 0.0;
            printf("%d, %lld, %lld, %lld, %lld, %.6f\n",
                   m, (long long)pm, (long long)n,
                   (long long)pan, (long long)s, cramer);
        } else {
            printf("%d, %lld, UNRESOLVED, -1, -1, -1\n",
                   m, (long long)pm);
        }
    }

    free(results);
    free(g_small_prime);
    free(g_base_bsieve);
    return 0;
}
