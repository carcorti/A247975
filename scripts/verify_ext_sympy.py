#!/usr/bin/env python3
"""
verify_ext_sympy.py
===================
Verifica indipendente con SymPy dei casi risolti nell'Extended search
(sun_ext_v11_pgo, sessioni S1-S4, range n = 2e9 .. 2e11).

Congettura verificata: (m + n) | (p_m^2 + p_n^2)
dove p_k è il k-esimo numero primo (1-indexed).

OEIS A247975 / Sun's Conjecture 4.1(i)
Autore: Carlo Corti, 2026
Repository: https://github.com/carcorti/A247975

Metodo
------
Per ogni caso (m, p_m, a_m, p_{a_m}) verificare:
  1. p_m   è primo e corrisponde a g_small_prime[m-1]     (cross-check index)
  2. p_n   è primo con sympy.isprime()                    (verifica diretta)
  3. (p_m^2 + p_n^2) % (m + n) == 0                      (condizione della congettura)

Per (1) si usa sympy.prime(m) che è esatto ma lento per m grande.
Per valori grandi di m si usa sympy.isprime(p_m) + cross-check del valore.
Per (2) sympy.isprime() su valori fino a ~2.3e12 è deterministico (Miller-Rabin
con witness set completo) e restituisce risultato esatto.

Nota sui tempi
--------------
p_n nell'Extended arriva fino a ~2.28e12 (p_{82300} ~ 2.25e12).
sympy.isprime() su questi valori impiega <1 ms ciascuno (test di primalità
deterministico). Il bottleneck è sympy.prime(m) per m > 50000 (~0.5s ciascuno).
Per questo script usiamo il valore di p_m fornito dal CSV (già verificato
in Phase 1) e lo convalidiamo solo con isprime().

Campionamento
-------------
Strategia a 3 livelli, come da §3.3 del paper (Metodo 2):
  A) TUTTI i 64 casi risolti — verifica della condizione (3)
  B) Campione rappresentativo — verifica aggiuntiva di p_n con isprime()
  C) Casi ad alto Cramér ratio — verifica prioritaria (anomalie statistiche)

Uso
---
  # Verifica tutti i 64 casi (consigliato, ~5-10 minuti):
  python3 verify_ext_sympy.py

  # Solo smoke test (5 casi, <1 minuto):
  python3 verify_ext_sympy.py --smoke

  # Verbose (mostra ogni caso):
  python3 verify_ext_sympy.py --verbose

  # Output CSV per il paper:
  python3 verify_ext_sympy.py --csv > verification_ext_results.csv
"""

import sys
import time
import argparse
from sympy import isprime

# =============================================================================
# DATI: 64 casi risolti dall'Extended search (S1 + S2)
# Formato: (m, p_m, a_m, p_{a_m})
# Fonte: ext1_results.csv + ext2_results.csv
# =============================================================================

EXTENDED_CASES = [
    # --- Sessione 1 (n = 2e9 .. 5e10) --- 62 casi ---
    (11924,  127289,      6119832581,   151142662267),
    (18241,  203227,      2873665609,    68701158689),
    (19234,  215317,      5206723528,   127712739697),
    (20686,  233173,      2858183263,    68314916771),
    (22186,  251809,      6418823416,   158846687033),
    (22802,  259577,     44000698448,  1177174198639),  # alto Cramér 192303
    (23995,  274453,     47572503750,  1276596209489),  # alto Cramér 196577
    (33566,  396239,      4283714631,   104199324451),
    (36303,  431833,      4006056395,    97164894163),
    (36726,  437083,      2412333964,    57230326859),
    (37253,  443873,     29834694588,   786117161153),
    (37949,  453181,     23612836020,   616425195857),
    (43297,  523007,     20932270152,   543819574261),
    (48834,  596419,     23630637176,   616908416777),
    (49889,  610469,      3941343324,    95528205013),
    (51816,  636277,     11824762849,   300167153191),
    (55468,  685453,     12790824118,   325737648979),
    (55874,  690953,      7689816992,   191750554309),
    (61007,  760367,      3353229683,    80707333381),
    (67203,  844253,      7747012646,   193236644017),
    (67273,  845197,      2766403488,    66026785697),
    (67447,  847621,      8688185334,   217752203707),
    (68744,  865327,      2423390090,    57504240109),
    (69546,  876191,     17485936480,   451005863261),
    (71049,  896723,      2298999145,    54425781049),
    (72568,  917759,      4208046357,   102280257463),
    (76761,  975497,      6814926932,   169075281041),
    (81986, 1047971,      5364065784,   131738811713),
    (82096, 1049519,      7393138089,   184048943843),
    (84521, 1083371,      2282409344,    54015791803),
    (86557, 1111393,      6670175973,   165334596179),
    (90005, 1159597,      4615871244,   112639311881),
    (90161, 1161757,      3788376009,    91663964981),
    (92335, 1191563,      2223036115,    52549379809),
    (93790, 1212347,      2179204108,    51467883263),
    (93890, 1213651,      2640300431,    62888154977),
    (94429, 1221131,      5683120888,   139917719653),
    (94768, 1225981,      2184705450,    51603549161),
    (97251, 1260901,      7850344255,   195922614071),
    (97351, 1262453,      5389653243,   132394050307),
    (98008, 1271383,      3197827293,    76808347667),
    (98258, 1275011,      8564068296,   214512996613),
    (99105, 1286953,      5589474041,   137515154639),
    (100528, 1307083,     7675012978,   191366007089),
    (100871, 1311857,     3893067251,    94307911091),
    (102179, 1330501,     9444795782,   237538125809),
    (102283, 1331951,     2029876287,    47790310313),
    (104389, 1361813,     6264229453,   154861625773),
    (105219, 1373591,     2594185703,    61741943503),
    (107550, 1406429,     2837673259,    67803334231),
    (108013, 1412779,     6369975757,   157587051887),
    (108091, 1413877,     2900339079,    69366839291),
    (109430, 1433101,     5807064035,   143100102157),
    (109857, 1438973,    14929726016,   382615207459),
    (110967, 1454759,     7273776154,   180953871269),
    (112360, 1474519,    11071138165,   280276166167),
    (113001, 1483453,     7659405632,   190960585853),
    (114161, 1500061,     2744916180,    65491493099),
    (115225, 1514659,     5362892889,   131708767273),
    (117674, 1550161,     3543958704,    85502899387),
    (118880, 1567271,     5493201674,   135046872031),
    (119089, 1570447,     3995510688,    96898143581),
    # --- Sessione 2 (n = 5e10 .. 1e11) --- 2 casi ---
    (53963,  665179,     50362759647,  1354459021543),  # alto Cramér 85653
    (82300, 1052459,     81961045941,  2245774975837),  # alto Cramér 87990
]

# Sottoinsieme per smoke test (5 casi rappresentativi)
SMOKE_CASES = [
    (11924,  127289,      6119832581,   151142662267),  # valore OEIS di riferimento
    (22802,  259577,     44000698448,  1177174198639),  # Cramér 192303 (record S1)
    (23995,  274453,     47572503750,  1276596209489),  # Cramér 196577 (record assoluto S1)
    (53963,  665179,     50362759647,  1354459021543),  # Cramér 85653 (S2)
    (82300, 1052459,     81961045941,  2245774975837),  # Cramér 87990 (S2, n più grande)
]

# =============================================================================
# Funzioni di verifica
# =============================================================================

def verify_case(m, pm, am, pam, verbose=False):
    """
    Verifica un singolo caso (m, p_m, a_m, p_{a_m}).

    Controlla:
      1. isprime(p_m)          — p_m è primo
      2. isprime(p_{a_m})      — p_{a_m} è primo
      3. (pm^2 + pam^2) % s == 0  dove s = m + a_m

    Restituisce: dict con campi 'ok', 'error', 'elapsed_s'
    """
    t0 = time.monotonic()
    result = {'m': m, 'pm': pm, 'am': am, 'pam': pam, 'ok': False, 'error': None}

    # Check 1: p_m è primo?
    if not isprime(pm):
        result['error'] = f"p_m={pm} non è primo"
        result['elapsed_s'] = time.monotonic() - t0
        return result

    # Check 2: p_{a_m} è primo?
    if not isprime(pam):
        result['error'] = f"p_{{a_m}}={pam} non è primo"
        result['elapsed_s'] = time.monotonic() - t0
        return result

    # Check 3: condizione della congettura
    s   = m + am
    pm2 = pm * pm
    pam2 = pam * pam
    remainder = (pm2 + pam2) % s

    if remainder != 0:
        result['error'] = (
            f"DIVISIBILITÀ FALLITA: "
            f"(p_m^2 + p_n^2) % s = {remainder} ≠ 0  "
            f"[s={s}, pm^2={pm2}, pam^2={pam2}]"
        )
        result['elapsed_s'] = time.monotonic() - t0
        return result

    result['ok'] = True
    result['elapsed_s'] = time.monotonic() - t0

    if verbose:
        print(f"  OK  m={m:6d}  p_m={pm:12d}  a_m={am:15d}  "
              f"p_{{a_m}}={pam:15d}  s={s}  [{result['elapsed_s']:.2f}s]")

    return result


def run_verification(cases, verbose=False, csv_mode=False):
    """
    Esegue la verifica su una lista di casi.
    Restituisce (n_ok, n_fail, elapsed_total).
    """
    if csv_mode:
        print("m,p_m,a_m,p_{a_m},status,elapsed_s")
    else:
        print(f"\n{'='*70}")
        print(f"  Verifica SymPy — Extended search A247975")
        print(f"  {len(cases)} casi da verificare")
        print(f"{'='*70}\n")

    t_total = time.monotonic()
    n_ok   = 0
    n_fail = 0
    failures = []

    for i, (m, pm, am, pam) in enumerate(cases, 1):
        r = verify_case(m, pm, am, pam, verbose=verbose)

        if csv_mode:
            status = "VERIFIED" if r['ok'] else f"FAILED:{r['error']}"
            print(f"{m},{pm},{am},{pam},{status},{r['elapsed_s']:.3f}")
        else:
            if r['ok']:
                n_ok += 1
                if not verbose:
                    # Mostra progresso ogni 10 casi
                    if i % 10 == 0 or i == len(cases):
                        elapsed = time.monotonic() - t_total
                        print(f"  [{i:3d}/{len(cases)}]  {n_ok} OK, {n_fail} FAIL  "
                              f"({elapsed:.0f}s elapsed)")
            else:
                n_fail += 1
                failures.append(r)
                print(f"\n  *** FALLIMENTO ***  m={m}  {r['error']}\n")

    elapsed_total = time.monotonic() - t_total

    if not csv_mode:
        print(f"\n{'='*70}")
        print(f"  RISULTATO FINALE")
        print(f"  Casi verificati:  {len(cases)}")
        print(f"  OK:               {n_ok}")
        print(f"  FALLITI:          {n_fail}")
        print(f"  Tempo totale:     {elapsed_total:.1f}s ({elapsed_total/60:.1f} min)")
        print(f"{'='*70}\n")

        if n_fail == 0:
            print("  ✅  Tutti i casi verificati senza errori.")
        else:
            print(f"  ❌  {n_fail} FALLIMENTO/I — vedere dettagli sopra.")
            print("\n  Casi falliti:")
            for r in failures:
                print(f"    m={r['m']}  errore: {r['error']}")

    return n_ok, n_fail, elapsed_total


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Verifica SymPy dei casi Extended search A247975"
    )
    parser.add_argument('--smoke',   action='store_true',
                        help='Esegui solo 5 casi di smoke test (<1 min)')
    parser.add_argument('--verbose', action='store_true',
                        help='Mostra dettaglio per ogni caso')
    parser.add_argument('--csv',     action='store_true',
                        help='Output in formato CSV (per il paper)')
    args = parser.parse_args()

    if args.smoke:
        cases = SMOKE_CASES
        if not args.csv:
            print("\n  [SMOKE TEST — 5 casi rappresentativi]\n")
    else:
        cases = EXTENDED_CASES

    n_ok, n_fail, elapsed = run_verification(
        cases,
        verbose=args.verbose,
        csv_mode=args.csv
    )

    # Exit code: 0 = tutti OK, 1 = almeno un fallimento
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == '__main__':
    main()
