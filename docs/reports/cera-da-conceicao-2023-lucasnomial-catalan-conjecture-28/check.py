"""Checks for Conjecture 28 of J. Cera Da Conceicao, "Multinomial Catalan Numbers and Lucas Analogues",
arXiv:2410.04857v1 (J. Integer Sequences 26 (2023), Article 23.9.7).

U = U(P, Q): U_0 = 0, U_1 = 1, U_{j+2} = P U_{j+1} - Q U_j.  t!_U = U_1 ... U_t.
C_{U,l}(n) = ((l+1)n)!_U / ((n!_U)^{l+1} U_{n+1}^l).  Conjecture 28: integral for a set of n of density 1.

(a) Exact rational arithmetic for n <= NEXACT agrees with the valuation criterion used below.
(b) Valuation criterion: by the source's Theorem 14 the multinomial is integral, and by Theorem 19 only primes
    p not dividing Q with rank rho_U(p) <= l can make C non-integral; these primes divide U_2 ... U_l.
(c) Lemma 25 / Table 3: for the listed pairs (U, l) the exceptional n are the ones the lemma states.
(d) Density readings: the proportion of 1 <= n <= N with C integral, for pairs outside Table 3.
"""
import sys
from fractions import Fraction
from math import gcd

NEXACT = 40


def lucas(P, Q, T):
    u = [0, 1]
    while len(u) <= T:
        u.append(P * u[-1] - Q * u[-2])
    return u


def exact_C(P, Q, l, n):
    u = lucas(P, Q, (l + 1) * n + 2)
    num = 1
    for j in range(1, (l + 1) * n + 1):
        num *= u[j]
    den = 1
    for j in range(1, n + 1):
        den *= u[j]
    den = den ** (l + 1) * u[n + 1] ** l
    return Fraction(num, den)


def prime_factors(x):
    x = abs(x)
    out = set()
    d = 2
    while d * d <= x:
        while x % d == 0:
            out.add(d)
            x //= d
        d += 1
    if x > 1:
        out.add(x)
    return out


def obstructing_primes(P, Q, l):
    u = lucas(P, Q, l + 1)
    ps = set()
    for j in range(2, l + 1):
        ps |= prime_factors(u[j])
    return sorted(p for p in ps if Q % p != 0)


def valuations(P, Q, p, T, K=60):
    mod = p ** K
    a, b = 0, 1
    v = [None, 0]
    for _ in range(2, T + 1):
        a, b = b, (P * b - Q * a) % mod
        if b == 0:
            v.append(K)
        else:
            e = 0
            x = b
            while x % p == 0:
                x //= p
                e += 1
            v.append(e)
    return v


def integral_flags(P, Q, l, N):
    T = (l + 1) * N + 2
    flags = [True] * (N + 1)
    for p in obstructing_primes(P, Q, l):
        v = valuations(P, Q, p, T)
        pre = [0] * (T + 1)
        for j in range(1, T + 1):
            pre[j] = pre[j - 1] + v[j]
        for n in range(1, N + 1):
            num = pre[(l + 1) * n] - (l + 1) * pre[n]
            if num < l * v[n + 1]:
                flags[n] = False
    return flags


def regular(P, Q):
    if P == 0 or Q == 0 or gcd(P, Q) != 1 or P * P - 4 * Q == 0:
        return False
    return all(x != 0 for x in lucas(P, Q, 12)[1:])


def main():
    bad = 0
    pairs = [(1, -1), (1, 2), (2, -1), (3, 2), (1, -3), (5, 3), (2, 3), (4, 3), (3, -1)]
    print("(a) exact versus valuation criterion, n <= %d" % NEXACT)
    for (P, Q) in pairs:
        assert regular(P, Q), (P, Q)
        for l in range(1, 7):
            flags = integral_flags(P, Q, l, NEXACT)
            mism = [n for n in range(1, NEXACT + 1) if (exact_C(P, Q, l, n).denominator == 1) != flags[n]]
            bad += len(mism) > 0
            if mism:
                print(f"  MISMATCH U({P},{Q}) l={l} n={mism}")
    print("  done")
    print("(c) Lemma 25 exceptions for Table 3 pairs, n <= 3000")
    # Lemma 25: n = p - 1 if p | Delta and l = p; n = p - 2 if rho_U(p) = p - 1 and l = p; n = 5 for U(1,-1), l = 6.
    table = [((1, -1), 4, []), ((1, -1), 5, [4]), ((1, -1), 6, [5]), ((1, 2), 4, []), ((1, 2), 6, []),
             ((1, 2), 7, [6]), ((1, 2), 8, []), ((1, 3), 5, [3]), ((2, 3), 3, []), ((1, -5), 3, [2])]
    for (P, Q), l, expected in table:
        flags = integral_flags(P, Q, l, 3000)
        fails = [n for n in range(1, 3001) if not flags[n]]
        ok = fails == expected
        bad += not ok
        print(f"  U({P},{Q}) l={l}: non-integral n = {fails} {'OK' if ok else 'MISMATCH'}")
    print("(d) proportion of integral C_{U,l}(n), 1 <= n <= N, pairs outside Table 3")
    for (P, Q), l in [((3, 2), 2), ((3, 2), 3), ((1, -1), 7), ((2, -1), 4), ((5, 3), 3), ((4, 3), 6), ((3, -1), 5)]:
        row = []
        flags = integral_flags(P, Q, l, 100000)
        for N in (1000, 10000, 100000):
            k = sum(flags[1:N + 1])
            row.append(f"N={N}: {k / N:.5f}")
        nonint = [n for n in range(1, 100001) if not flags[n]]
        print(f"  U({P},{Q}) l={l} primes={obstructing_primes(P, Q, l)}: " + ", ".join(row) +
              f"; non-integral count {len(nonint)}, first {nonint[:6]}")
    print(f"EXIT={1 if bad else 0}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
