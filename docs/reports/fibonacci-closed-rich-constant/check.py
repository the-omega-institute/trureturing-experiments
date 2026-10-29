#!/usr/bin/env python3
"""Exact finite checks for Maity--Puzynina arXiv:2605.19535v1, C5.3.

Only integers and Fraction arithmetic are used. No finite test proves C5.3.
"""
from fractions import Fraction
from functools import cache


# Source indexing: f[-1]=b, f[0]=a, f[n]=f[n-1]f[n-2].
f = {-1: "b", 0: "a"}
for k in range(1, 17):
    f[k] = f[k - 1] + f[k - 2]
F = {k: len(w) for k, w in f.items()}


@cache
def closed(w):
    """Literal prefix/suffix and no internal occurrence definition."""
    if len(w) <= 1:
        return True
    for size in range(1, len(w)):
        v = w[:size]
        if w.endswith(v) and w.find(v, 1, len(w) - 1) == -1:
            return True
    return False


def closed_by_occurrences(w):
    """Independent count of overlapping occurrences of each proper border."""
    if len(w) <= 1:
        return True
    for size in range(1, len(w)):
        v = w[:size]
        positions = [i for i in range(len(w) - size + 1)
                     if w[i:i + size] == v]
        if positions == [0, len(w) - size]:
            return True
    return False


def subfactors(w):
    return {w[i:j] for i in range(len(w) + 1)
            for j in range(i, len(w) + 1)}


@cache
def cl(w):
    return sum(closed(v) for v in subfactors(w))


def factors(n):
    """Complete factor set from the three legal pairs of substitution tiles.

    Under a->ab,b->a, level-k tiles are f[k], f[k-1]. Choose both at
    least n long. An n-window meets at most two consecutive tiles. The
    legal pairs aa,ab,ba all occur in f[3]=abaab; bb never occurs. Thus
    their level-k images cover every n-factor and all occur in f[k+3].
    """
    if n == 0:
        return {""}
    k = 1
    while F[k - 1] < n:
        k += 1
    a, b = f[k], f[k - 1]
    result = {w[i:i + n] for w in (a + a, a + b, b + a)
              for i in range(len(w) - n + 1)}
    prefix = f[k + 3]
    direct = {prefix[i:i + n] for i in range(len(prefix) - n + 1)}
    assert result == direct, (n, "tile/prefix disagreement")
    assert len(result) == n + 1, (n, "Sturmian complexity mismatch")
    return result


def qadd(x, y):
    return x[0] + y[0], x[1] + y[1]


def qmul(x, y):
    return x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def qinv(x):
    den = x[0] * x[0] - 5 * x[1] * x[1]
    assert den != 0
    return x[0] / den, -x[1] / den


def qpow(x, n):
    result = (Fraction(1), Fraction(0))
    for _ in range(n):
        result = qmul(result, x)
    return result


def above_conjectured_constant(r):
    # r >= (2+5*sqrt(5))/121 iff t=121*r-2 >= 0 and t^2 >= 125.
    t = 121 * r - 2
    return t >= 0 and t * t >= 125


def main():
    phi = (Fraction(1, 2), Fraction(1, 2))
    five = (Fraction(5), Fraction(0))
    three = (Fraction(3), Fraction(0))
    two = (Fraction(2), Fraction(0))
    numerator = qadd(qmul(five, phi), three)
    denominator = qmul(qpow(phi, 2), qpow(qadd(qpow(phi, 3), two), 2))
    constant = qmul(numerator, qinv(denominator))
    assert constant == (Fraction(2, 121), Fraction(5, 121))
    assert qmul(qpow(phi, 3), qinv(qpow(qadd(qpow(phi, 3), two), 2))) == constant
    print("EXACT_CONSTANT (2 + 5*sqrt(5))/121")
    assert closed("") and closed("a") and closed("aaa")
    assert not closed("ab")
    table_p = [2,3,4,5,7,9,11,13,16,19,22,25,28,31,36,41,46,
               51,56,61,66,71,76,81,89,97,105,113,121,129,137,145,153]
    table_m = [2,3,4,5,6,9,10,13,15,17,22,24,26,31,34,37,40,
               48,56,59,62,65,73,81,86,91,96,101,106,119,132,145,150]
    all_words = set()
    best = None
    for n in range(1, 65):
        words = factors(n)
        all_words.update(words)
        counts = {w: cl(w) for w in words}
        m = min(counts.values())
        ratio = Fraction(m, n * n)
        assert above_conjectured_constant(ratio), (n, m, "counterexample")
        if best is None or ratio < best[0]:
            best = ratio, n, m
        if n <= 33:
            assert cl(f[16][:n]) == table_p[n - 1], (n, "Table 1 prefix")
            assert m == table_m[n - 1], (n, "Table 1 minimum")
        if n in (1, 5, 6, 13, 21, 33, 48, 64):
            print(f"SMALL n={n} factors={len(words)} M_n={m} ratio={ratio}")
        if n == 48:
            minimizers = [w for w in words if counts[w] == m]
            special = f[4] + f[7] + f[4][:-2]
            assert minimizers == [special]
            assert special == special[::-1]
    for w in all_words:
        assert closed(w) == closed_by_occurrences(w), (w, "predicate disagreement")
    print(f"FIDELITY Table_1_rows=33 exact_prefix_counts_and_minima=OK")
    print(f"PREDICATES distinct_words={len(all_words)} overlapping_occurrence_crosscheck=OK")
    print(f"SMALL_BOUND all_factors_lengths=1..64 checked_factors={len(all_words)} "
          f"least_ratio={best[0]} at_length={best[1]} Cl={best[2]}")
    for k in range(7, 13):
        w = f[k - 3] + f[k] + f[k - 3][:-2]
        count = cl(w)
        predicted = 1 + F[k - 3] * F[k]
        assert count == predicted, (k, count, predicted)
        assert len(w) == F[k] + 2 * F[k - 3] - 2
        assert w == w[::-1]
        assert w in f[k + 3], (k, "source factor membership")
        assert above_conjectured_constant(Fraction(count, len(w) ** 2))
        print(f"SPECIAL k={k} length={len(w)} Cl={count} ratio={Fraction(count,len(w)**2)}")
    print("SCOPE finite checks support C5.3; the universal lower bound and exact infimum remain unproved")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
