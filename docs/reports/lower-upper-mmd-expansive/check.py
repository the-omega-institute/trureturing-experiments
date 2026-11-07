#!/usr/bin/env python3
"""Exact finite checks for the lower/upper mmd-expansive separation.

Only integers and fractions are used. These checks validate finite instances
of the estimates in the proposed general argument; they do not establish
any infinite-dimensional limit or constitute a Lean proof.
"""
from fractions import Fraction as F
from itertools import combinations, product


def N(j):
    return 2 ** (2 ** j)


def active(j):
    r = 0
    while N(2 * r) < j:
        if j <= N(2 * r + 1):
            return True
        r += 1
    return False


def A(k):
    total, r = 0, 0
    while N(2 * r) < k:
        total += min(k, N(2 * r + 1)) - N(2 * r)
        r += 1
    return total


def prefixes(q, m, k):
    assert 0 <= m <= k
    pts = [F(0)]
    for j in range(m + 1, k + 1):
        if active(j):
            digit = F(q - 1, q ** j)
            pts += [x + digit for x in pts]
    return sorted(pts)


def residue_weights(p):
    """Exact sum of 2^(-|i|) in each residue class modulo p."""
    den = 2 ** p - 1
    return [F(1) + F(2, den)] + [
        F(2 ** (p - r) + 2 ** r, den) for r in range(1, p)
    ]


def rho_periodic(x, y, t):
    p = len(x)
    assert len(y) == p
    return sum((w * abs(x[(r + t) % p] - y[(r + t) % p])
                for r, w in enumerate(residue_weights(p))), F(0))


def rho_word(x, y, t):
    """Sequences equal to these finite words on [0,n), and zero elsewhere."""
    return sum((F(1, 2 ** abs(j - t)) * abs(a - b)
                for j, (a, b) in enumerate(zip(x, y))), F(0))


def fidelity():
    """Source Lemma 2.1: both cylinder inclusions, on exact periodic data."""
    cases = 0
    for p in range(1, 5):
        assert sum(residue_weights(p), F(0)) == 3
        for eps in [F(1, 3), F(1), F(3, 2)]:
            for x in product([F(0), F(1, 2), F(1)], repeat=p):
                choices = [[max(F(0), a - eps / 3), a,
                            min(F(1), a + eps / 3)] for a in x]
                for y in product(*choices):
                    assert all(abs(a - b) <= eps / 3 for a, b in zip(x, y))
                    assert all(rho_periodic(x, y, t) <= eps for t in range(p))
                    cases += 1
            grid = list(product([F(0), F(1, 2), F(1)], repeat=p))
            for x in grid:
                for y in grid:
                    if all(rho_periodic(x, y, t) <= eps for t in range(p)):
                        assert all(abs(a - b) <= eps for a, b in zip(x, y))
                        cases += 1
    print(f'FIDELITY: source Lemma 2.1, periods 1..4, {cases} exact cases; weight sum = 3')


def endpoint_densities():
    for k in range(301):
        assert A(k) == sum(active(j) for j in range(1, k + 1))
    cases = 0
    for r in range(4):
        lo, hi = N(2 * r), N(2 * r + 1)
        previous = 0 if r == 0 else N(2 * r - 1)
        assert F(A(lo), lo) <= F(previous, lo)
        assert F(A(hi), hi) >= 1 - F(lo, hi)
        for m in range(min(lo, 100) + 1):
            assert 0 <= F(A(lo) - A(m), lo) <= F(previous, lo)
            assert F(A(hi) - A(m), hi) >= 1 - F(lo + A(m), hi)
            cases += 1
        print(f'ENDPOINT: r={r}, A(N_even)/N_even={F(A(lo),lo)}, '
              f'A(N_odd)/N_odd={F(A(hi),hi)}')
    print(f'DENSITY: A(k) checked for k=0..300; {cases} tail endpoint cases')


def digit_geometry():
    cases, points = 0, 0
    for q in range(3, 10):
        for k in range(1, 25):
            for m in [v for v in [0, 3, 8, 18] if v <= k]:
                pts = prefixes(q, m, k)
                assert len(pts) == 2 ** (A(k) - A(m))
                assert len(set(pts)) == len(pts)
                assert all(b - a >= F(q - 2, q ** k)
                           for a, b in zip(pts, pts[1:]))
                assert all(0 <= x <= F(1, q ** m) for x in pts)
                # Every prefix cylinder is contained in [prefix,prefix+q^(-k)].
                for suffix in prefixes(q, k, k + 3):
                    assert 0 <= suffix <= F(1, q ** k)
                # The entire infinite geometric tail, not a floating approximation.
                first = F(q - 1, q ** (k + 1))
                assert first / (1 - F(1, q)) == F(1, q ** k)
                cases += 1
                points += len(pts)
    print(f'DIGITS: q=3..9, k=1..24, m in {{0,3,8,18}} with m<=k; '
          f'{cases} prefix families, {points} points')


def shift_geometry():
    pairs, cover_cases = 0, 0
    for q in range(3, 6):
        for k in range(3, 7):
            for m in range(4):
                alphabet = prefixes(q, m, k)
                eps = F(q - 2, 2 * q ** k)
                for n in range(1, 4):
                    words = list(product(alphabet, repeat=n))
                    assert len(words) == 2 ** ((A(k) - A(m)) * n)
                    for x, y in combinations(words, 2):
                        assert max(rho_word(x, y, t) for t in range(n)) > eps
                        pairs += 1
    # Scale band q^(-(k+1)) < epsilon <= q^(-k): prefix depth k+4.
    for q in range(3, 10):
        for k in range(1, 25):
            for eps in [F(1, q ** k),
                        (F(1, q ** k) + F(1, q ** (k + 1))) / 2]:
                L = 1
                while F(2, 2 ** L) >= eps / 2:
                    L += 1
                cylinder_diam = F(1, q ** (k + 4))
                assert cylinder_diam < eps / 12
                assert 3 * cylinder_diam + F(2, 2 ** L) < eps
                assert F(q - 2, q ** (k - 1)) > eps
                for n in range(1, 9):
                    assert len(range(-L, n + L)) == n + 2 * L
                cover_cases += 1
    print(f'SHIFT: {pairs} independent-word separated pairs, lengths 1..3; '
          f'{cover_cases} exact finite-window cover budgets, lengths 1..8')


def small_bowen_balls():
    cases = 0
    for q in range(3, 10):
        for delta in [F(1), F(1, 2), F(1, 7), F(1, 100), F(1, 1000)]:
            m = 0
            while F(3, q ** m) > delta:
                m += 1
            # This symbolic tail alphabet is nontrivial: the first future active digit exists.
            j = m + 1
            while not active(j):
                j += 1
            assert j > m and F(q - 1, q ** j) <= F(1, q ** m)
            alphabet = prefixes(q, m, max(j, 18))[:4]
            for p in range(1, 4):
                zero = (F(0),) * p
                for word in product(alphabet, repeat=p):
                    for t in range(p):
                        assert rho_periodic(word, zero, t) <= F(3, q ** m) <= delta
                        cases += 1
    print(f'LOCAL: {cases} exact periodic tail-orbit checks, q=3..9, '
          'delta in {1,1/2,1/7,1/100,1/1000}; no loss of future active digits')


def main():
    fidelity()
    endpoint_densities()
    digit_geometry()
    shift_geometry()
    small_bowen_balls()
    print('LIMITS: finite checks only; unbounded limits and compactness require the general argument')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
