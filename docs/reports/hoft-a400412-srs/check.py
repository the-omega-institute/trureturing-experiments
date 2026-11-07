#!/usr/bin/env python3
"""Exact finite checks for OEIS A400412, not a universal proof or Lean certificate.

Run with Python 3. The only arithmetic used as evidence is integer arithmetic and
fractions.Fraction. A missing bounded-search witness is printed as unknown.
"""
from array import array
from collections import Counter
from fractions import Fraction
from math import isqrt

BOUND = 2_000_000
ROWS = 60


def primes_upto(n):
    return [p for p in range(2, n + 1)
            if all(p % d for d in range(2, isqrt(p) + 1))]


def trial_divisors(n):
    ds = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            ds.append(d)
            if d * d != n:
                ds.append(n // d)
    return sorted(ds)


def block_count(ds, threshold=2):
    assert ds and ds == sorted(set(ds))
    return 1 + sum(b > threshold * a for a, b in zip(ds, ds[1:]))


def full_binary_divisors(odd_ds, e):
    return sorted(d * 2**j for d in odd_ds for j in range(e + 1))


def check_source_fidelity():
    # A400412 DATA, first nine rows (search to 3025, as in its Mathematica code).
    data = [
        [1], [6, 3], [18, 0, 9], [30, 78, 15, 21],
        [162, 0, 0, 0, 81], [90, 666, 45, 75, 63, 147],
        [1458, 0, 0, 0, 0, 0, 729],
        [210, 1830, 135, 105, 165, 189, 357, 903],
        [450, 0, 225, 0, 1225, 0, 441, 0, 3025],
    ]
    minima = {}
    for v in range(1, 3026):
        ds = trial_divisors(v)
        n = sum(d % 2 for d in ds)
        k = block_count(ds)
        if n <= 9:
            minima.setdefault((n, k), v)
    for n, row in enumerate(data, 1):
        for k, expected in enumerate(row, 1):
            if expected:
                assert minima[(n, k)] == expected, (n, k, expected)
            else:
                assert (n, k) not in minima
    assert block_count(trial_divisors(10)) == 2  # A237271 example
    assert trial_divisors(10) == [1, 2, 5, 10]
    print('FIDELITY: all 45 first-nine-row DATA entries matched source search to 3025')
    print('FIDELITY: zero entries mean no witness to 3025, not an experimental emptiness proof')


def check_interlacing():
    qs = [p for p in primes_upto(43) if p != 2]
    rs = [p for p in primes_upto(199) if p != 2]
    order_checks = full_checks = 0
    for h in range(2, 17):
        for q in qs:
            for r in rs:
                if q == r:
                    continue
                A = [q**j for j in range(h)]
                R = [r*q**j for j in range(h)]
                odd_ds = sorted(A + R)
                assert len(set(odd_ds)) == 2*h
                ratios = Counter(Fraction(b, a) for a, b in zip(odd_ds, odd_ds[1:]))
                if r > A[-1]:
                    expected = Counter({Fraction(q): 2*h-2})
                    expected[Fraction(r, A[-1])] += 1
                    t = None
                else:
                    t = next(t for t in range(h-1) if q**t < r < q**(t+1))
                    alpha, beta = Fraction(r, q**t), Fraction(q**(t+1), r)
                    expected = Counter()
                    expected[Fraction(q)] += 2*t
                    expected[alpha] += h-t
                    expected[beta] += h-t-1
                    expected += Counter()  # Drop zero multiplicities.
                assert ratios == expected, (h, q, r, ratios, expected)
                order_checks += 1
                for e in range(7):
                    C = 2**(e+1)
                    ds = full_binary_divisors(odd_ds, e)
                    B = block_count(ds)
                    gap = block_count(odd_ds, C)
                    formula = 1 + sum(mult for ratio, mult in expected.items() if ratio > C)
                    assert B == gap == formula, (h, q, r, e, B, gap, formula)
                    assert not (B % 2 == 0 and 2 < B < h), (h, q, r, e, B)
                    full_checks += 1
    power_checks = 0
    for p in [2, 3, 5, 7, 11, 13, 17]:
        for q in rs:
            ds = [q**j for j in range(2*p)]
            for e in range(7):
                B = block_count(full_binary_divisors(ds, e))
                assert B == (1 if q <= 2**(e+1) else 2*p)
                power_checks += 1
    print(f'INTERLACING: {order_checks} exact adjacent-ratio multisets, {full_checks} full divisor lists passed')
    print(f'PRIME_POWER: {power_checks} full divisor lists passed')


def bounded_triangle():
    # Smallest-prime-factor sieve, then each odd part m and every 2^e*m <= BOUND.
    # This covers every positive integer <= BOUND exactly once.
    spf = array('I', range(BOUND + 1))
    for p in range(2, isqrt(BOUND) + 1):
        if spf[p] == p:
            for j in range(p*p, BOUND+1, p):
                if spf[j] == j:
                    spf[j] = p
    minima = {}
    covered = formula_checks = 0
    for m in range(1, BOUND + 1, 2):
        z, ds = m, [1]
        while z > 1:
            p = spf[z]
            old, power = ds[:], 1
            while z % p == 0:
                z //= p
                power *= p
                ds.extend(d*power for d in old)
        ds.sort()
        n = len(ds)
        v, C, e = m, 2, 0
        while v <= BOUND:
            actual_ds = full_binary_divisors(ds, e)
            k = block_count(actual_ds)
            assert k == block_count(ds, C), (v, n, k)
            if n <= ROWS:
                key = (n, k)
                if key not in minima or v < minima[key]:
                    minima[key] = v
            # Independently compare actual divisor lists on a second finite window.
            if v <= 20000:
                assert block_count(trial_divisors(v)) == k, (v, n, k)
                formula_checks += 1
            covered += 1
            v *= 2
            C *= 2
            e += 1
    assert covered == BOUND
    assert formula_checks == 20000
    prime_rows = primes_upto(ROWS // 2)
    for p in prime_rows:
        for k in range(4, p, 2):
            assert (2*p, k) not in minima, (p, k, minima[(2*p, k)])
    counterexamples = [(n, k, v) for (n, k), v in sorted(minima.items())
                       if k >= 3 and v % 2 == 0]
    print(f'TRIANGLE: exhaustive 1 <= v <= {BOUND}, 1 <= n <= {ROWS}; {len(minima)} exact positive minima')
    print(f'GAP_FORMULA: actual full divisor lists agreed for every v <= {BOUND}')
    print(f'DIVISOR_GENERATION: independent trial divisor lists agreed for every v <= {formula_checks}')
    print('CONJECTURE_2: no forbidden witness in the exhaustive integer window')
    print('CONJECTURE_1: even exact positive minima with k>=3:', counterexamples)
    print('TABLE: each row lists k=1..n; ? denotes uncomputed (>BOUND or empty), never zero')
    for n in range(1, ROWS+1):
        print(f'n={n}: ' + ','.join(str(minima[(n,k)]) if (n,k) in minima else '?' for k in range(1,n+1)))
    print(f'FAMILY_CHECK: {len(counterexamples)} bounded counterexamples; no universal minimum-parity conclusion')
    return counterexamples


def main():
    check_source_fidelity()
    check_interlacing()
    counterexamples = bounded_triangle()
    # A counterexample would be reported rather than hidden by an assertion.
    print('EXPECTED_MAIN_SETTLEMENT: Proved (informal universal argument, finite checks only here)')
    print('FAMILY_STATUS:', 'refutation candidate' if counterexamples else 'no counterexample found; proof route unresolved')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
