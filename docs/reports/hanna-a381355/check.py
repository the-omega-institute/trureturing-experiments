#!/usr/bin/env python3
"""Exact finite checks for Hanna's OEIS A381353/A381355 conjecture.

Uses only Python's standard library. No floating-point calculations, network
access, or external packages are needed. Finite checks are not a universal proof.
"""
from fractions import Fraction

# OEIS DATA, A381353 (offset 0) and A381355 (offset 1).
# Paul D. Hanna; https://oeis.org/A381353 and https://oeis.org/A381355.
# These short source excerpts retain OEIS attribution and CC BY-SA 4.0 terms.
SOURCE_F = [1, 1, -1, 2, -5, 13, -31, 48, 129, -2035, 12963, -20703,
            -782282, 14675113, -177056253, 1716591959]
SOURCE_A = [1, -3, 10, -35, 121, -390, 1037, -1083, -14030, 137837,
            -382106, -8791718, 199408912, -2701500413, 28888970650]


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def next_prime(n):
    n += 1
    while not is_prime(n):
        n += 1
    return n


def prime_list(count):
    ps, n = [], 1
    while len(ps) < count:
        n = next_prime(n)
        ps.append(n)
    return ps


def multiply(a, b, degree):
    out = [0] * (degree + 1)
    for i, ai in enumerate(a[:degree + 1]):
        for j, bj in enumerate(b[:degree + 1 - i]):
            out[i + j] += ai * bj
    return out


def power(a, exponent, degree):
    out = [1] + [0] * degree
    base = list(a[:degree + 1])
    while exponent:
        if exponent % 2:
            out = multiply(out, base, degree)
        exponent //= 2
        if exponent:
            base = multiply(base, base, degree)
    return out


def construct(schedule, slope, degree):
    f = [1, slope]
    for n in range(2, degree + 1):
        # Compute the rational correction before testing any integrality claim.
        q = Fraction(-power(f, schedule[n], n)[n], schedule[n])
        f.append(q.numerator if q.denominator == 1 else q)
    return f


def inverse(f, degree):
    assert f[0] == 1
    g = [1]
    for n in range(1, degree + 1):
        g.append(-sum(f[j] * g[n - j] for j in range(1, n + 1)))
    return g


def log_coefficients(f, degree):
    # Direct formal multiplication F' * inverse(F), shifted by X.
    g = inverse(f, degree)
    return [0] + [sum(j * f[j] * g[n - j] for j in range(1, n + 1))
                  for n in range(1, degree + 1)]


def verify(schedule, slope, degree):
    assert all(is_prime(schedule[n]) and schedule[n] > n
               for n in range(2, degree + 1))
    f = construct(schedule, slope, degree)
    assert all(Fraction(c).denominator == 1 for c in f)
    a = log_coefficients(f, degree)
    assert multiply(f, inverse(f, degree), degree) == [1] + [0] * degree
    # Independent logarithmic-derivative recurrence from F * (XF'/F) = XF'.
    recursive_a = [0]
    for n in range(1, degree + 1):
        recursive_a.append(n * f[n] - sum(recursive_a[j] * f[n - j]
                                          for j in range(1, n)))
    assert a == recursive_a
    for n in range(2, degree + 1):
        p = schedule[n]
        h = power(f, p, n)
        assert h[n] == 0
        assert all(h[j] % p == 0 for j in range(1, n + 1))
        assert all(p * sum(a[i] * h[k - i] for i in range(1, k + 1))
                   == k * h[k] for k in range(1, n + 1))
        assert a[n] % p == 0
        # Nonzero rational perturbations check the coefficient-correction law.
        prefix = f[:n]
        for delta in [Fraction(1), Fraction(-2), Fraction(1, 2)]:
            changed = prefix + [f[n] + delta]
            assert power(changed, p, n)[n] == p * delta
    return f, a


def main():
    degree = 64
    ps = [0] + prime_list(degree)
    f, a = verify(ps, 1, degree)
    assert f[:len(SOURCE_F)] == SOURCE_F
    assert a[1:len(SOURCE_A) + 1] == SOURCE_A
    print('FIDELITY A381353: 16 DATA coefficients match, indices 0..15')
    print('FIDELITY A381355: 15 DATA coefficients match, indices 1..15')
    # Also reproduce the first displayed coefficient table in A381353.
    expected_rows = {
        3: [1, 3, 0, 1, -3, 6, -1, -93, 834],
        5: [1, 5, 5, 0, 0, -4, 30, -190, 1145],
        7: [1, 7, 14, 7, 0, -7, 35, -188, 1050],
    }
    assert all(power(f, p, 8) == row for p, row in expected_rows.items())
    print('FIDELITY A381353: displayed power rows p=3,5,7 match through degree 8')
    print('SOURCE: integral F, diagonal zero equations, low-degree Frobenius residues,')
    print('        logarithmic-derivative identity and prime divisibility pass for n=2..64')
    print('SOURCE initial quotients:', [a[n] // ps[n] for n in range(2, 13)])
    family_degree = 24
    schedules = {
        'nth-prime': [0] + prime_list(family_degree),
        'least-prime-above-n': [0] + [next_prime(n) for n in range(1, family_degree + 1)],
        'least-prime-above-2n': [0] + [next_prime(2 * n) for n in range(1, family_degree + 1)],
    }
    slopes = [-3, -1, 0, 1, 2, 3]
    for name, schedule in schedules.items():
        for slope in slopes:
            verify(schedule, slope, family_degree)
        print('FAMILY:', name, 'slopes', slopes, 'pass for n=2..24')
    print('BOUNDARY: a(1)=1 is not divisible by p_1=2; the source excludes n=1')
    assert a[1] == 1 and a[1] % ps[1] != 0
    print('No universal proof or Lean proof is asserted by these finite checks.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
