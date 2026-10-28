#!/usr/bin/env python3
"""Exact finite checks for Peter Bala's A000073 prime-power conjecture.

Source: https://oeis.org/A000073, %N and %F (Bala, Dec 28 2022),
revision #742; eligible primes: https://oeis.org/A106282, %N, revision #18.
Only Python integer arithmetic is used. These checks are not a proof.
"""
from math import isqrt

SOURCE_TERMS = (
    0, 0, 1, 1, 2, 4, 7, 13, 24, 44, 81, 149, 274, 504, 927, 1705,
    3136, 5768, 10609, 19513, 35890, 66012, 121415, 223317, 410744,
    755476, 1389537, 2555757, 4700770, 8646064, 15902591, 29249425,
    53798080, 98950096, 181997601, 334745777, 615693474, 1132436852,
)
SOURCE_PRIMES = (
    3, 5, 23, 31, 37, 59, 67, 71, 89, 97, 113, 137, 157, 179, 181,
    191, 223, 229, 251, 313, 317, 331, 353, 367, 379, 383, 389, 433,
    443, 449, 463, 467, 487, 509, 521, 577, 619,
)


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def eligible(p):
    return prime(p) and all((x**3 - x**2 - x - 1) % p for x in range(p))


def recurrence_terms(length, modulus=None):
    a = [0, 0, 1][:length]
    while len(a) < length:
        value = sum(a[-3:])
        a.append(value if modulus is None else value % modulus)
    return a


def qmul(a, b, q):
    """Multiplication in (Z/qZ)[t]/(t^3-t^2-t-1), low coefficients first."""
    c = [0] * 5
    for i in range(3):
        for j in range(3):
            c[i + j] += a[i] * b[j]
    for d in (4, 3):
        for shift in (1, 2, 3):
            c[d - shift] += c[d]
    return tuple(x % q for x in c[:3])


def qpow(a, n, q):
    result = (1, 0, 0)
    while n:
        if n & 1:
            result = qmul(result, a, q)
        a = qmul(a, a, q)
        n //= 2
    return result


def jump(n, q):
    # The t^2 coefficient has initial values 0,0,1 and the source recurrence.
    return qpow((0, 1, 0), n, q)[2]


def mmul(a, b, q):
    return tuple(tuple(sum(a[i][h] * b[h][j] for h in range(3)) % q
                       for j in range(3)) for i in range(3))


def matrix_jump(n, q):
    # State (a(n),a(n+1),a(n+2)); independent companion-matrix computation.
    result = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    base = ((0, 1, 0), (0, 0, 1), (1, 1, 1))
    while n:
        if n & 1:
            result = mmul(result, base, q)
        base = mmul(base, base, q)
        n //= 2
    return result[0][2]


def indices(p, n, k):
    return (n * p**(k - 1), n * p**k, n * p**(k + 1))


def main():
    assert tuple(recurrence_terms(len(SOURCE_TERMS))) == SOURCE_TERMS
    primes = [p for p in range(2, 620) if eligible(p)]
    assert tuple(primes) == SOURCE_PRIMES
    print(f'SOURCE_FIDELITY terms={len(SOURCE_TERMS)} eligible_primes_below_620={len(primes)} OK')

    algorithm_cases = 0
    for q in (2, 3, 5, 9, 25, 49, 125):
        terms = recurrence_terms(301, q)
        for n, value in enumerate(terms):
            assert jump(n, q) == matrix_jump(n, q) == value, (q, n)
            algorithm_cases += 1
    print(f'ALGORITHM_CROSSCHECK direct_recurrence=polynomial=matrix cases={algorithm_cases} OK')

    # Full rectangular grid; includes n divisible by p, k=1 and p=3.
    rectangular = 0
    for p in primes:
        for k in range(1, 7):
            q = p**k
            for n in range(1, 21):
                ix = indices(p, n, k)
                assert sum(jump(m, q) for m in ix) % q == 0, (p, n, k)
                rectangular += 1
    print(f'CONJECTURE_GRID primes=all_eligible_below_620 n=1..20 k=1..6 cases={rectangular} OK')

    # Additional large n, including p^2, p^3; two independent jump methods.
    extra = 0
    for p in primes:
        for n in (37, 128, p, p**2, p**3, p + 1):
            for k in (1, 2, 4, 6):
                q = p**k
                vals = [jump(m, q) for m in indices(p, n, k)]
                assert vals == [matrix_jump(m, q) for m in indices(p, n, k)]
                assert sum(vals) % q == 0, (p, n, k)
                extra += 1
    print(f'EXTRA_GRID two_algorithms n=37,128,p,p^2,p^3,p+1 k=1,2,4,6 cases={extra} OK')

    direct = 0
    for p in (3, 5):
        for k in range(1, 4):
            q = p**k
            terms = recurrence_terms(5 * p**(k + 1) + 1, q)
            for n in range(1, 6):
                assert sum(terms[m] for m in indices(p, n, k)) % q == 0
                direct += 1
    print(f'INDEPENDENT_DIRECT_RECURRENCE p=3,5 n=1..5 k=1..3 cases={direct} OK')

    a = recurrence_terms(28)
    assert (a[3], a[9], a[27]) == (1, 44, 2555757)
    assert a[3] + a[9] + a[27] == 9 * 283978
    print('EXAMPLE p=3 n=1 k=2: 1+44+2555757=2555802=9*283978 OK')

    # Controls catch the index convention and the essential prime hypothesis.
    assert a[1:4] == [0, 1, 1]
    bad_p2 = sum(jump(m, 2) for m in indices(2, 1, 1)) % 2
    assert not eligible(2) and bad_p2 == 1
    bad_p7 = sum(jump(m, 7) for m in indices(7, 1, 1)) % 7
    assert not eligible(7) and bad_p7 != 0
    print(f'NEGATIVE_CONTROLS ineligible_p=2 residual={bad_p2}; ineligible_p=7 residual={bad_p7} OK')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
