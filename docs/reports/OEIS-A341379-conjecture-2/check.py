#!/usr/bin/env python3
"""Exact, finite preregistration checks for Paul D. Hanna's OEIS A341379 (2).

No floating point and no proof assistant result are used. A direct triangular
expansion checks the source equation. An independent modulo-three recurrence
uses only the general Hurwitz Frobenius identity B^3=B(0)^3, not the proposed
six-period answer. The ordinary rational-series calculation checks the source
exponential formula without using the double-sum transformation.
"""
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
import argparse
import hashlib
import random
import urllib.request

SOURCE_INITIAL = [1, 1, 1, 4, 19, 161, 1486, 18922, 268045, 4609873, 88133356]
SOURCE_SQUARE = [1, 2, 4, 14, 76, 592, 5794, 70730, 1019456, 17176400, 329143714]
BFILE_URL = 'https://oeis.org/A341379/b341379.txt'


def pascal(nmax, modulus=None):
    rows = [[1]]
    for n in range(1, nmax + 1):
        row = [1] + [rows[-1][k - 1] + rows[-1][k] for k in range(1, n)] + [1]
        rows.append(bytearray(x % modulus for x in row) if modulus else row)
    return rows


def hmul(a, b, rows, modulus=None):
    out = []
    for n, row in enumerate(rows):
        z = sum(c * a[k] * b[n - k] for k, c in enumerate(row) if c)
        out.append(z % modulus if modulus else z)
    return out


def source_sequence(nmax, modulus=None):
    """Solve A^2=sum_(r,s) e_r e_s A^(rs), without reducing exponents.

    powers[m][n] is computed by ordinary binomial convolution with powers[m-1].
    Source degree n depends on A only below n. Division by two is exact over Z,
    and multiplication by 2 is its inverse in F_3.
    """
    rows = pascal(nmax, modulus)
    max_power = max(1, nmax * nmax // 4)
    powers = [[1] + [0] * nmax for _ in range(max_power + 1)]
    a = [1] + [0] * nmax
    square = [1] + [0] * nmax
    for n in range(1, nmax + 1):
        rhs = sum(rows[n][r] * rows[n - r][s] * powers[r * s][n - r - s]
                  for r in range(n + 1) for s in range(n - r + 1))
        interior = sum(rows[n][k] * a[k] * a[n - k] for k in range(1, n))
        if modulus:
            assert modulus == 3
            a[n] = (rhs - interior) * 2 % 3
            square[n] = rhs % 3
        else:
            assert (rhs - interior) % 2 == 0, ('integrality', n)
            a[n] = (rhs - interior) // 2
            square[n] = rhs
        powers[1][n] = a[n]
        for m in range(2, max_power + 1):
            v = sum(c * a[k] * powers[m - 1][n - k]
                    for k, c in enumerate(rows[n]) if c)
            powers[m][n] = v % modulus if modulus else v
    assert hmul(a, a, rows, modulus) == square
    return a, square


def omul(a, b, nmax):
    return [sum((a[k] * b[n - k] for k in range(n + 1)), Fraction(0))
            for n in range(nmax + 1)]


def oexp(u, nmax):
    assert u[0] == 0
    out = [Fraction(1)] + [Fraction(0)] * nmax
    for n in range(1, nmax + 1):
        out[n] = sum((k * u[k] * out[n - k] for k in range(1, n + 1)), Fraction(0)) / n
    return out


def rational_source_check(a):
    """Literally sum x^r/r! exp(x A^r), using rational ordinary coefficients."""
    nmax = len(a) - 1
    A = [Fraction(z, factorial(n)) for n, z in enumerate(a)]
    power = [Fraction(1)] + [Fraction(0)] * nmax
    rhs = [Fraction(0)] * (nmax + 1)
    for r in range(nmax + 1):
        exp_part = oexp([Fraction(0)] + power[:-1], nmax)
        for n in range(r, nmax + 1):
            rhs[n] += exp_part[n - r] / factorial(r)
        power = omul(power, A, nmax)
    assert omul(A, A, nmax) == rhs
    return rhs


def large_mod_three(nmax):
    """Evaluate all source monomials grouped by exponent modulo three.

    W[q][t] = sum_(r+s=t, rs=q mod 3) binom(t,r) mod 3.
    This grouping plus general Frobenius is valid for every normalized input.
    No formula for W, F, the answer, or the square-zero correction is used to
    generate a. In particular the residues being tested are not input data.
    """
    rows = pascal(nmax, 3)
    weights = [[0] * (nmax + 1) for _ in range(3)]
    for t, row in enumerate(rows):
        for r, c in enumerate(row):
            q = r * (t - r) % 3
            weights[q][t] = (weights[q][t] + c) % 3
    a = [1] + [0] * nmax
    square = [1] + [0] * nmax
    for n in range(1, nmax + 1):
        # W[1][0]=W[2][0]=0, so this reads only previously determined a[k].
        rhs = weights[0][n]
        for t in range(1, n + 1):
            rhs += rows[n][t] * (weights[1][t] * a[n - t] +
                                 weights[2][t] * square[n - t])
        inner = sum(rows[n][k] * a[k] * a[n - k] for k in range(1, n))
        a[n] = 2 * (rhs - inner) % 3
        square[n] = rhs % 3
    assert hmul(a, a, rows, 3) == square
    assert all(z == 0 for z in weights[2])
    E = [1] * (nmax + 1)
    E2 = [pow(2, n, 3) for n in range(nmax + 1)]
    F = weights[1]
    assert hmul(E, E, rows, 3) == E2
    assert F == [pow(2, (n - 2) // 3 + 1, 3) if n % 3 == 2 else 0
                 for n in range(nmax + 1)]
    assert hmul(F, F, rows, 3) == [0] * (nmax + 1)
    E2F = hmul(E2, F, rows, 3)
    assert E2F == [pow(2, 2 * ((n - 2) // 3) + 1, 3) if n % 3 == 2 else 0
                   for n in range(nmax + 1)]
    assert a == [(1 + 2 * (F[n] - E2F[n])) % 3 for n in range(nmax + 1)]
    assert a == [2 if n % 6 == 5 else 1 for n in range(nmax + 1)]
    assert all(z != 0 for z in a)
    assert hmul(square, a, rows, 3) == [1] + [0] * nmax
    # The nonlinear source reduction is also checked independently as a series identity.
    correction = hmul([(z - 1) % 3 if n == 0 else z for n, z in enumerate(a)], F, rows, 3)
    assert square == [(E2[n] + correction[n]) % 3 for n in range(nmax + 1)]
    return a


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bfile', type=Path, help='Optional local OEIS b-file; default downloads it')
    args = parser.parse_args()
    data = args.bfile.read_bytes() if args.bfile else urllib.request.urlopen(BFILE_URL, timeout=60).read()
    bfile = {}
    for line in data.decode().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            n, v = line.split()
            bfile[int(n)] = int(v)
    assert sorted(bfile) == list(range(126))
    exact, exact_square = source_sequence(24)
    assert exact[:11] == SOURCE_INITIAL
    assert exact_square[:11] == SOURCE_SQUARE
    assert all(exact[n] == bfile[n] for n in range(25))
    literal_rhs = rational_source_check(exact)
    assert all(literal_rhs[n] * factorial(n) == exact_square[n] for n in range(25))
    print('SOURCE_FIDELITY: all 11 displayed a(n) and all 11 displayed A(x)^2 coefficients agree')
    print('EXACT_SOURCE_EQUATION: degrees 0..24, literal rational exponential expansion agrees')
    print('EXACT_BFILE_FIDELITY: degrees 0..24 agree; b-file indices 0..125 parsed')
    direct_mod, _ = source_sequence(64, 3)
    assert [z % 3 for z in exact] == direct_mod[:25]
    rows = pascal(64, 3)
    rng = random.Random(341379)
    for trial in range(8):
        b = [rng.randrange(3) for _ in range(65)]
        assert hmul(hmul(b, b, rows, 3), b, rows, 3) == [b[0] ** 3 % 3] + [0] * 64
    print('FROBENIUS_CHECK: eight deterministic arbitrary inputs through degree 64')
    mod3 = large_mod_three(2000)
    assert mod3[:65] == direct_mod
    assert all(mod3[n] == bfile[n] % 3 for n in range(126))
    print('MOD3_SOURCE_CROSSCHECK: unreduced exponents and grouped recurrence agree at 0..64')
    print('MOD3_BFILE: all 126 residues agree')
    print('MOD3_ALL_SMALL_CASES: all 2001 indices 0..2000 have residue 2 at n=5 mod 6, else 1')
    print('AUXILIARY_IDENTITIES: E^2, A^3=1, nonlinear reduction, F, F^2=0, E^2*F, correction at 0..2000')
    print('NONVANISHING: 2001/2001 labelled coefficients are nonzero mod 3')
    parity_bad = [n for n in range(4,126) if (bfile[n] % 2) != (1 if n % 4 in (0,1) else 0)]
    assert not parity_bad
    print('CONJECTURE_1_FINITE_CHECK: published range 4..125 agrees; no all-index parity proof claimed')
    print('BFILE_SHA256:', hashlib.sha256(data).hexdigest())
    print('EVIDENCE_LIMIT: finite exact checks; no Lean proof of claim or negation')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
